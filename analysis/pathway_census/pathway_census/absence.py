"""Absence check: is a step missing from a genome's annotation also missing from its DNA?

Spec: section 2 anti-goal "Annotation gap read as biology", section 5.4 step 6, section 5.6.
A zero in the copy matrix can mean "gene not annotated". `check_absence` searches the genome
nucleotide sequence with the step's reference proteins (tblastn) for the (genome, step) pairs that
`select_absences` picks.

Outcomes of `check_absence`
---------------------------
present_unannotated  a reference has a tblastn hit with bitscore >= min_bitscore whose HSPs cover
                     >= min_query_cover of the query (union of query intervals on one subject).
fragmented           hits reach min_bitscore for some reference, but no reference reaches
                     min_query_cover. Best hit(s) = all HSPs of one (query, subject) pair with
                     bitscore >= min_bitscore.
absent               no HSP reaches min_bitscore.
unconfirmed          genome_fna is None: no BFD genome FASTA (see genomes.downloaded). This is not
                     evidence of absence, and it is not evidence of presence either.
Out of scope: a hit that covers the query but has a stop codon in frame (pseudogene). MEASURED
2026-11-20 on tests/data/mini_genome_stop.fna (refA, TAA inserted after residue 75): tblastn reports one
HSP over query 1-150 with bitscore 310, so the result is present_unannotated. (One synthetic case only.
Behaviour on real pseudogenes is UNVERIFIED.) A pseudogene call needs a separate gene-structure check.

PROPOSED thresholds (no measured data behind them; set after the first real run):
  min_query_cover = 0.7   min_bitscore = 50 (tblastn bitscore; chance hits are far below this for
                          the 150-aa test queries, MEASURED only on the synthetic fixtures).
`min_bitscore` has no default. `min_query_cover` defaults to the PROPOSED 0.7.
`evalue` (tblastn -evalue, default 1e-3) is PROPOSED, with no measured basis. HSPs must pass both the
e-value and min_bitscore. A nonexistent genome_fna raises FileNotFoundError (never "absent").

Selection (`select_absences`) uses relatives supplied by the caller. The species tree is a separate
plan. For the pilot use `same_genus_relatives`. Genomes with has_proteins False are NOT selected:
their zero in the matrix is not evidence (no proteome was searched). `unevaluated_genomes` lists them.
tblastn input: a plain FASTA path. A .gz/.fna.gz/.fa.gz input is streamed into `scratch_dir`
(default: $SCRATCH, else the system temp dir) and removed after the search.
"""
from __future__ import annotations

import gzip
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

BLAST_MODULE = "ncbi-blast/2.17.0+"
DEFAULT_MIN_QUERY_COVER = 0.7  # PROPOSED
OUTCOMES = ("present_unannotated", "fragmented", "absent", "unconfirmed")


def _blast_cmd(cmd: list[str]) -> list[str]:
    if shutil.which(cmd[0]):
        return cmd
    mod = os.environ.get("PATHWAY_CENSUS_BLAST_MODULE", BLAST_MODULE)
    return ["bash", "-lc", f"module load {mod} >/dev/null 2>&1 && " + shlex.join(cmd)]


def blast_available() -> bool:
    if shutil.which("tblastn"):
        return True
    mod = os.environ.get("PATHWAY_CENSUS_BLAST_MODULE", BLAST_MODULE)
    r = subprocess.run(["bash", "-lc", f"module load {mod} >/dev/null 2>&1; command -v tblastn"],
                       capture_output=True, text=True)
    return r.returncode == 0 and "tblastn" in r.stdout


def select_absences(matrix: pd.DataFrame, relatives: dict[str, list[str]]) -> pd.DataFrame:
    """(locustag, step_id) pairs where the genome has copy_raw == 0, has_proteins is True, and at
    least one listed relative has copy_raw > 0 for the step.

    copy_raw is used: a copy that is only flagged is not an absence, and a relative whose only copy
    is flagged still carries the step. has_proteins comes from the matrix (same source as
    `unevaluated_genomes`). Relatives absent from the matrix count as no copy.
    """
    cols = ["locustag", "step_id"]
    m = matrix[["locustag", "step_id", "copy_raw", "has_proteins"]]
    carriers = m[m.copy_raw > 0].groupby("step_id")["locustag"].agg(set).to_dict()
    cand = m[(m.copy_raw == 0) & (m.has_proteins.astype(bool))]
    keep = [
        any(r in carriers.get(sid, ()) for r in relatives.get(lt, ()))
        for lt, sid in zip(cand.locustag, cand.step_id)
    ]
    return cand.loc[keep, cols].sort_values(cols).reset_index(drop=True)


def unevaluated_genomes(matrix: pd.DataFrame) -> list[str]:
    """Sorted locustags with has_proteins False in the matrix (the column `select_absences` uses).
    Their zero copies are not evidence of absence, so selection skips them. Report them separately."""
    g = matrix.loc[~matrix.has_proteins.astype(bool), "locustag"]
    return sorted(set(g))


def chunk_id(explicit: str | None, locustags_path: Path | str | None) -> str:
    """Chunk id for output names. Explicit id (for example SLURM_ARRAY_TASK_ID) wins; else the
    basename of the --locustags file minus its last suffix (chunk.001.txt -> chunk_001); else "". Characters outside [A-Za-z0-9_-] become "_"."""
    import re
    raw = explicit or (Path(locustags_path).stem if locustags_path else "")
    return re.sub(r"[^A-Za-z0-9_-]", "_", raw)


def chunk_stem(stem: str, chunk: str) -> str:
    """"carotenoid_absence" + "7" -> "carotenoid_absence.chunk7"; empty chunk leaves the stem."""
    return f"{stem}.chunk{chunk}" if chunk else stem


def same_genus_relatives(genomes: pd.DataFrame) -> dict[str, list[str]]:
    """{locustag: [other locustags of the same genus]}. Only genus_source != "unresolved" forms a
    group. A genome is never its own relative. Every genome has a key (empty list if none)."""
    out: dict[str, list[str]] = {lt: [] for lt in genomes.locustag}
    ok = genomes[genomes.genus_source != "unresolved"]
    for _, grp in ok.groupby("genus"):
        members = sorted(grp.locustag)
        for lt in members:
            out[lt] = [x for x in members if x != lt]
    return out


def _covered(intervals: list[tuple[int, int]]) -> int:
    """Number of query positions covered by the union of 1-based inclusive intervals."""
    merged: list[list[int]] = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s + 1 for s, e in merged)


def check_absence(genome_fna: Path | None, step_refs: Path, min_bitscore: float,
                  min_query_cover: float = DEFAULT_MIN_QUERY_COVER,
                  scratch_dir: Path | None = None, cpus: int = 1, evalue: float = 1e-3) -> str:
    """Search genome_fna with the step reference proteins (tblastn). See module docstring."""
    if not 0.0 < min_query_cover <= 1.0:
        raise ValueError(f"min_query_cover must be in (0, 1], got {min_query_cover}")
    if genome_fna is None:
        return "unconfirmed"
    genome_fna, step_refs = Path(genome_fna), Path(step_refs)
    if not genome_fna.is_file():
        raise FileNotFoundError(f"genome FASTA not found: {genome_fna}")
    base = scratch_dir or os.environ.get("SCRATCH") or None
    with tempfile.TemporaryDirectory(dir=base) as td:
        subject = genome_fna
        if genome_fna.suffix == ".gz":
            subject = Path(td) / "subject.fna"
            with gzip.open(genome_fna, "rb") as src, open(subject, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
        out = Path(td) / "hits.tsv"
        cmd = ["tblastn", "-query", str(step_refs), "-subject", str(subject),
               "-evalue", str(evalue), "-num_threads", str(cpus), "-out", str(out),
               "-outfmt", "6 qseqid sseqid bitscore qstart qend qlen"]
        r = subprocess.run(_blast_cmd(cmd), capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"tblastn failed: {r.stderr.strip()[:500]}")
        rows = [ln.split("\t") for ln in out.read_text().splitlines() if ln.strip()]
    best_cover = None
    pairs: dict[tuple[str, str], list[tuple[int, int]]] = {}
    qlen: dict[str, int] = {}
    for q, s, bits, qs, qe, ql in rows:
        if float(bits) < min_bitscore:
            continue
        a, b = sorted((int(qs), int(qe)))
        pairs.setdefault((q, s), []).append((a, b))
        qlen[q] = int(ql)
    for (q, _s), iv in pairs.items():
        cov = _covered(iv) / qlen[q]
        best_cover = cov if best_cover is None else max(best_cover, cov)
    if best_cover is None:
        return "absent"
    return "present_unannotated" if best_cover >= min_query_cover else "fragmented"
