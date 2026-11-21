"""Candidate retrieval: step HMM build, hmmsearch, domtbl parsing (spec section 5.4, steps 1-3).

Tools are run through `run_tool`. If a tool is not on PATH, the command runs under
`bash -lc 'module load <modules> && ...'`. Set PATHWAY_CENSUS_MODULES to change the module list.
Inside a SLURM job, load the modules in the batch script; `run_tool` then uses PATH directly.
"""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

COLUMNS = ["protein_id", "hmm_name", "score", "evalue", "ali_from", "ali_to", "hmm_from", "hmm_to"]
DEFAULT_MODULES = "hmmer/3.4 mafft/7.505 mmseqs2/17-b804f"


class UnverifiedPathwayError(RuntimeError):
    """Raised when a pathway with verified: false is used without --provisional."""


def run_tool(cmd: list[str], stdout_path: Path | None = None) -> None:
    """Run cmd. Raise CalledProcessError on failure. Wrap in `module load` if the tool is not on PATH."""
    if shutil.which(cmd[0]) is None:
        mods = os.environ.get("PATHWAY_CENSUS_MODULES", DEFAULT_MODULES)
        line = f"module load {mods} >/dev/null 2>&1 && " + shlex.join(cmd)
        cmd = ["bash", "-lc", line]
    if stdout_path is None:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    else:
        with open(stdout_path, "wb") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.PIPE)


def tool_available(name: str) -> bool:
    """True if `name` is on PATH or becomes available after `module load`."""
    if shutil.which(name):
        return True
    mods = os.environ.get("PATHWAY_CENSUS_MODULES", DEFAULT_MODULES)
    r = subprocess.run(["bash", "-lc", f"module load {mods} >/dev/null 2>&1; command -v {shlex.quote(name)}"],
                       capture_output=True, text=True)
    return r.returncode == 0


def provisional_tag(pathway, provisional: bool) -> str:
    """Return "PROVISIONAL" or "". Refuse an unverified pathway unless provisional is True."""
    if pathway.verified:
        return ""
    if not provisional:
        raise UnverifiedPathwayError(
            f"pathway '{pathway.id}' has verified: false (not reviewed by the user, spec 5.8). "
            "Pass --provisional to proceed; all outputs are then marked PROVISIONAL.")
    return "PROVISIONAL"


def tagged_name(stem: str, ext: str, tag: str) -> str:
    return f"{stem}.{tag}.{ext}" if tag else f"{stem}.{ext}"


def parse_domtbl(path: Path) -> pd.DataFrame:
    """Parse `hmmsearch --domtblout`. score = domain (best-1-domain-style) score, evalue = i-Evalue.

    Columns of the file (0-based): 0 target, 3 query, 12 i-Evalue, 13 domain score,
    15-16 hmm from/to, 17-18 ali from/to. For hmmsearch the target is the protein, the query is the HMM.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            rows.append((f[0], f[3], float(f[13]), float(f[12]),
                         int(f[17]), int(f[18]), int(f[15]), int(f[16])))
    return pd.DataFrame(rows, columns=COLUMNS)


def build_hmm(ref_fasta: Path, out_hmm: Path) -> None:
    """MAFFT alignment of the reference proteins, then hmmbuild. One reference gives a single-sequence HMM."""
    ref_fasta, out_hmm = Path(ref_fasta), Path(out_hmm)
    name = out_hmm.name.split(".")[0]
    with tempfile.TemporaryDirectory(dir=os.environ.get("SCRATCH")) as td:
        aln = Path(td) / "aln.faa"
        run_tool(["mafft", "--auto", "--quiet", str(ref_fasta)], stdout_path=aln)
        run_tool(["hmmbuild", "-n", name, "--amino", str(out_hmm), str(aln)])


def search_step(hmm: Path, proteins: Path, out_domtbl: Path, cpus: int) -> None:
    """hmmsearch --domtblout with the report threshold -E 10 (loose; the score cutoff is applied later).

    A .gz protein file is streamed to hmmsearch through zcat.
    """
    hmm, proteins, out_domtbl = Path(hmm), Path(proteins), Path(out_domtbl)
    base = ["hmmsearch", "--cpu", str(cpus), "--noali", "-o", os.devnull,
            "--domtblout", str(out_domtbl), "-E", "10", str(hmm)]
    if proteins.suffix == ".gz":
        line = f"zcat {shlex.quote(str(proteins))} | " + shlex.join(base + ["-"])
        if shutil.which("hmmsearch") is None:
            mods = os.environ.get("PATHWAY_CENSUS_MODULES", DEFAULT_MODULES)
            line = f"module load {mods} >/dev/null 2>&1 && " + line
        subprocess.run(["bash", "-lc", "set -o pipefail; " + line], check=True, capture_output=True, text=True)
    else:
        run_tool(base + [str(proteins)])


def candidates_for_step(domtbl: pd.DataFrame, min_score: float) -> pd.DataFrame:
    """One row per protein: its best-scoring domain, kept if score >= min_score. min_score has no default."""
    if domtbl.empty:
        return pd.DataFrame(columns=COLUMNS)
    best = domtbl.sort_values(["score", "evalue"], ascending=[False, True]).drop_duplicates("protein_id")
    best = best[best.score >= min_score]
    return best.sort_values("score", ascending=False).reset_index(drop=True)[COLUMNS]


def dereplicate(fasta_in: Path, fasta_out: Path, identity: float, coverage: float = 0.8) -> None:
    """Cluster with `mmseqs easy-cluster` and write one representative per cluster.

    identity: --min-seq-id (0-1). coverage: -c, bidirectional (--cov-mode 0).
    """
    if not 0.0 < identity <= 1.0:
        raise ValueError(f"identity must be in (0, 1], got {identity}")
    with tempfile.TemporaryDirectory(dir=os.environ.get("SCRATCH")) as td:
        prefix = Path(td) / "clu"
        run_tool(["mmseqs", "easy-cluster", str(fasta_in), str(prefix), str(Path(td) / "tmp"),
                  "--min-seq-id", str(identity), "-c", str(coverage), "--cov-mode", "0", "-v", "1"])
        shutil.copyfile(f"{prefix}_rep_seq.fasta", fasta_out)
