#!/usr/bin/env python3.12
"""Absence check (spec 5.4 step 6): tblastn each selected (genome, step) absence against the genome DNA.

Selection: genome has copy_raw == 0 for a step, has_proteins True, and a same-genus relative has
copy_raw > 0 (pilot relatives = same_genus_relatives). Genomes with has_proteins False are not
selected; they are written to <pathway>_unevaluated_genomes.tsv.
A genome with downloaded == False has no BFD FASTA. It gets outcome "unconfirmed", reason
"no_bfd_fasta". The CLI does not look for a FASTA anywhere else.
Genome FASTA path = genomes.scaffolds_path(proteins path from seqdb/genomes.tsv).
Inputs (plain, .gz or .zst TSV): --matrix (build_matrix output), --genomes (genome table).
Output: <out-dir>/<pathway>_absence[.chunk<ID>].tsv (columns locustag, step_id, outcome, reason); names carry
.PROVISIONAL, and the file starts with a "# PROVISIONAL" line, when the pathway is verified: false and
--provisional is given. Refuses a verified: false pathway without --provisional.
--min-bitscore has no default. --min-query-cover default 0.7 is PROPOSED.
Chunk ID: --chunk-id, else the --locustags file basename stem; it keeps chunk jobs from overwriting each other.
Optional --locustags FILE restricts the run to listed genomes (for chunking into SLURM jobs).
"""
import argparse
import io
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.absence import (DEFAULT_MIN_QUERY_COVER, check_absence, chunk_id,  # noqa: E402
                                    chunk_stem, same_genus_relatives, select_absences,
                                    unevaluated_genomes)
from pathway_census.defs import load_pathway  # noqa: E402
from pathway_census.genomes import scaffolds_path  # noqa: E402
from pathway_census.matrix import parse_bool_col  # noqa: E402


def read_tsv(path: Path, **kw) -> pd.DataFrame:
    """Read a TSV (plain or .zst). Only LEADING lines starting with "#" are skipped (stamps);
    "#" inside a field is kept."""
    if str(path).endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout
    else:
        raw = Path(path).read_bytes()
    lines = raw.splitlines(keepends=True)
    n = 0
    while n < len(lines) and lines[n].startswith(b"#"):
        n += 1
    return pd.read_csv(io.BytesIO(b"".join(lines[n:])), sep="\t", low_memory=False, **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--genomes", type=Path, required=True)
    ap.add_argument("--seqdb-genomes", type=Path, required=True, help="seqdb/genomes.tsv (LOCUSTAG, out, proteins, cds)")
    ap.add_argument("--refs-dir", type=Path, required=True, help="directory with <step_id>.faa")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-bitscore", type=float, required=True)
    ap.add_argument("--min-query-cover", type=float, default=DEFAULT_MIN_QUERY_COVER)
    ap.add_argument("--locustags", type=Path, default=None)
    ap.add_argument("--chunk-id", default="", help="suffix for output names (sbatch passes SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--evalue", type=float, default=1e-3, help="tblastn -evalue (PROPOSED)")
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")

    m = read_tsv(a.matrix)
    m["has_proteins"] = parse_bool_col(m.has_proteins, "has_proteins", allow_na=False).astype(bool)
    gen = read_tsv(a.genomes)
    gen["has_proteins"] = parse_bool_col(gen.has_proteins, "has_proteins", allow_na=False).astype(bool)
    gen["downloaded"] = parse_bool_col(gen.downloaded, "downloaded", allow_na=False).astype(bool)
    seqdb = pd.read_csv(a.seqdb_genomes, sep="\t", header=None, usecols=[0, 2], names=["locustag", "proteins"])
    fasta = {lt: scaffolds_path(p) for lt, p in zip(seqdb.locustag, seqdb.proteins)}
    dl = dict(zip(gen.locustag, gen.downloaded))

    sel = select_absences(m, same_genus_relatives(gen))
    if a.locustags:
        keep = set(a.locustags.read_text().split())
        sel = sel[sel.locustag.isin(keep)]
    refs = {s.id: a.refs_dir / f"{s.id}.faa" for s in pw.steps}
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for lt, sid in zip(sel.locustag, sel.step_id):
        path = fasta.get(lt) if dl.get(lt) else None
        if path is None:
            rows.append((lt, sid, "unconfirmed", "no_bfd_fasta"))
            continue
        outcome = check_absence(Path(path), refs[sid], a.min_bitscore, a.min_query_cover, cpus=a.cpus, evalue=a.evalue)
        rows.append((lt, sid, outcome, ""))
    res = pd.DataFrame(rows, columns=["locustag", "step_id", "outcome", "reason"])
    chunk = chunk_id(a.chunk_id, a.locustags)
    out = a.out_dir / hits.tagged_name(chunk_stem(f"{pw.id}_absence", chunk), "tsv", tag)
    with open(out, "w") as fh:
        if tag:
            fh.write(f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n")
        fh.write(f"# min_bitscore={a.min_bitscore} min_query_cover={a.min_query_cover}\n")
        res.to_csv(fh, sep="\t", index=False)
    uneval = a.out_dir / hits.tagged_name(chunk_stem(f"{pw.id}_unevaluated_genomes", chunk), "tsv", tag)
    with open(uneval, "w") as fh:
        if tag:
            fh.write(f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n")
        fh.write("locustag\n" + "".join(f"{x}\n" for x in unevaluated_genomes(m)))
    print(f"{len(res)} pairs -> {out}; outcomes {res.outcome.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
