#!/usr/bin/env python3.12
"""Build the genome x step copy matrix and the phylum summary.

Inputs: assignment TSV, flags TSV, genome table TSV (each plain, .gz or .zst), pathway YAML.
Outputs: <out-dir>/<pathway>_matrix.tsv.zst and <out-dir>/<pathway>_phylum_summary.tsv
(names carry .PROVISIONAL when the pathway is verified: false and --provisional is given).
Refuses a verified: false pathway unless --provisional.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402
from pathway_census.matrix import DEFAULT_REQUIRE, build_matrix, phylum_summary  # noqa: E402

BOOL = {"True": True, "False": False}


def read_tsv(path: Path) -> pd.DataFrame:
    if str(path).endswith(".zst"):
        out = subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout
        import io
        return pd.read_csv(io.BytesIO(out), sep="\t", comment="#", low_memory=False)
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--assigned", type=Path, required=True)
    ap.add_argument("--flags", type=Path, required=True)
    ap.add_argument("--genomes", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("results"))
    ap.add_argument("--require", nargs="+", default=list(DEFAULT_REQUIRE))
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")
    asg = read_tsv(a.assigned)
    fl = read_tsv(a.flags)
    for c in ("domain_ok", "motif_ok", "truncated", "contam_suspect", "hap_dup_suspect"):
        if c in fl.columns:
            fl[c] = fl[c].map(lambda v: BOOL.get(str(v), pd.NA)).astype("boolean")  # empty/NA stays NA
    gen = read_tsv(a.genomes)
    gen["has_proteins"] = gen.has_proteins.map(lambda v: BOOL.get(str(v), v)).astype(bool)
    m = build_matrix(asg, fl, gen, pw, require=tuple(a.require))
    a.out_dir.mkdir(parents=True, exist_ok=True)
    mpath = a.out_dir / hits.tagged_name(f"{pw.id}_matrix", "tsv", tag)
    m.to_csv(mpath, sep="\t", index=False)
    subprocess.run(["zstd", "-q", "-f", "-T0", "--rm", str(mpath), "-o", str(mpath) + ".zst"], check=True)
    spath = a.out_dir / hits.tagged_name(f"{pw.id}_phylum_summary", "tsv", tag)
    with open(spath, "w") as fh:
        if tag:
            fh.write(f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n")
        fh.write(f"# count_col=copy_filtered require={','.join(a.require)}\n")
        phylum_summary(m, gen).to_csv(fh, sep="\t", index=False)
    print(f"{len(m)} rows -> {mpath}.zst; summary -> {spath}")


if __name__ == "__main__":
    main()
