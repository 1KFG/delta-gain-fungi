#!/usr/bin/env python3.12
"""Paired deposited-versus-BFD comparison (thin wrapper over pathway_census.paired).

Inputs: deposited copy matrix, BFD copy matrix, genome table (each plain, .gz or .zst), pathway YAML.
The genome table needs: locustag, prior_status, has_proteins, has_deposited_proteins,
annotation_provider, busco_complete. Both matrices are keyed by the same locustag.
Outputs in --out-dir (names carry .PROVISIONAL when the pathway is verified: false and --provisional is
given): <pathway>_paired_table.tsv.zst, _paired_summary.tsv, _paired_dropped_providers.tsv,
_paired_by_busco.tsv, _paired_unpaired.json.
Refuses a verified: false pathway unless --provisional.
"""
import argparse
import gzip
import io
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402
from pathway_census.matrix import parse_bool_col  # noqa: E402
from pathway_census.paired import (DEFAULT_MIN_PAIRS, DEFAULT_N_BOOT, DEFAULT_SEED, paired_table,  # noqa: E402
                                   summarize_paired, summarize_paired_by_busco)


def read_tsv(path: Path, **kw) -> pd.DataFrame:
    p = str(path)
    if p.endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout
        return pd.read_csv(io.BytesIO(raw), sep="\t", comment="#", low_memory=False, **kw)
    if p.endswith(".gz"):
        with gzip.open(p, "rb") as fh:
            return pd.read_csv(io.BytesIO(fh.read()), sep="\t", comment="#", low_memory=False, **kw)
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False, **kw)


def read_matrix(path: Path) -> pd.DataFrame:
    m = read_tsv(path)
    m["has_proteins"] = parse_bool_col(m.has_proteins.astype(object).where(m.has_proteins.notna(), None),
                                       "has_proteins", allow_na=False).astype(bool)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--deposited-matrix", type=Path, required=True)
    ap.add_argument("--bfd-matrix", type=Path, required=True)
    ap.add_argument("--genomes", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("results"))
    ap.add_argument("--min-pairs", type=int, default=DEFAULT_MIN_PAIRS, help="PROPOSED default")
    ap.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT, help="PROPOSED default")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED, help="PROPOSED default")
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")
    gen = read_tsv(a.genomes, dtype={"has_proteins": str, "has_deposited_proteins": str},
                   keep_default_na=False, na_values=[""])
    gen["has_proteins"] = parse_bool_col(gen.has_proteins, "has_proteins", allow_na=False).astype(bool)
    gen["has_deposited_proteins"] = parse_bool_col(gen.has_deposited_proteins, "has_deposited_proteins",
                                                   allow_na=True)
    dep = read_matrix(a.deposited_matrix)
    bfd = read_matrix(a.bfd_matrix)
    t, n_unpaired = paired_table(dep, bfd, gen, return_unpaired=True)
    summ, dropped = summarize_paired(t, a.min_pairs, a.n_boot, a.seed)
    bus = summarize_paired_by_busco(t)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    nm = lambda stem, ext: a.out_dir / hits.tagged_name(f"{pw.id}_{stem}", ext, tag)
    tpath = nm("paired_table", "tsv")
    t.to_csv(tpath, sep="\t", index=False)
    subprocess.run(["zstd", "-q", "-f", "-T0", "--rm", str(tpath), "-o", str(tpath) + ".zst"], check=True)
    note = f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n" if tag else ""
    params = f"# min_pairs={a.min_pairs} n_boot={a.n_boot} seed={a.seed} bootstrap=genomes ci=2.5-97.5\n"
    for stem, df, hdr in (("paired_summary", summ, params), ("paired_dropped_providers", dropped, ""),
                          ("paired_by_busco", bus, "")):
        with open(nm(stem, "tsv"), "w") as fh:
            fh.write(note + hdr)
            df.to_csv(fh, sep="\t", index=False)
    n_pairs = int(t.locustag.nunique())
    nm("paired_unpaired", "json").write_text(json.dumps(
        dict(n_pairs=n_pairs, n_unpaired=n_unpaired, provisional=bool(tag)), indent=2))
    print(f"{n_pairs} pairs; unpaired {n_unpaired}; outputs in {a.out_dir}")


if __name__ == "__main__":
    main()
