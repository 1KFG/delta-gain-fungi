#!/usr/bin/env python3
"""Share of searched genomes with at least one search hit, per taxonomic class, at several score floors, for one step.

A hit is a row of the candidate files (bin/search_chunk.py) with step_id == --step and score >= floor. Genomes are those of
results/genomes.tsv.zst with has_proteins true (the logic of reports/carotenoid/figures/bin/fig_presence_heatmap.py). Classes with fewer
than --min-n genomes are listed as `(other, small classes)`. Writes <out> (class, n_genomes, one share column per floor, ALL row).
These are search hits, not placed calls.
"""
from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import io
import subprocess
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("pvh", Path(__file__).parent / "placed_vs_hits_by_class.py")
pvh = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pvh)


def hits_at_floors(search_dir, step, floors):
    """floor -> set of locustags with a hit of the step at score >= floor."""
    out = {f: set() for f in floors}
    for p in sorted(glob.glob(str(Path(search_dir) / "candidates.chunk*.tsv.zst"))):
        for r in pvh.read_tsv(p):
            if r["step_id"] != step:
                continue
            sc = float(r["score"])
            g = r["protein_id"].split("_")[0]
            for f in floors:
                if sc >= f:
                    out[f].add(g)
    return out


def class_table(genomes, classes, sets, min_n):
    rows = pvh.by_class(genomes, classes, sets)
    big = [r for r in rows if r[1] >= min_n]
    small = [r for r in rows if r[1] < min_n]
    if small:
        n = sum(r[1] for r in small)
        agg = {k: sum(r[2][k] * r[1] for r in small) / n for k in sets}
        big.append(("(other, small classes)", n, agg))
    return big


def main(argv=None):
    import duckdb
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", required=True)
    ap.add_argument("--search-dir", required=True, type=Path)
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--species", required=True, type=Path)
    ap.add_argument("--floors", required=True, nargs="+", type=float)
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    genomes = [r["locustag"] for r in pvh.read_tsv(a.genomes) if r["has_proteins"].strip().lower() == "true"]
    classes = {r[0]: r[1] for r in duckdb.connect().execute(
        f"select LOCUSTAG, \"CLASS\" from read_parquet('{a.species}')").fetchall()}
    hits = hits_at_floors(a.search_dir, a.step, a.floors)
    sets = {f"ge_{f:g}": hits[f] for f in a.floors}
    rows = class_table(genomes, classes, sets, a.min_n)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("class\tn_genomes\t" + "\t".join(f"share_{k}" for k in sets) + "\n")
        for c, n, d in rows:
            fh.write(f"{c}\t{n}\t" + "\t".join(f"{d[k]:.4f}" for k in sets) + "\n")
        gs = set(genomes)
        fh.write(f"ALL\t{len(genomes)}\t" + "\t".join(f"{len(sets[k] & gs) / len(genomes):.4f}" for k in sets) + "\n")
    print(f"{len(genomes)} genomes -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
