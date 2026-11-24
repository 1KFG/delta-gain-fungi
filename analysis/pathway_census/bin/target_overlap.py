#!/usr/bin/env python3
"""Overlap of the target sets of several step-5 runs (proteins and genomes), and of raw search hits at a score floor.

--run LABEL=DIR (repeatable): DIR is a step-5 output folder. The target set = proteins with clade_call == target in
DIR/placements_detail.tsv. Writes a TSV with, for each protein set size and each pair, the number of shared proteins and genomes,
and for the whole group the number of proteins in k target sets (k = 1..n).
--hits-dir D --hit-steps a,b,c --hit-floor F: also the overlap of proteins with a search hit of each step at score >= F
(all steps of the same search folder). A protein can have a row for several steps in the candidate files.
"""
from __future__ import annotations

import argparse
import csv
import glob
import itertools
import sys
from collections import Counter
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def target_set(run_dir):
    return {r["protein_id"] for r in read_tsv(Path(run_dir) / "placements_detail.tsv") if r["clade_call"] == "target"}


def genomes(proteins):
    return {p.split("_")[0] for p in proteins}


def overlaps(sets):
    """sets: ordered {label: set}. Return (rows, depth) where rows = [(a, b, shared proteins, shared genomes)] for each pair
    and depth = {k: number of proteins that are in exactly k sets}."""
    rows = [(a, b, len(sets[a] & sets[b]), len(genomes(sets[a] & sets[b]))) for a, b in itertools.combinations(sets, 2)]
    depth = Counter(sum(p in s for s in sets.values()) for p in set().union(*sets.values())) if sets else Counter()
    return rows, dict(sorted(depth.items()))


def hit_sets(hits_dir, steps, floor):
    import subprocess
    out = {s: set() for s in steps}
    for p in sorted(glob.glob(str(Path(hits_dir) / "candidates.chunk*.tsv.zst"))):
        txt = subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout.decode().splitlines()
        rows = csv.DictReader((l for l in txt if not l.startswith("#")), delimiter="\t")
        for r in rows:
            if r["step_id"] in out and float(r["score"]) >= floor:
                out[r["step_id"]].add(r["protein_id"])
    return out


def write(path, title, sets, rows, depth):
    with open(path, "a") as fh:
        fh.write(f"# {title}\n")
        fh.write("set\tn_proteins\tn_genomes\n")
        for k, s in sets.items():
            fh.write(f"{k}\t{len(s)}\t{len(genomes(s))}\n")
        fh.write("pair_a\tpair_b\tshared_proteins\tshared_genomes\n")
        for a, b, n, g in rows:
            fh.write(f"{a}\t{b}\t{n}\t{g}\n")
        fh.write("in_k_sets\tn_proteins\n")
        for k, n in depth.items():
            fh.write(f"{k}\t{n}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="append", default=[], metavar="LABEL=DIR")
    ap.add_argument("--hits-dir", type=Path)
    ap.add_argument("--hit-steps", default="")
    ap.add_argument("--hit-floor", type=float, default=200.0)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    a.out.write_text("")
    if a.run:
        sets = {k: target_set(v) for k, v in (x.split("=", 1) for x in a.run)}
        write(a.out, "target sets (clade_call == target)", sets, *overlaps(sets))
    if a.hits_dir and a.hit_steps:
        sets = hit_sets(a.hits_dir, a.hit_steps.split(","), a.hit_floor)
        write(a.out, f"search hits at score >= {a.hit_floor:g}", sets, *overlaps(sets))
    print(a.out.read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
