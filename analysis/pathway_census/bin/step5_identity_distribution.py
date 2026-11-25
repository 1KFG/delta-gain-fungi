#!/usr/bin/env python3.12
"""Identity structure of a set of domain sequences (redundancy check before choosing a clustering identity).

Input: an all-against-all mmseqs2 table (columns query, target, pident in percent; others ignored) and the FASTA of the same
sequences. Output (stdout and --out TSV):
  exact duplicates      sequences, distinct sequences, largest group of identical sequences
  nearest-neighbour     for each sequence, the highest percent identity to any OTHER sequence; counts per bin
  pairs                 counts per bin of all pairs a < b found in the table (pairs the search did not report are not counted)
Bins (percent identity): <50, 50-70, 70-80, 80-90, 90-95, 95-99, 99-<100, 100.
mmseqs2 reports the best alignment per pair, not the global identity of unaligned regions. Plain or .gz/.zst FASTA.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import subprocess
import sys
from collections import Counter

BINS = [(0, 50, "<50"), (50, 70, "50-70"), (70, 80, "70-80"), (80, 90, "80-90"), (90, 95, "90-95"),
        (95, 99, "95-99"), (99, 100, "99-<100"), (100, 101, "100")]


def open_text(path):
    p = str(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"))
    if p.endswith(".zst"):
        return io.StringIO(subprocess.run(["zstd", "-dc", p], check=True, capture_output=True, text=True).stdout)
    return open(p)


def read_fasta(path):
    seqs, cur = {}, None
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def bin_of(pident: float) -> str:
    for lo, hi, name in BINS:
        if lo <= pident < hi:
            return name
    return BINS[-1][2]


def duplicate_summary(seqs: dict) -> dict:
    groups = Counter(seqs.values())
    return {"sequences": len(seqs), "distinct": len(groups), "largest_group": max(groups.values()) if groups else 0}


def nearest_and_pairs(rows, ids):
    """rows: iterable of (query, target, pident). -> (nearest {id: best identity to another id}, pairs {(a,b): best identity}).
    Sequences with no other hit are absent from `nearest`."""
    nearest, pairs = {}, {}
    ids = set(ids)
    for q, t, p in rows:
        if q == t or q not in ids or t not in ids:
            continue
        if p > nearest.get(q, -1):
            nearest[q] = p
        key = (q, t) if q < t else (t, q)
        if p > pairs.get(key, -1):
            pairs[key] = p
    return nearest, pairs


def read_m8(path):
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            yield f[0], f[1], float(f[2])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--m8", required=True, help="all-against-all table: query, target, pident (percent)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seqs = read_fasta(a.fasta)
    dup = duplicate_summary(seqs)
    nearest, pairs = nearest_and_pairs(read_m8(a.m8), seqs)
    nn_bins = Counter(bin_of(v) for v in nearest.values())
    nn_bins["no other hit"] = len(seqs) - len(nearest)
    pair_bins = Counter(bin_of(v) for v in pairs.values())
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["table", "bin", "count"])
        for k, v in dup.items():
            w.writerow(["duplicates", k, v])
        for _, _, name in BINS:
            w.writerow(["nearest_neighbour", name, nn_bins.get(name, 0)])
        w.writerow(["nearest_neighbour", "no other hit", nn_bins["no other hit"]])
        for _, _, name in BINS:
            w.writerow(["pairs", name, pair_bins.get(name, 0)])
    print(dup)
    print("nearest neighbour:", {n: nn_bins.get(n, 0) for _, _, n in BINS}, "no other hit", nn_bins["no other hit"])
    print("pairs:", {n: pair_bins.get(n, 0) for _, _, n in BINS})


if __name__ == "__main__":
    sys.exit(main())
