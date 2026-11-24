#!/usr/bin/env python3
"""Sweep score cutoffs for one step over a low-floor search of a set of genomes.

Input: candidates.chunk*.tsv.zst from bin/search_chunk.py run with a low floor, the searched-genome manifest of the same run,
and the taxonomy table. For each cutoff: the share of searched genomes with at least one call per class, and the distribution
of calls per genome. Also the median i-Evalue of the hits in each score band (the i-Evalue depends on the size of the searched
set, so a bit-score cutoff is the portable one; the E-value column is for reading only).

Writes sweep_by_class.tsv, sweep_copies.tsv and score_bands.tsv into --out-dir.
"""
from __future__ import annotations

import argparse
import csv
import glob
import io
import statistics
import subprocess
from pathlib import Path


def read_tsv(path):
    path = str(path)
    if path.endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    else:
        raw = Path(path).read_text()
    return list(csv.DictReader(io.StringIO("\n".join(l for l in raw.splitlines() if not l.startswith("#"))), delimiter="\t"))


def load_hits(search_dir, step):
    out = []
    for f in sorted(glob.glob(str(Path(search_dir) / "candidates.chunk*.tsv.zst"))):
        out += [(r["protein_id"].split("_")[0], float(r["score"]), float(r["evalue"])) for r in read_tsv(f) if r["step_id"] == step]
    return out


def searched_genomes(search_dir):
    g = []
    for f in sorted(glob.glob(str(Path(search_dir) / "searched.chunk*.tsv"))):
        g += [r["locustag"] for r in read_tsv(f) if r["status"] == "searched"]
    return g


def by_class_table(hits, genomes, classes, cutoffs, min_n=1):
    """Return rows (class, n_genomes, {cutoff: fraction with >=1 hit at or above}) for classes with at least min_n genomes."""
    top = {}
    for g, sc, _ in hits:
        if sc > top.get(g, -1):
            top[g] = sc
    members = {}
    for g in genomes:
        members.setdefault(classes.get(g) or "(no label)", []).append(g)
    rows = []
    for cl, gs in members.items():
        if len(gs) < min_n:
            continue
        rows.append((cl, len(gs), {c: sum(top.get(g, -1) >= c for g in gs) / len(gs) for c in cutoffs}))
    rows.sort(key=lambda r: -r[1])
    return rows


def copies_table(hits, genomes, cutoffs):
    """cutoff -> {0: n genomes with 0 calls, 1: ..., 2: ..., '3+': ...}."""
    out = {}
    for c in cutoffs:
        cnt = {}
        for g, sc, _ in hits:
            if sc >= c:
                cnt[g] = cnt.get(g, 0) + 1
        d = {0: 0, 1: 0, 2: 0, "3+": 0}
        for g in genomes:
            n = cnt.get(g, 0)
            d[n if n < 3 else "3+"] += 1
        out[c] = d
    return out


def score_bands(hits, edges):
    """Rows (low, high, n hits, median i-Evalue) for score bands [low, high)."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ev = [e for _, s, e in hits if lo <= s < hi]
        rows.append((lo, hi, len(ev), statistics.median(ev) if ev else None))
    return rows


def main(argv=None):
    import duckdb
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--search-dir", required=True, type=Path)
    ap.add_argument("--step", required=True)
    ap.add_argument("--species", required=True, type=Path)
    ap.add_argument("--cutoffs", required=True, help="comma-separated bit scores")
    ap.add_argument("--min-n", type=int, default=15)
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args(argv)
    cutoffs = [float(x) for x in a.cutoffs.split(",")]
    hits = load_hits(a.search_dir, a.step)
    genomes = searched_genomes(a.search_dir)
    classes = {r[0]: r[1] for r in duckdb.connect().execute(f"select LOCUSTAG, \"CLASS\" from read_parquet('{a.species}')").fetchall()}
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rows = by_class_table(hits, genomes, classes, cutoffs, a.min_n)
    with open(a.out_dir / "sweep_by_class.tsv", "w") as fh:
        fh.write("class\tn_genomes\t" + "\t".join(f"frac_ge_{c:g}" for c in cutoffs) + "\n")
        tot = [sum(1 for g in genomes)]
        for cl, n, d in rows:
            fh.write(f"{cl}\t{n}\t" + "\t".join(f"{d[c]:.3f}" for c in cutoffs) + "\n")
        allfrac = {c: len({g for g, s, _ in hits if s >= c}) / len(genomes) for c in cutoffs}
        fh.write(f"ALL\t{len(genomes)}\t" + "\t".join(f"{allfrac[c]:.3f}" for c in cutoffs) + "\n")
    cp = copies_table(hits, genomes, cutoffs)
    with open(a.out_dir / "sweep_copies.tsv", "w") as fh:
        fh.write("cutoff\tgenomes_0\tgenomes_1\tgenomes_2\tgenomes_3plus\n")
        for c in cutoffs:
            d = cp[c]
            fh.write(f"{c:g}\t{d[0]}\t{d[1]}\t{d[2]}\t{d['3+']}\n")
    edges = sorted(set([0.0, 100.0, 150.0, 200.0, 250.0, 300.0, 369.3, 450.0, 600.0, 10000.0] + cutoffs))
    with open(a.out_dir / "score_bands.tsv", "w") as fh:
        fh.write("low\thigh\tn_hits\tmedian_i_evalue\n")
        for lo, hi, n, ev in score_bands(hits, edges):
            fh.write(f"{lo:g}\t{hi:g}\t{n}\t{'' if ev is None else f'{ev:.2g}'}\n")
    print(f"{len(hits)} hits in {len(genomes)} genomes; wrote {a.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
