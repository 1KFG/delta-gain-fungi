#!/usr/bin/env python3
"""Merge the per-chunk BFD placement calls and compare them with the old upo_sites
nearest-neighbour calls. Standard library only (zstd via the `zstd` binary).

  compare_bfd.py --placed bfd_place/chunk_*.tsv --old upo_sites_bfd_all.tsv.zst \
                 --out-merged bfd_place.tsv.zst --out-report compare.md

The merged table is the old upo_sites row (all columns, including the HspUPO-numbered residues
and lineage) with the placement columns appended (`place_*`). Rows are joined on the protein ID
(first whitespace token). A protein present in one input but not the other is counted and
reported, never silently dropped.

The report gives: old clade x confidence tier against the placement clade; placement clade
counts; Figure 1 clade counts; Glu at position 180 by placement clade and, for old Clade I calls,
by where placement puts them; and the share of old `extrapolated` calls that placement confirms.
"""
from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from collections import Counter, defaultdict


def open_text(path, mode="r"):
    if path.endswith(".zst"):
        if mode == "r":
            p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)
            return io.TextIOWrapper(p.stdout)
        p = subprocess.Popen(["zstd", "-q", "-f", "-o", path], stdin=subprocess.PIPE)
        return _Closer(io.TextIOWrapper(p.stdin), p)
    return open(path, mode, newline="")


class _Closer:
    def __init__(self, fh, proc):
        self.fh, self.proc = fh, proc

    def write(self, s):
        return self.fh.write(s)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.fh.close()
        self.proc.wait()


def read_placed(paths):
    out = {}
    for p in paths:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[r["query"].split()[0]] = r
    return out


def tier(conf: str) -> str:
    for t in ("validated_range", "extrapolated", "below_identity_floor"):
        if conf.startswith(t):
            return t
    return conf or "NA"


def old_clade(c: str) -> str:
    return "unassigned" if c in ("", "unassigned") else f"Clade {c}"


def table(rows, title, cols=None):
    """rows: Counter keyed by (row, col) -> markdown crosstab."""
    rkeys = sorted({k[0] for k in rows})
    ckeys = cols or sorted({k[1] for k in rows})
    out = [f"### {title}", "", "| | " + " | ".join(ckeys) + " | total |",
           "|---" * (len(ckeys) + 2) + "|"]
    for r in rkeys:
        vals = [rows.get((r, c), 0) for c in ckeys]
        out.append(f"| {r} | " + " | ".join(str(v) for v in vals) + f" | {sum(vals)} |")
    tot = [sum(rows.get((r, c), 0) for r in rkeys) for c in ckeys]
    out.append("| total | " + " | ".join(str(v) for v in tot) + f" | {sum(tot)} |")
    return "\n".join(out) + "\n"


def compare(placed, old_path, merged_path):
    stats = defaultdict(Counter)
    seen = set()
    pcols = ["our_clade", "our_lwr", "our_stem_lwr", "fig1_clade", "fig1_lwr", "best_edge_lwr",
             "pendant_length", "flags"]
    missing_place = 0
    with open_text(old_path) as fin, open_text(merged_path, "w") as fout:
        lines = (l for l in fin if not l.startswith("#"))
        rd = csv.DictReader(lines, delimiter="\t")
        cols = rd.fieldnames + [f"place_{c}" for c in pcols]
        w = csv.DictWriter(fout, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rd:
            pid = r["protein_id"].split()[0]
            p = placed.get(pid)
            if p is None:
                missing_place += 1
                continue
            seen.add(pid)
            for c in pcols:
                r[f"place_{c}"] = p[c]
            w.writerow(r)
            oc, t, nc, fc = old_clade(r["clade"]), tier(r["clade_confidence"]), p["our_clade"], p["fig1_clade"]
            glu = "E" if r.get("p180_Acid") == "E" else "notE"
            stats["old_vs_new"][(f"{oc} / {t}", nc)] += 1
            stats["new"][(nc, "n")] += 1
            stats["fig1"][(fc, "n")] += 1
            stats["glu_by_new"][(nc, glu)] += 1
            stats["glu_by_fig1"][(fc, glu)] += 1
            if oc == "Clade I":
                stats["oldI_fig1"][(fc, glu)] += 1
            if p["flags"]:
                stats["flags"][(nc, p["flags"])] += 1
            stats["lineage_new"][(r.get("PHYLUM", "") or "NA", nc)] += 1
    return stats, missing_place, len(set(placed) - seen)


def report(stats, missing_place, missing_old, n_placed):
    s = ["# BFD clade calls: placement vs nearest neighbour", "",
         f"Placed proteins: {n_placed}. In old table but not placed: {missing_place}. "
         f"Placed but not in old table: {missing_old}.", ""]
    s.append(table(stats["old_vs_new"], "Old call (clade / tier) x placement clade"))
    s.append(table(stats["new"], "Placement clade (Clades I-V)"))
    s.append(table(stats["fig1"], "Placement clade (Zamocky Fig 1)"))
    s.append(table(stats["glu_by_new"], "Residue 180 (HspUPO) by placement clade",
                   ["E", "notE"]))
    s.append(table(stats["glu_by_fig1"], "Residue 180 by Fig 1 clade", ["E", "notE"]))
    s.append(table(stats["oldI_fig1"], "Old Clade I calls: where placement puts them, and Glu180",
                   ["E", "notE"]))
    s.append(table(stats["flags"], "Flags by placement clade"))
    s.append(table(stats["lineage_new"], "Phylum x placement clade"))
    return "\n".join(s)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--placed", nargs="+", required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--out-merged", required=True)
    ap.add_argument("--out-report", required=True)
    a = ap.parse_args(argv)
    placed = read_placed(a.placed)
    stats, mp, mo = compare(placed, a.old, a.out_merged)
    open(a.out_report, "w").write(report(stats, mp, mo, len(placed)))
    print(f"placed {len(placed)}; old-not-placed {mp}; placed-not-old {mo}", file=sys.stderr)


if __name__ == "__main__":
    main()
