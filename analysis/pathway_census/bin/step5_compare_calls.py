#!/usr/bin/env python3.12
"""Compare the clade calls of two step 5 runs of the same step (for example the 300-representative run and a v2 run).

Input: the placements_detail.tsv of each run (columns protein_id, clade_call, long_pendant). Output (stdout and --out TSV):
  transitions   old call -> new call, count of proteins (proteins in both runs)
  only_old / only_new   proteins that are in one run only
  direction     changed calls grouped by direction: old_call>new_call, count
  long_pendant  number of proteins with long_pendant True in each run
Nothing here is a call; it counts differences between two placements. Reads plain or .gz/.zst files.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import subprocess
import sys
from collections import Counter


def open_text(path):
    p = str(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"))
    if p.endswith(".zst"):
        return io.StringIO(subprocess.run(["zstd", "-dc", p], check=True, capture_output=True, text=True).stdout)
    return open(p, newline="")


def read_calls(path):
    """-> ({protein_id: clade_call}, n_long_pendant)"""
    calls, lp = {}, 0
    with open_text(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            calls[r["protein_id"]] = r["clade_call"]
            lp += str(r.get("long_pendant", "")).lower() == "true"
    return calls, lp


def compare(old, new):
    """old, new: {protein_id: call}. -> dict with transitions Counter, only_old, only_new, changed Counter."""
    both = old.keys() & new.keys()
    trans = Counter((old[p], new[p]) for p in both)
    changed = Counter({k: v for k, v in trans.items() if k[0] != k[1]})
    return {"transitions": trans, "changed": changed, "n_both": len(both),
            "only_old": len(old.keys() - new.keys()), "only_new": len(new.keys() - old.keys())}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--out", required=True, help="TSV of transitions")
    a = ap.parse_args()
    old, lp_old = read_calls(a.old)
    new, lp_new = read_calls(a.new)
    r = compare(old, new)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["old_call", "new_call", "proteins"])
        for (o, n), c in sorted(r["transitions"].items(), key=lambda kv: (-kv[1], kv[0])):
            w.writerow([o, n, c])
    print(f"proteins in both: {r['n_both']}; only in old: {r['only_old']}; only in new: {r['only_new']}")
    print(f"changed calls: {sum(r['changed'].values())}")
    for (o, n), c in sorted(r["changed"].items(), key=lambda kv: -kv[1]):
        print(f"  {o}>{n}: {c}")
    print(f"long_pendant: old {lp_old}, new {lp_new}")


if __name__ == "__main__":
    sys.exit(main())
