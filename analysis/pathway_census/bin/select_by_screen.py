#!/usr/bin/env python3
"""Write the FASTA records that pass an alignability screen (bin/screen_outgroups.py output) for one HMM.

--mode best uses pass_best; --mode union uses pass_union (union passes include every best pass).
--ids limits the records to a list of ids (comma list); --exclude removes ids (comma list).
Record order follows the input FASTA. Prints the kept ids. Nothing under refs/ is changed.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_fasta(path):
    name, hdr, seq = None, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, hdr, "".join(seq)
            hdr = line[1:]
            name, seq = hdr.split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, hdr, "".join(seq)


def passing(screen_tsv, hmm, mode):
    col = {"best": "pass_best", "union": "pass_union"}[mode]
    with open(screen_tsv) as fh:
        return {r["protein"] for r in csv.DictReader(fh, delimiter="\t") if r["hmm"] == hmm and r[col] == "True"}


def select(records, ok, ids=None, exclude=()):
    return [(n, h, s) for n, h, s in records if n in ok and (ids is None or n in ids) and n not in set(exclude)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--screen", required=True, type=Path)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--mode", choices=["best", "union"], required=True)
    ap.add_argument("--ids", help="comma list: keep only these")
    ap.add_argument("--exclude", default="", help="comma list: drop these")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    ids = set(a.ids.split(",")) if a.ids else None
    kept = select(list(read_fasta(a.fasta)), passing(a.screen, a.hmm, a.mode), ids, [x for x in a.exclude.split(",") if x])
    with open(a.out, "w") as fh:
        for n, h, s in kept:
            fh.write(f">{h}\n{s}\n")
    print(f"{a.hmm} {a.mode}: kept {len(kept)}: {' '.join(n for n, _, _ in kept)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
