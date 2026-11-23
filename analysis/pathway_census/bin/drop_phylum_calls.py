#!/usr/bin/env python3
"""Drop every tip of the named phyla: set all steps to NA in calls.tsv and mark the tip (`dropped_phylum` = 1).

Usage: drop_phylum_calls.py --calls calls.tsv --phylum Microsporidia [--phylum X ...] --out calls_dropped.tsv --dropped-list tips.txt
The phylum comes from the PHYLUM column of calls.tsv (species.parquet). Tips already NA stay NA. The list file holds the
LOCUSTAGs of the dropped tips (used by bin/prune_tree_tips.R to prune the tree before reconstruction).
"""
import argparse
import csv
import sys
from pathlib import Path


def drop_phyla(rows, phyla):
    """rows: list of dict rows of calls.tsv. Returns (new rows, list of dropped locustags)."""
    phyla = set(phyla)
    out, dropped = [], []
    for r in rows:
        r = dict(r)
        hit = r.get("PHYLUM", "") in phyla
        r["dropped_phylum"] = "1" if hit else "0"
        if hit:
            for k in list(r):
                if k.endswith("_source"):
                    r[k[: -len("_source")]] = "NA"
            r["searched"] = "0"
            dropped.append(r["locustag"])
        out.append(r)
    return out, dropped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--phylum", action="append", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dropped-list", required=True)
    a = ap.parse_args(argv)
    with open(a.calls, newline="") as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    new, dropped = drop_phyla(rows, a.phylum)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(new[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(new)
    Path(a.dropped_list).write_text("\n".join(dropped) + "\n")
    print(f"dropped {len(dropped)} tips of {', '.join(a.phylum)}", file=sys.stderr)


if __name__ == "__main__":
    main()
