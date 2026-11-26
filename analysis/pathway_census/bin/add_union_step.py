#!/usr/bin/env python3
"""Add a derived "any of" step to calls.tsv.

Usage: add_union_step.py --calls calls.tsv --name any_dhn_pks --members thn_synthase hexaketide_pks heptaketide_pks --out calls_union.tsv

Value per genome: 1 if any member is 1; NA if any member is NA and no member is 1; else 0. The new column
`<name>_source` is `placed` when every member is placed, else `hits`. Member columns are kept.
Tip rows with NA for one member are NA for all members in practice (the same searched set), so the NA rule is a guard.
"""
import argparse
import csv


def union_value(vals):
    if "1" in vals:
        return "1"
    if "NA" in vals:
        return "NA"
    return "0"


def add_union(rows, name, members):
    """rows: list of dicts (calls.tsv). Returns (new rows, new column names in order)."""
    out = []
    for r in rows:
        r = dict(r)
        missing = [m for m in members if m not in r]
        if missing:
            raise KeyError(f"member columns not in calls table: {missing}")
        r[name] = union_value([r[m] for m in members])
        srcs = {r.get(f"{m}_source", "") for m in members}
        r[f"{name}_source"] = "placed" if srcs == {"placed"} else "hits"
        out.append(r)
    return out, [name, f"{name}_source"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--members", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--drop-members", action="store_true",
                    help="remove the member columns and their _source columns from the output")
    a = ap.parse_args(argv)
    with open(a.calls, newline="") as fh:
        rd = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        cols = list(rd.fieldnames)
        rows = list(rd)
    new, extra = add_union(rows, a.name, a.members)
    if a.drop_members:
        gone = set(a.members) | {f"{m}_source" for m in a.members}
        cols = [c for c in cols if c not in gone]
        new = [{k: v for k, v in r.items() if k not in gone} for r in new]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + extra, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(new)
    n1 = sum(1 for r in new if r[a.name] == "1")
    n0 = sum(1 for r in new if r[a.name] == "0")
    nn = sum(1 for r in new if r[a.name] == "NA")
    print(f"{a.name}: present {n1}, absent {n0}, NA {nn}")


if __name__ == "__main__":
    main()
