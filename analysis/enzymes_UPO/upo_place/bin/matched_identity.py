#!/usr/bin/env python3
"""Matched-identity validation of placement clade calls on the 187-leaf reference.

The same design as the skylark session's test (PLACEMENT_FINDINGS.md, job 31868181), so the two
reference sets can be compared on equal terms: for each query and identity ceiling, remove the
query AND every other reference whose identity to it exceeds the ceiling, place the query on what
is left, and score the call against identity to the NEAREST SURVIVING reference. Plain
leave-one-out (ceiling = none) keeps close relatives in the tree and so tests an easier regime.

  cases      reference MSA + clade table -> cases.tsv (one row per query x ceiling)
  summarize  cases.tsv + per-case assign outputs -> precision table (markdown)

Identity is over the reference alignment's columns where both sequences have a residue. A case is
`impossible` when the removal leaves no member of the query's expected clade in the tree (the
clade cannot be called); skylark counted 71 of 114 such cases with 19 references. Expected clades
for the 19 characterized enzymes are their published clades; for Zamocky leaves they are the
tree-derived labels (name_clades.py), so the characterized subset is the independent test.
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from upo_place import read_fasta  # noqa: E402

BINS = [(70, 101, "70-100"), (55, 70, "55-70"), (45, 55, "45-55"), (35, 45, "35-45"),
        (0, 35, "<35")]


def identity(a, b):
    n = m = 0
    for x, y in zip(a, b):
        if x != "-" and y != "-":
            n += 1
            m += x == y
    return 100.0 * m / n if n else 0.0


def load_expected(path):
    with open(path, newline="") as fh:
        return {r["label"]: (r["our_clade"] or "unassigned") for r in csv.DictReader(fh, delimiter="\t")}


def cases(msa, clades, ceilings):
    seqs = dict(read_fasta(msa))
    names = list(seqs)
    exp = load_expected(clades)
    members = defaultdict(set)
    for k, c in exp.items():
        if c != "unassigned":
            members[c].add(k)
    ident = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ident[a, b] = ident[b, a] = identity(seqs[a], seqs[b])
    rows = []
    for q in names:
        for ceil in ceilings:
            removed = [o for o in names if o != q and ceil is not None and ident[q, o] > ceil]
            surv = [o for o in names if o != q and o not in removed]
            near = max(surv, key=lambda o: ident[q, o])
            e = exp[q]
            possible = e == "unassigned" or bool((members[e] - {q}) - set(removed))
            rows.append({
                "case": f"{q}@{'loo' if ceil is None else ceil}", "query": q,
                "ceiling": "loo" if ceil is None else ceil, "set":
                    "characterized" if q.startswith("B__") else "zamocky",
                "expected": e, "possible": "yes" if possible else "no",
                "n_removed": len(removed), "nearest_surviving": near,
                "nearest_identity": f"{ident[q, near]:.1f}", "removed": ",".join(removed)})
    return rows


def bin_of(x):
    for lo, hi, name in BINS:
        if lo <= x < hi:
            return name
    return "<35"


def summarize(cases_tsv, assign_dir):
    with open(cases_tsv, newline="") as fh:
        cs = list(csv.DictReader(fh, delimiter="\t"))
    tab = defaultdict(lambda: defaultdict(int))
    missing = 0
    for c in cs:
        key = (c["set"], c["ceiling"] == "loo")
        if c["possible"] == "no":
            tab[key, "impossible"]["n"] += 1
            continue
        p = os.path.join(assign_dir, c["case"].replace("@", "__at__"), "assign.tsv")
        if not os.path.exists(p):
            missing += 1
            continue
        with open(p, newline="") as fh:
            got = next(csv.DictReader(fh, delimiter="\t"))["our_clade"]
        b = bin_of(float(c["nearest_identity"]))
        t = tab[key, b]
        t["n"] += 1
        if got == c["expected"]:
            t["correct"] += 1
        elif got == "unassigned":
            t["declined"] += 1          # expected a clade, abstained
        else:
            t["wrong"] += 1             # wrong clade, or a clade for an outside leaf
    out = [f"Cases: {len(cs)}; missing outputs: {missing}.\n"]
    for setname in ("characterized", "zamocky"):
        for loo in (False, True):
            title = f"{setname}, {'plain leave-one-out' if loo else 'matched identity (ceilings 70/55/45/35)'}"
            out += [f"### {title}\n", "| identity to nearest surviving reference | n | correct | "
                    "declined | wrong | precision when committing |", "|---|---|---|---|---|---|"]
            tot = defaultdict(int)
            for _, _, b in BINS:
                t = tab[(setname, loo), b]
                if not t["n"]:
                    continue
                com = t["correct"] + t["wrong"]
                out.append(f"| {b} | {t['n']} | {t['correct']} | {t['declined']} | {t['wrong']} | "
                           f"{100 * t['correct'] / com:.0f} % |" if com else
                           f"| {b} | {t['n']} | {t['correct']} | {t['declined']} | {t['wrong']} | - |")
                for k in ("n", "correct", "declined", "wrong"):
                    tot[k] += t[k]
            com = tot["correct"] + tot["wrong"]
            out.append(f"| **all** | {tot['n']} | {tot['correct']} | {tot['declined']} | "
                       f"{tot['wrong']} | {100 * tot['correct'] / com:.0f} % |" if com else "")
            out.append(f"\nImpossible cases (expected clade removed entirely): "
                       f"{tab[(setname, loo), 'impossible']['n']}\n")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cases")
    c.add_argument("--msa", required=True)
    c.add_argument("--clades", required=True)
    c.add_argument("--ceilings", default="loo,70,55,45,35")
    c.add_argument("--out", required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--cases", required=True)
    s.add_argument("--assign-dir", required=True)
    s.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "cases":
        ceils = [None if x == "loo" else float(x) for x in a.ceilings.split(",")]
        rows = cases(a.msa, a.clades, ceils)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"cases: {len(rows)} ({sum(r['possible'] == 'no' for r in rows)} impossible)",
              file=sys.stderr)
    else:
        md = summarize(a.cases, a.assign_dir)
        open(a.out, "w").write(md + "\n")
        print(md)


if __name__ == "__main__":
    main()
