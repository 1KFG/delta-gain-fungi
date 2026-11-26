#!/usr/bin/env python3.12
"""Survey of the epa-ng placement cap (default: at most 7 placements per query) in step 5 runs.

survey: for each run, count leaves and distinct leaf sequences of the trimmed reference alignment (what epa-ng sees),
  queries, calls (target / other side / unassigned, from placements_detail), and queries that hit the cap
  (n_placements == --cap in assign.tsv). Reads result files only.
compare: old vs new placements_detail of one step: counts per class, proteins whose target-vs-other call changed.
Run spec for survey: LABEL=RUNDIR[:STEP]. With :STEP the files are assign.STEP.tsv and placements_detail.STEP.tsv.
Nothing here makes a call. Reads plain or .gz/.zst input where noted.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step5_compare_calls import open_text  # noqa: E402


def fasta_seqs(path):
    seqs, cur = [], None
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = []
                seqs.append(cur)
            elif cur is not None:
                cur.append(line.strip())
    return ["".join(s).upper() for s in seqs]


def leaf_counts(afa):
    """-> (n leaves, n distinct sequences) of an aligned FASTA (gap characters kept)."""
    s = fasta_seqs(afa)
    return len(s), len(set(s))


def call_class(call):
    if call == "target":
        return "target"
    if call in ("unassigned", "", None):
        return "unassigned"
    return "other"


def read_classes(detail):
    """-> {protein_id: 'target'|'other'|'unassigned'} from a placements_detail.tsv"""
    out = {}
    with open_text(detail) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["protein_id"]] = call_class(r["clade_call"])
    return out


def cap_hits(assign, proteins, cap):
    """-> (queries in assign among proteins, queries with n_placements >= cap, those among them unassigned)"""
    n = hit = hit_un = 0
    with open_text(assign) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if proteins is not None and r["query"] not in proteins:
                continue
            n += 1
            if int(r["n_placements"]) >= cap:
                hit += 1
                hit_un += r["our_clade"] == "unassigned"
    return n, hit, hit_un


def survey_one(label, rundir, step, cap=7):
    d = Path(rundir)
    sfx = f".{step}" if step else ""
    cls = read_classes(d / f"placements_detail{sfx}.tsv")
    c = Counter(cls.values())
    afa = d / "ref_set.trim.afa"
    if not afa.exists():  # sub-folder of a run that reuses the parent tree (for example per_step_decoys/)
        afa = d.parent / "ref_set.trim.afa"
    leaves, distinct = leaf_counts(afa)
    n_assign, hit, hit_un = cap_hits(d / f"assign{sfx}.tsv", set(cls), cap)
    n = len(cls)
    return {"run": label, "leaves": leaves, "distinct_leaf_seqs": distinct, "queries": n, "queries_in_assign": n_assign,
            "target": c["target"], "other": c["other"], "unassigned": c["unassigned"],
            "unassigned_frac": f"{c['unassigned'] / n:.4f}" if n else "NA", "cap_hit_queries": hit,
            "cap_hit_frac": f"{hit / n:.4f}" if n else "NA", "cap_hit_unassigned": hit_un}


def compare_classes(old, new):
    both = old.keys() & new.keys()
    trans = Counter((old[p], new[p]) for p in both)
    tvo = sum(v for (a, b), v in trans.items() if a != b and "unassigned" not in (a, b))
    return {"n_both": len(both), "only_old": len(old.keys() - new.keys()), "only_new": len(new.keys() - old.keys()),
            "transitions": trans, "target_other_flips": tvo, "changed": sum(v for (a, b), v in trans.items() if a != b)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("survey")
    s.add_argument("runs", nargs="+")
    s.add_argument("--cap", type=int, default=7)
    s.add_argument("--out", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--old", required=True, help="old placements_detail.tsv")
    c.add_argument("--new", required=True)
    c.add_argument("--label", required=True)
    c.add_argument("--out", required=True, help="TSV: label, old_class, new_class, n")
    a = ap.parse_args()
    if a.cmd == "survey":
        rows = []
        for spec in a.runs:
            label, rest = spec.split("=", 1)
            rundir, _, step = rest.partition(":")
            rows.append(survey_one(label, rundir, step or None, a.cap))
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    else:
        old, new = read_classes(a.old), read_classes(a.new)
        r = compare_classes(old, new)
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["label", "old_class", "new_class", "n"])
            for (x, y), v in sorted(r["transitions"].items()):
                w.writerow([a.label, x, y, v])
        print({k: v for k, v in r.items() if k != "transitions"})


if __name__ == "__main__":
    main()
