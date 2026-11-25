#!/usr/bin/env python3.12
"""Combine the per-run report tables of step 5 runs into three small files (one-off helper, not tested code).

usage: collate_tables.py OUT_DIR RUN_DIR [RUN_DIR ...]   (run from analysis/pathway_census; each RUN_DIR has report/)
Writes OUT_DIR/other_steps_assigned_by_class.tsv (rows with at least one assigned genome), OUT_DIR/other_steps_nearest_summary.tsv,
OUT_DIR/other_steps_leaves.tsv (labelled leaves of each run with role).
"""
import csv, sys
from pathlib import Path

out = Path(sys.argv[1])
by_class, nearest, leaves = [], [], []
for d in map(Path, sys.argv[2:]):
    run = d.name
    for r in csv.DictReader(open(d / "report" / "assigned_by_class.tsv"), delimiter="\t"):
        if int(r["n_with_assigned"]) > 0:
            by_class.append([run, r["level"], r["group"], r["n_genomes"], r["n_with_assigned"], r["share"]])
    for r in csv.DictReader(open(d / "report" / "nearest_summary.tsv"), delimiter="\t"):
        nearest.append([run, r["nearest_role"], r["n_representatives"], r["n_candidates_in_their_clusters"], r["median_distance"]])
    for r in csv.DictReader(open(d / "clade_members.tsv"), delimiter="\t"):
        if r["role"] not in ("representative", "other"):
            leaves.append([run, r["label"], r["role"], r["our_clade"] or "(none)"])
for name, cols, rows in (("other_steps_assigned_by_class.tsv", ["run", "level", "group", "n_genomes", "n_with_assigned", "share"], by_class),
                         ("other_steps_nearest_summary.tsv", ["run", "nearest_role", "n_representatives", "n_candidates_in_their_clusters", "median_distance"], nearest),
                         ("other_steps_leaves.tsv", ["run", "leaf", "role", "clade"], leaves)):
    with open(out / name, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)
    print(name, len(rows))
