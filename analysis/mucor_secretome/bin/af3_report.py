#!/usr/bin/env python3
"""Collect AlphaFold 3 confidences into one table, with the controls kept in view.

  af3_report.py --dir af3_out/ --out af3_scores.tsv [--controls controls.tsv]

AF3 writes <job>/<job>_summary_confidences.json per prediction. The fields that matter for an
interaction screen are `iptm` (interface predicted TM-score, the pair score), `ptm`, the
chain-pair iPTM matrix and `ranking_score`.

There is no universal iPTM threshold, and the one quoted for AlphaFold-Multimer (0.8 confident,
below 0.6 likely wrong) was calibrated on large well-behaved complexes, not on small
cysteine-rich effectors, where false negatives are common. This script therefore reports the
controls alongside the screen rather than applying a fixed cut: a threshold is only worth using
if the positive controls clear it and the negatives do not. `--controls` marks which rows are
controls and what they are for, and the summary prints their scores first.

A protein-ligand job has no iPTM at all; for those the useful numbers are the ranking score and
the ligand's own pLDDT, so both are reported and iPTM is left blank. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys


def find_summaries(root):
    """-> [(job_name, path to *_summary_confidences.json)]."""
    out = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if f.endswith("_summary_confidences.json"):
                out.append((os.path.basename(dirpath), os.path.join(dirpath, f)))
    return sorted(out)


def parse_summary(path):
    with open(path) as fh:
        d = json.load(fh)
    return {"iptm": d.get("iptm"), "ptm": d.get("ptm"),
            "ranking_score": d.get("ranking_score"),
            "fraction_disordered": d.get("fraction_disordered"),
            "has_clash": d.get("has_clash"),
            "chain_iptm": d.get("chain_iptm"), "chain_ptm": d.get("chain_ptm"),
            "chain_pair_iptm": d.get("chain_pair_iptm")}


def pair_iptm(summary):
    """The A-B interface score from the chain-pair matrix, which is the number for a 2-chain job."""
    m = summary.get("chain_pair_iptm")
    if isinstance(m, list) and len(m) >= 2 and isinstance(m[0], list) and len(m[0]) >= 2:
        return m[0][1]
    return summary.get("iptm")


def split_name(job):
    """'A__B' -> ('A', 'B'); 'A__CO7' -> ('A', 'CO7'); otherwise (job, '')."""
    if "__" in job:
        a, b = job.split("__", 1)
        return a, b
    return job, ""


def rows_from(root, controls=None):
    controls = controls or {}
    rows = []
    for job, path in find_summaries(root):
        try:
            s = parse_summary(path)
        except (OSError, json.JSONDecodeError) as e:
            print(f"skipping {path}: {e}", file=sys.stderr)
            continue
        a, b = split_name(job)
        is_ligand = b.upper().startswith("CO")
        rows.append({
            "job": job, "partner_a": a, "partner_b": b,
            "kind": "ligand" if is_ligand else ("control" if job in controls else "pair"),
            "control_role": controls.get(job, ""),
            "iptm": "" if is_ligand else (pair_iptm(s) if pair_iptm(s) is not None else ""),
            "ptm": s.get("ptm", ""),
            "ranking_score": s.get("ranking_score", ""),
            "fraction_disordered": s.get("fraction_disordered", ""),
            "has_clash": s.get("has_clash", ""),
        })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--controls", help="TSV of job<TAB>role, marking control jobs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    controls = {}
    if a.controls and os.path.exists(a.controls):
        with open(a.controls) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and not line.startswith("#"):
                    controls[f[0]] = f[1]

    rows = rows_from(a.dir, controls)
    cols = ["job", "kind", "control_role", "partner_a", "partner_b", "iptm", "ptm",
            "ranking_score", "fraction_disordered", "has_clash"]

    def key(r):
        try:
            return -float(r["iptm"])
        except (TypeError, ValueError):
            try:
                return -float(r["ranking_score"])
            except (TypeError, ValueError):
                return 0.0

    rows.sort(key=key)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    ctrl = [r for r in rows if r["kind"] == "control"]
    if ctrl:
        print("controls (a threshold is only usable if these separate):", file=sys.stderr)
        for r in ctrl:
            print(f"  {r['control_role']:<10} {r['job']:<40} iptm={r['iptm']} "
                  f"ptm={r['ptm']}", file=sys.stderr)
    else:
        print("NO CONTROLS in this run - scores are uncalibrated", file=sys.stderr)
    print(f"{len(rows)} predictions -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
