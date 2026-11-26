#!/usr/bin/env python3
"""Evidence figure for the omt step: HMM score of the candidates by tree call, with the two cutoffs marked.

Reads the placement folder of the single-reference omt run (candidates_step.tsv for the score, placements_detail.tsv for
the call `clade_call`: target, decoy, outgroup or unassigned). Draws a histogram of the score, one series per call group
(target; other side = decoy + outgroup; unassigned), and marks the cutoff used (100.6) and the highest decoy score (120.7).
Writes PNG (150 dpi) and PDF. Run: /usr/bin/python3.12 fig_omt_scores.py --omt-dir <results/step5/mycosporine/omt_v2> --out-dir <dir>
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, GREY, INK, INK2, SURFACE = "#2a78d6", "#eb6834", "#8a8985", "#0b0b0b", "#52514e", "#fcfcfb"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def group(call):
    return "target" if call == "target" else ("unassigned" if call == "unassigned" else "other side")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omt-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--cutoff", type=float, default=100.6)
    ap.add_argument("--decoy-max", type=float, default=120.7)
    ap.add_argument("--name", default="fig_omt_score_by_call")
    a = ap.parse_args(argv)
    score = {}
    for r in read_tsv(a.omt_dir / "candidates_step.tsv"):
        score[r["protein_id"]] = float(r["score"])
    calls = {r["protein_id"]: group(r["clade_call"]) for r in read_tsv(a.omt_dir / "placements_detail.tsv")}
    series = {"target": [], "other side": [], "unassigned": []}
    for pid, c in calls.items():
        if pid in score:
            series[c].append(score[pid])
    plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                         "axes.edgecolor": GREY, "axes.spines.top": False, "axes.spines.right": False,
                         "font.size": 9, "pdf.fonttype": 42, "axes.grid": True, "grid.color": "#e6e5e1"})
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    lo = a.cutoff
    bins = [lo + 4 * i for i in range(int((450 - lo) // 4) + 2)]
    ax.hist([series["target"], series["other side"]], bins=bins, stacked=False, color=[BLUE, ORANGE],
            label=[f"target (n={len(series['target']):,})", f"other side: decoy or outgroup clade (n={len(series['other side']):,})"],
            edgecolor=SURFACE, linewidth=0.4, alpha=0.9)
    ax.set_yscale("log")
    ax.axvline(a.cutoff, color=INK, lw=1.5, label=f"cutoff used ({a.cutoff:g})")
    ax.axvline(a.decoy_max, color=INK2, lw=1.5, ls="--", label=f"highest decoy score ({a.decoy_max:g}, Q8NKC1)")
    ax.set_xlabel("omt HMM bit score (4-bit bins)")
    ax.set_ylabel("candidate proteins (log scale)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("omt: score of each candidate, by call in the single-reference tree", loc="left", fontweight="bold")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out_dir / f"{a.name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(a.out_dir / f"{a.name}.pdf", bbox_inches="tight")
    cnt = {k: len(v) for k, v in series.items()}
    print("counts", cnt, "scored", len(score))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
