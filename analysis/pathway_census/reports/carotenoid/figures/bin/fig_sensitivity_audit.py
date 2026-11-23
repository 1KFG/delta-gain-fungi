#!/usr/bin/env python3
"""Figure for the sensitivity audit: score of every reviewed fungal entry of the step's EC against the step HMM, with the cutoff (PNG and PDF).

Reads sensitivity_entries.tsv and sensitivity_summary.tsv from bin/sensitivity_audit.py. A blue circle is an entry at or above the
cutoff, an orange diamond an entry below it (also shown when the HMM reports no hit: drawn at the left edge as an open diamond).
Steps with fewer than --min-entries entries are listed as "too few to judge" and not drawn.
"""
from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401

BLUE, ORANGE, INK, INK2, GREY, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8985", "#fcfcfb"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def split_steps(summary, min_entries):
    """Return (steps to draw, steps with too few entries) in summary order."""
    draw = [s for s in summary if int(float(s["n_entries"] or 0)) >= min_entries]
    few = [s for s in summary if int(float(s["n_entries"] or 0)) < min_entries]
    return draw, few


def draw(entries, summary, out_dir, name="fig12_sensitivity_audit", min_entries=5):
    draw_steps, few = split_steps(summary, min_entries)
    fig, ax = plt.subplots(figsize=(7.4, 0.9 * max(len(draw_steps), 1) + 1.8))
    rnd = random.Random(1)
    for i, s in enumerate(draw_steps):
        ee = [e for e in entries if e["step"] == s["step"]]
        for e in ee:
            sc = num(e["score"])
            y = i + rnd.uniform(-0.28, 0.28)
            if sc is None:
                ax.scatter(1, y, marker="D", s=22, facecolors="none", edgecolors=ORANGE, linewidth=1)
            elif e["at_cutoff"] == "True":
                ax.scatter(sc, y, marker="o", s=22, color=BLUE, edgecolor=SURFACE, linewidth=0.5, zorder=3)
            else:
                ax.scatter(sc, y, marker="D", s=24, color=ORANGE, edgecolor=SURFACE, linewidth=0.5, zorder=3)
        c = num(s["cutoff"])
        ax.plot([c, c], [i - 0.4, i + 0.4], color=INK, lw=2)
        rec = num(s["recall_at_cutoff"])
        ax.text(1.02, i, f'{s["n_at_cutoff"]}/{s["n_entries"]}  ({100 * rec:.0f}%)', transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(draw_steps)))
    ax.set_yticklabels([f'{s["step"]}\nEC {s["ec"]}' for s in draw_steps], fontsize=8)
    ax.set_ylim(len(draw_steps) - 0.5, -0.6)
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200, 500])
    ax.set_xticklabels(["20", "50", "100", "200", "500"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("best HMM domain score of the entry (log scale)")
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([], [], marker="o", ls="", color=BLUE, label="entry at or above the cutoff"),
               plt.Line2D([], [], marker="D", ls="", color=ORANGE, label="entry below the cutoff"),
               plt.Line2D([], [], marker="|", ls="", color=INK, markersize=12, markeredgewidth=2, label="cutoff used")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3, frameon=False, fontsize=8)
    ax.set_title("Do the step cutoffs find reviewed fungal entries of the same EC?", loc="left", fontsize=10, fontweight="bold")
    if few:
        ax.text(0, -0.42, "Too few reviewed entries to judge (n): " + ", ".join(f'{s["step"]} ({s["n_entries"]})' for s in few),
                transform=ax.transAxes, fontsize=7, color=INK2, ha="left", va="top")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return [s["step"] for s in draw_steps], [s["step"] for s in few]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--audit-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-entries", type=int, default=5)
    a = ap.parse_args(argv)
    drawn, few = draw(read_tsv(a.audit_dir / "sensitivity_entries.tsv"), read_tsv(a.audit_dir / "sensitivity_summary.tsv"), a.out_dir, min_entries=a.min_entries)
    print(f"drew {drawn}; too few: {few}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
