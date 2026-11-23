#!/usr/bin/env python3
"""Share of genomes per class with a call: rule-cutoff hits, lower-cutoff hits, and proteins placed in the target clade (PNG and PDF).

Reads by_class.tsv from bin/placed_vs_hits_by_class.py. Classes with at least --min-n genomes. Grey: hits at the rule cutoff; orange: hits at
the lower cutoff; blue: proteins placed in the target clade (blue and orange are validated categorical slots; every series is also in the legend).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GREY, ORANGE, BLUE, INK, INK2, SURFACE = "#b9b8b3", "#eb6834", "#2a78d6", "#0b0b0b", "#52514e", "#fcfcfb"


def read_by_class(path, min_n):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["class"] != "ALL" and int(r["n_genomes"]) >= min_n:
                rows.append((r["class"], int(r["n_genomes"]), float(r["share_rule_hits"]), float(r["share_low_hits"]), float(r["share_placed"])))
    return rows


def draw(rows, out_dir, name, rule_cutoff, low_cutoff):
    fig, ax = plt.subplots(figsize=(7.4, 0.5 * len(rows) + 1.6))
    h = 0.26
    for i, (cl, n, a, b, c) in enumerate(rows):
        ax.barh(i - h, 100 * a, height=h, color=GREY)
        ax.barh(i, 100 * b, height=h, color=ORANGE)
        ax.barh(i + h, 100 * c, height=h, color=BLUE)
        ax.text(101.5, i + h, f"{100 * c:.0f}%", va="center", fontsize=7, color=INK2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{cl} (n={n:,})" for cl, n, *_ in rows], fontsize=8)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.set_xlim(0, 108)
    ax.set_xlabel("share of genomes with a call (%)")
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=GREY, label=f"hit at the rule cutoff ({rule_cutoff:g} bits)"),
               plt.Rectangle((0, 0), 1, 1, color=ORANGE, label=f"hit at the lower cutoff ({low_cutoff:g} bits)"),
               plt.Rectangle((0, 0), 1, 1, color=BLUE, label=f"protein placed in the target clade (hits of {low_cutoff:g} bits or more)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.4, -0.12), frameon=False, fontsize=7.5)
    ax.set_title("GGPP synthase: hits against placement in the GGPP synthase clade", loc="left", fontsize=10, fontweight="bold")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--by-class", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-n", type=int, default=150)
    ap.add_argument("--rule-cutoff", type=float, default=369.3)
    ap.add_argument("--low-cutoff", type=float, default=150)
    ap.add_argument("--name", default="fig15_ggpp_placed_vs_hits")
    a = ap.parse_args(argv)
    draw(read_by_class(a.by_class, a.min_n), a.out_dir, a.name, a.rule_cutoff, a.low_cutoff)
    print("wrote", a.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
