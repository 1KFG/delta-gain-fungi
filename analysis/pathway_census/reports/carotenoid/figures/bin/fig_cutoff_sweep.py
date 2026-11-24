#!/usr/bin/env python3
"""Share of genomes with a call against the score cutoff, by class, for a current and a wider HMM (PNG and PDF).

Reads sweep_by_class.tsv from bin/cutoff_sweep.py (one folder per HMM). Three classes are highlighted with the validated categorical
slots (blue, orange, aqua) and direct labels; other classes are thin grey lines; the all-genome line is dashed black. A vertical
line marks the cutoff in use and one marks the cutoff proposed in the report.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401

BLUE, ORANGE, AQUA, GREY, INK, INK2, SURFACE = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b3", "#0b0b0b", "#52514e", "#fcfcfb"
HIGHLIGHT = {"Saccharomycetes": BLUE, "Pichiomycetes": ORANGE, "Agaricomycetes": AQUA}


def read_sweep(path):
    """Return (cutoffs, {class: [fractions]}, {class: n}) from sweep_by_class.tsv."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    head = rows[0]
    cut = [float(c.replace("frac_ge_", "")) for c in head[2:]]
    order = sorted(range(len(cut)), key=lambda i: cut[i])
    cutoffs = [cut[i] for i in order]
    data, ns = {}, {}
    for r in rows[1:]:
        vals = [float(x) for x in r[2:]]
        data[r[0]] = [vals[i] for i in order]
        ns[r[0]] = int(r[1])
    return cutoffs, data, ns


def panel(ax, path, title, in_use, proposed):
    cutoffs, data, ns = read_sweep(path)
    for cl, v in data.items():
        if cl == "ALL" or cl in HIGHLIGHT:
            continue
        ax.plot(cutoffs, [100 * x for x in v], color=GREY, lw=0.9, zorder=1)
    for cl, col in HIGHLIGHT.items():
        if cl in data:
            ax.plot(cutoffs, [100 * x for x in data[cl]], color=col, lw=2.0, zorder=3, label=f"{cl} (n={ns[cl]})")
    ax.plot(cutoffs, [100 * x for x in data["ALL"]], color=INK, lw=1.6, ls="--", zorder=4, label=f"all genomes (n={ns['ALL']:,})")
    ax.plot([], [], color=GREY, lw=0.9, label="other classes")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=7.5)
    ax.axvline(in_use, color=INK2, lw=1, zorder=2)
    ax.axvline(proposed, color=INK2, lw=1, ls=":", zorder=2)
    ax.text(in_use, 3, " in use", fontsize=7, color=INK2, rotation=90, va="bottom", ha="right")
    ax.text(proposed, 3, " tried", fontsize=7, color=INK2, rotation=90, va="bottom", ha="right")
    ax.set_xscale("log")
    ax.set_xticks([100, 150, 200, 300, 500])
    ax.set_xticklabels(["100", "150", "200", "300", "500"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(-2, 104)
    ax.set_xlabel("score cutoff (bits, log scale)")
    ax.set_title(title, fontsize=9, loc="left", fontweight="bold")
    ax.grid(axis="y", color="#e6e5e1", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    return ns["ALL"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--current", required=True, type=Path)
    ap.add_argument("--wide", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--name", default="fig14_ggpp_cutoff_sweep")
    a = ap.parse_args(argv)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    n = panel(axes[0], a.current, "Current HMM (2 references)", 369.3, 150)
    panel(axes[1], a.wide, "Wider HMM (references and named GGPP synthases)", 314.8, 200)
    axes[0].set_ylabel(f"share of genomes with a call (%), {n:,} genomes")
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
    fig.tight_layout()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out_dir / f"{a.name}.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(a.out_dir / f"{a.name}.pdf", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", a.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
