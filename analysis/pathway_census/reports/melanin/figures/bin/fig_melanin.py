#!/usr/bin/env python3
"""Figures of the melanin report (PNG and PDF) from the tables of bin/summarise_calls.py.

fig1  placed-call presence of each step by class (percent of searched genomes)
fig2  DHN only / L-DOPA only / both / neither per species, by class and by order (headline rule pair)
fig3  the four classes under each rule pair, genome level and species level (rule sensitivity)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
BLUE_RAMP = LinearSegmentedColormap.from_list("blue", ["#eef3fb", "#9ecae1", "#2171b5", "#08306b"])
GREY_RAMP = LinearSegmentedColormap.from_list("grey", ["#f3f2ef", "#c9c8c2", "#8a8985", "#52514e"])
# four classes: blue, orange, aqua from the report palette, grey for neither
CLASS_COL = {"DHN_only": "#2a78d6", "LDOPA_only": "#eb6834", "both": "#1baf7a", "neither": "#c9c8c2"}
CLASS_LAB = {"DHN_only": "DHN only", "LDOPA_only": "L-DOPA only", "both": "both", "neither": "neither"}
DHN = ["thn_synthase", "hexaketide_pks", "heptaketide_pks", "ayg1", "t4hn_reductase", "t3hn_reductase",
       "scytalone_dehydratase", "laccase_abr2", "mco_abr1"]
LDOPA = ["tyrosinase", "laccase_diphenol", "aspulvinone_e_synthase", "asp_melanin_tyrosinase"]
FAM = ["hppd", "hmga"]
ANY_COLS = ["D_pks", "L_route"]


def save(fig, out_dir, name):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{name}.{ext}", dpi=160, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def cell(v):
    if v != v:
        return ""
    if v >= 0.005:
        return f"{int(100 * v + 0.5)}"
    return "<0.5" if v > 0 else ""


def fig1(step_by_class, out_dir, min_n):
    d = step_by_class[step_by_class.n_genomes >= min_n].copy()
    d = d.sort_values(["PHYLUM", "n_genomes"], ascending=[True, False])
    d["order_key"] = d.PHYLUM.map(d.groupby("PHYLUM").n_genomes.sum())
    d = d.sort_values(["order_key", "PHYLUM", "n_genomes"], ascending=[False, True, False])
    groups = [("DHN steps (placed)", DHN, BLUE_RAMP), ("L-DOPA-pathway steps (placed)", LDOPA, BLUE_RAMP),
              ("family present, not a placed step", FAM, GREY_RAMP)]
    ratios = [len(g[1]) for g in groups]
    fig, axes = plt.subplots(1, 3, figsize=(1.0 + 0.62 * sum(ratios) + 3.0, 0.3 * len(d) + 2.2),
                             gridspec_kw={"width_ratios": ratios, "wspace": 0.05}, facecolor=SURFACE)
    for i, (ax, (title, cols, ramp)) in enumerate(zip(axes, groups)):
        m = np.array([[r[c] / r["n_genomes"] for c in cols] for _, r in d.iterrows()])
        ax.imshow(m, cmap=ramp, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8, color=INK)
        ax.set_title(title, fontsize=8.5, color=INK2)
        for (y, x), v in np.ndenumerate(m):
            ax.text(x, y, cell(v), ha="center", va="center", fontsize=7, color="white" if v > 0.55 else INK)
        ax.set_yticks(range(len(d)) if i == 0 else [])
        if i == 0:
            ax.set_yticklabels([f"{c} ({n:,})" for c, n in zip(d.CLASS, d.n_genomes)], fontsize=8, color=INK)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    save(fig, out_dir, "fig1_step_presence_by_class")


def stacked(ax, d, label_col, species=True, title=""):
    pre = "species_" if species else "genome_"
    n = d[pre.rstrip("_").replace("species", "n_species").replace("genome", "n_genomes")] if False else d["n_species" if species else "n_genomes"]
    left = np.zeros(len(d))
    for k in ("DHN_only", "LDOPA_only", "both", "neither"):
        v = d[pre + k].to_numpy() / n.to_numpy()
        ax.barh(range(len(d)), v, left=left, color=CLASS_COL[k], label=CLASS_LAB[k], height=0.78, edgecolor=SURFACE, linewidth=0.4)
        for y, (l0, w) in enumerate(zip(left, v)):
            if w >= 0.08:
                ax.text(l0 + w / 2, y, f"{int(100 * w + 0.5)}", ha="center", va="center", fontsize=7, color=INK)
        left = left + v
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f"{x} ({int(k):,})" for x, k in zip(d[label_col], n)], fontsize=8, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=8)
    ax.set_xlabel("percent of species" if species else "percent of genomes", fontsize=8.5, color=INK2)
    ax.set_title(title, fontsize=9, color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig2(by_class, by_order, out_dir, min_class, min_order):
    c = by_class[by_class.n_species >= min_class].sort_values("n_species", ascending=False)
    o = by_order[by_order.n_species >= min_order].sort_values("n_species", ascending=False)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 0.27 * max(len(c), len(o)) + 1.8), facecolor=SURFACE,
                             gridspec_kw={"wspace": 0.75})
    stacked(axes[0], c, "CLASS", True, f"Class (species >= {min_class}; n species in brackets)")
    stacked(axes[1], o, "ORDER", True, f"Order (species >= {min_order})")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=8.5)
    save(fig, out_dir, "fig2_four_class_by_clade")


def fig3(sens, out_dir):
    note = {"D_pks|L_route": "headline: any DHN PKS | L-DOPA route", "D_any|L_any": "wide: any placed step | any placed step",
            "D_route|L_route": "strict: complete yaml PKS route | L-DOPA route", "D_pks|L_nontyr": "any DHN PKS | L-DOPA without tyrosinase",
            "D_route_ox|L_route": "strict route + oxidase | L-DOPA route", "D_hit_any|L_hit_any": "search-hit level (not placed)"}
    sens = sens.copy()
    sens["pair"] = sens.DHN_rule + "|" + sens.LDOPA_rule
    pairs = list(dict.fromkeys(sens.pair))
    rows = [(p, lvl) for p in pairs for lvl in ("species_any", "genome")]
    fig, ax = plt.subplots(figsize=(9.5, 0.3 * len(rows) + 1.6), facecolor=SURFACE)
    for y, (p, lvl) in enumerate(rows):
        r = sens[(sens.pair == p) & (sens.level == lvl)].iloc[0]
        left = 0.0
        for k in ("DHN_only", "LDOPA_only", "both", "neither"):
            w = r[k] / r["n"]
            ax.barh(y, w, left=left, color=CLASS_COL[k], height=0.8, edgecolor=SURFACE, linewidth=0.4,
                    label=CLASS_LAB[k] if y == 0 else None)
            if w >= 0.06:
                ax.text(left + w / 2, y, f"{int(r[k]):,}", ha="center", va="center", fontsize=7, color=INK)
            left += w
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{note.get(p, p)} - {'species' if lvl == 'species_any' else 'genomes'}" for p, lvl in rows], fontsize=7.5, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of species (any genome positive) or of searched genomes; counts in the bars", fontsize=8.5, color=INK2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=4, frameon=False, fontsize=8.5)
    save(fig, out_dir, "fig3_rule_sensitivity")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-genomes-class", type=int, default=25)
    ap.add_argument("--min-species-class", type=int, default=20)
    ap.add_argument("--min-species-order", type=int, default=50)
    a = ap.parse_args(argv)
    t = Path(a.tables)
    fig1(pd.read_csv(t / "step_presence_by_class.tsv", sep="\t"), a.out_dir, a.min_genomes_class)
    fig2(pd.read_csv(t / "four_class_by_class.tsv", sep="\t"), pd.read_csv(t / "four_class_by_order.tsv", sep="\t"),
         a.out_dir, a.min_species_class, a.min_species_order)
    fig3(pd.read_csv(t / "sensitivity_four_class.tsv", sep="\t"), a.out_dir)


if __name__ == "__main__":
    main()
