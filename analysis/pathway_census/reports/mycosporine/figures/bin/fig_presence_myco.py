#!/usr/bin/env python3
"""Presence of the mycosporine steps by phylum and class: fraction of searched genomes with at least one call.

Input: two calls tables written by bin/make_calls_table.py (run with --prefix-aliases and a reps file that lists every
searched genome; see make_myco_calls.sh):
  --calls-hits    source `hits`   for ddgs, omt, atp_grasp_ligase (search hit at the cutoff; NOT placed)
  --calls-placed  source `placed` for ddgs (ddgs_allplace run) and omt (omt_v2 run, single reference); `assigned` rows only
Writes presence_by_phylum.tsv and presence_by_class.tsv (numbers behind the figure and the report tables) and a heatmap
of classes with at least --min-n genomes (PNG 150 dpi and PDF). Every cell is percent of the searched genomes of the group.
Genomes without a label for the rank form their own row. A call is at least one protein; absence is not shown by this table.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

BLUE_RAMP = LinearSegmentedColormap.from_list("blue", ["#eef3fb", "#9ecae1", "#2171b5", "#08306b"])
GREY_RAMP = LinearSegmentedColormap.from_list("grey", ["#f3f2ef", "#c9c8c2", "#8a8985", "#52514e"])
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
# column key, source table, step, label, group (placed = blue panel, hit = grey panel)
COLUMNS = [("ddgs_placed", "placed", "ddgs", "ddgs\nplaced", "placed"),
           ("omt_placed", "placed", "omt", "omt target\n(single ref.)", "placed"),
           ("ddgs_hit", "hits", "ddgs", "ddgs\nsearch hit", "hit"),
           ("omt_hit", "hits", "omt", "omt\nsearch hit", "hit"),
           ("atp_grasp_ligase_hit", "hits", "atp_grasp_ligase", "atp_grasp_ligase\nsearch hit", "hit")]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table(hits, placed, rank):
    """rank -> list of dict(group, n, <column key>: count) over searched genomes."""
    src = {"hits": {r["locustag"]: r for r in hits}, "placed": {r["locustag"]: r for r in placed}}
    agg = defaultdict(lambda: defaultdict(int))
    for g, r in src["hits"].items():
        if r["searched"] != "1":
            continue
        grp = (r["PHYLUM"] or "(no label)", r[rank] or "(no label)") if rank == "CLASS" else r[rank] or "(no label)"
        agg[grp]["n"] += 1
        for key, tab, step, _, _ in COLUMNS:
            if src[tab][g][step] == "1":
                agg[grp][key] += 1
    # for the class rank the key is (phylum, class): a missing class label can sit in several phyla
    rows = [({"phylum": g[0], "group": g[1]} if isinstance(g, tuple) else {"group": g}) | dict(v) for g, v in agg.items()]
    for r in rows:
        for key, *_ in COLUMNS:
            r.setdefault(key, 0)
    return sorted(rows, key=lambda r: -r["n"])


def write(rows, path, parent=None):
    cols = (["phylum"] if parent else []) + ["group", "n_genomes"] + [k for k, *_ in COLUMNS]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow(([r.get("phylum", "")] if parent else []) + [r["group"], r["n"]] + [r[k] for k, *_ in COLUMNS])


def heatmap(rows, out_dir, name, min_n):
    rows = [r for r in rows if r["n"] >= min_n]
    rows.sort(key=lambda r: (r["phylum"], -r["n"]))
    pct = np.array([[100.0 * r[k] / r["n"] for k, *_ in COLUMNS] for r in rows])
    fig, ax = plt.subplots(figsize=(7.6, 0.27 * len(rows) + 1.9))
    ax.set_facecolor(SURFACE)
    fig.patch.set_facecolor(SURFACE)
    for j, (_, _, _, _, grp) in enumerate(COLUMNS):
        cmap = BLUE_RAMP if grp == "placed" else GREY_RAMP
        ax.imshow(np.where(np.arange(len(COLUMNS)) == j, pct, np.nan), cmap=cmap, vmin=0, vmax=100, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(COLUMNS)):
            v = pct[i, j]
            ax.text(j, i, "<0.5" if 0 < v < 0.5 else f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                    color="white" if v > 55 else INK)
    ax.set_xticks(range(len(COLUMNS)))
    ax.set_xticklabels([c[3] for c in COLUMNS], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f'{r["group"]} ({r["n"]:,}) [{r["phylum"][:4]}]' for r in rows], fontsize=8)
    ax.axvline(1.5, color=SURFACE, lw=4)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(out_dir) / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(Path(out_dir) / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls-hits", required=True, type=Path)
    ap.add_argument("--calls-placed", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-n", type=int, default=25)
    ap.add_argument("--name", default="fig_presence_by_class")
    a = ap.parse_args(argv)
    hits, placed = read_tsv(a.calls_hits), read_tsv(a.calls_placed)
    ph = table(hits, placed, "PHYLUM")
    write(ph, a.out_dir / "presence_by_phylum.tsv")
    cl = table(hits, placed, "CLASS")
    write(cl, a.out_dir / "presence_by_class.tsv", parent=True)
    heatmap(cl, a.out_dir, a.name, a.min_n)
    print("genomes searched:", sum(r["n"] for r in ph), "phyla:", len(ph), "classes:", len(cl))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
