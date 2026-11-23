#!/usr/bin/env python3
"""Presence of each pathway step across taxonomic groups: fraction of genomes with at least one call (PNG and PDF).

A genome counts as having a step when the step has at least one call in that genome. A call is
  - a protein placed in the target clade (`assigned` rows from bin/assign_steps.py) for a step that has a placement run
    (--placed STEP=assigned.tsv, repeatable), or
  - otherwise a search hit at the cutoff used (candidates.chunk*.tsv.zst of bin/search_chunk.py). This is NOT a placed call.
Steps whose cutoff is a decoy floor (`basis` says single reference in the min_scores file) hit nearly every genome: they are drawn
in a separate grey panel and must be read as "family present", not "step present".
By default every searched genome is used (--min-busco 0). A BUSCO filter of 90 removes whole lineages with real gene loss or
a poor fit to the BUSCO set (MEASURED 2026-11-22: Saccharomycetes 2,560 genomes to 0, Pichiomycetes 2,929 to 2), so it is an
option and not the default. Groups with fewer than --min-n genomes are dropped. Rank: phylum, subphylum, class, order or family.
"""
from __future__ import annotations

import argparse
import csv
import glob
import io
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

BLUE_RAMP = LinearSegmentedColormap.from_list("blue", ["#eef3fb", "#9ecae1", "#2171b5", "#08306b"])
GREY_RAMP = LinearSegmentedColormap.from_list("grey", ["#f3f2ef", "#c9c8c2", "#8a8985", "#52514e"])
INK, INK2, SURFACE, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
# the three carotenoid-specific steps of the core route; GGPP synthase serves all terpenoid metabolism and is left out
CORE_STEPS = ("phytoene_synthase", "phytoene_desaturase", "lycopene_cyclase")
STEP_ORDER = ["ggpp_synthase", "phytoene_synthase", "phytoene_desaturase", "lycopene_cyclase", "carx", "cart", "ylo1",
              "crtr", "astaxanthin_synthase", "cars", "acaa", "tsp3", "tsp2", "tdh"]
RANKS = {"phylum": "PHYLUM", "subphylum": "SUBPHYLUM", "class": "CLASS", "order": "ORDER", "family": "FAMILY"}
LINEAGE = ["PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY"]


def read_tsv(path):
    path = str(path)
    if path.endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    else:
        raw = Path(path).read_text()
    lines = [l for l in raw.splitlines() if not l.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))


def genome_of(protein_id):
    return protein_id.split("_")[0]


def hits_by_step(search_dir):
    """step -> set of locustags with at least one search hit."""
    out = {}
    for f in sorted(glob.glob(str(Path(search_dir) / "candidates.chunk*.tsv.zst"))):
        for r in read_tsv(f):
            out.setdefault(r["step_id"], set()).add(genome_of(r["protein_id"]))
    return out


def placed_genomes(assigned_tsv):
    """locustags with at least one `assigned` call in an assign_steps output."""
    return {genome_of(r["protein_id"]) for r in read_tsv(assigned_tsv) if r["label"] == "assigned"}


def floor_steps(min_scores_tsv):
    """steps whose cutoff basis is a decoy floor (single reference)."""
    return {r["step"] for r in read_tsv(min_scores_tsv) if "single reference" in (r.get("basis") or "")}


def load_taxonomy(species_parquet):
    import duckdb
    rows = duckdb.connect().execute(
        f'select LOCUSTAG, PHYLUM, SUBPHYLUM, "CLASS", "ORDER", FAMILY from read_parquet(\'{species_parquet}\')').fetchall()
    return {r[0]: dict(zip(LINEAGE, [x or "" for x in r[1:]])) for r in rows}


def presence_table(genomes, taxonomy, calls, steps):
    """Return list of dict(locustag, <lineage>, <step>: 0/1 ...) for the searched genomes."""
    out = []
    for g in genomes:
        tax = taxonomy.get(g)
        if tax is None:
            continue
        row = dict(locustag=g, **tax)
        for s in steps:
            row[s] = int(g in calls.get(s, set()))
        out.append(row)
    return out


def summarise(rows, rank_col, steps, min_n, core=CORE_STEPS):
    """Group rows by rank. Return list of dict(group, lineage, n, frac{step}, core_complete) with n >= min_n, in lineage order."""
    groups = {}
    for r in rows:
        key = r[rank_col] or "(no label)"
        groups.setdefault(key, []).append(r)
    out = []
    for key, rs in groups.items():
        if len(rs) < min_n:
            continue
        lineage = {c: rs[0][c] for c in LINEAGE}
        frac = {s: sum(r[s] for r in rs) / len(rs) for s in steps}
        core_ok = [s for s in core if s in steps]
        out.append(dict(group=key, lineage=lineage, n=len(rs), frac=frac,
                        core_complete=sum(all(r[s] for s in core_ok) for r in rs) / len(rs) if core_ok else None))
    phylum_n = {}
    for r in rows:
        phylum_n[r["PHYLUM"]] = phylum_n.get(r["PHYLUM"], 0) + 1
    out.sort(key=lambda d: (-phylum_n.get(d["lineage"]["PHYLUM"], 0), d["lineage"]["PHYLUM"], d["lineage"]["SUBPHYLUM"],
                            d["lineage"]["CLASS"], d["lineage"]["ORDER"], d["lineage"]["FAMILY"], d["group"]))
    return out


def draw(summary, steps_interp, steps_floor, out_dir, name, title, note):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_rows = len(summary)
    cols_a = steps_interp + ["core_complete"]
    fig_h = max(3.2, 0.26 * n_rows + 2.4)
    fig = plt.figure(figsize=(11.5, fig_h))
    ratios = [len(cols_a), max(len(steps_floor), 1) * 0.9]
    gs = fig.add_gridspec(1, 2, width_ratios=ratios, wspace=0.04)
    axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    labels = [f'{d["group"]}  (n={d["n"]:,})' for d in summary]

    def panel(ax, cols, cmap, header):
        mat = np.array([[ (d["core_complete"] if c == "core_complete" else d["frac"][c]) for c in cols] for d in summary], dtype=float)
        ax.pcolormesh(mat, cmap=cmap, vmin=0, vmax=1, edgecolors=SURFACE, linewidth=0.8)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                ax.text(j + 0.5, i + 0.5, f"{round(100 * v)}" if v >= 0.005 else "", ha="center", va="center", fontsize=6.3,
                        color="#ffffff" if v > 0.55 else INK)
        ax.set_xlim(0, mat.shape[1]); ax.set_ylim(mat.shape[0], 0)
        ax.set_xticks(np.arange(len(cols)) + 0.5)
        ax.set_xticklabels([c.replace("_", " ") if c != "core_complete" else "3 core carotenoid steps" for c in cols], rotation=60, ha="left", fontsize=7.5)
        ax.xaxis.tick_top()
        ax.tick_params(length=0)
        ax.set_title(header, fontsize=8, loc="left", pad=92, color=INK2)
        for sp in ax.spines.values():
            sp.set_visible(False)
    panel(axa, cols_a, BLUE_RAMP, "steps with a held-out cutoff (percent of genomes with a call)")
    axa.set_yticks(np.arange(n_rows) + 0.5)
    axa.set_yticklabels(labels, fontsize=7.5)
    if steps_floor:
        panel(axb, steps_floor, GREY_RAMP, "family-level steps: cutoff is a decoy floor, read as 'family present'")
        axb.set_yticks([])
    else:
        axb.axis("off")
    # phylum separators
    prev = None
    for i, d in enumerate(summary):
        ph = d["lineage"]["PHYLUM"]
        if prev is not None and ph != prev:
            for ax in (axa, axb):
                if ax.axison:
                    ax.axhline(i, color=INK2, lw=1.1)
        prev = ph
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, fontweight="bold", y=1.075)
    fig.text(0.01, 0.005, note, fontsize=7, color=INK2, ha="left", va="bottom", wrap=True)
    fig.savefig(out_dir / f"{name}.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", required=True, type=Path, help="results/genomes.tsv.zst (locustag, has_proteins, busco_complete)")
    ap.add_argument("--species", required=True, type=Path, help="species.parquet (taxonomy)")
    ap.add_argument("--search-dir", required=True, type=Path)
    ap.add_argument("--min-scores", required=True, type=Path)
    ap.add_argument("--placed", action="append", default=[], metavar="STEP=assigned.tsv")
    ap.add_argument("--rank", default="class", choices=list(RANKS))
    ap.add_argument("--within", default="", help="restrict to a lineage value, for example SUBPHYLUM=Pezizomycotina")
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--min-busco", type=float, default=0.0, help="keep genomes with BUSCO completeness >= this (0 keeps all)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--name", default="fig7_presence_by_class")
    a = ap.parse_args(argv)

    tax = load_taxonomy(a.species)
    genomes = [r for r in read_tsv(a.genomes) if r["has_proteins"].strip().lower() == "true"]
    n_searched = len(genomes)
    if a.min_busco > 0:
        genomes = [r for r in genomes if r["busco_complete"] not in ("", "<NA>", "NA") and float(r["busco_complete"]) >= a.min_busco]
    loc = [r["locustag"] for r in genomes]
    calls = hits_by_step(a.search_dir)
    placed_steps = []
    for item in a.placed:
        step, _, path = item.partition("=")
        calls[step] = placed_genomes(path)
        placed_steps.append(step)
    floors = floor_steps(a.min_scores)
    steps = [s for s in STEP_ORDER if s in {r["step"] for r in read_tsv(a.min_scores)}]
    rows = presence_table(loc, tax, calls, steps)
    if a.within:
        col, _, val = a.within.partition("=")
        rows = [r for r in rows if r.get(col) == val]
    summary = summarise(rows, RANKS[a.rank], steps, a.min_n)
    interp = [s for s in steps if s not in floors]
    flo = [s for s in steps if s in floors]
    where = f" within {a.within}" if a.within else ""
    title = f"Fraction of genomes with a call per step, by {a.rank}{where}"
    busco_txt = f"with BUSCO completeness >= {a.min_busco:g} " if a.min_busco > 0 else "(no BUSCO filter) "
    note = (f"{len(rows):,} genomes {busco_txt}of {n_searched:,} searched; groups with fewer than {a.min_n} genomes dropped "
            f"({len(summary)} groups shown). Placed calls: {', '.join(placed_steps) or 'none'}; every other step is a search hit at its cutoff, not a placed call.")
    draw(summary, interp, flo, a.out_dir, a.name, title, note)
    print(f"{len(summary)} groups, {len(rows)} genomes; wrote {a.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
