#!/usr/bin/env python3
"""Figures for a pathway census report: steps 2 to 4 (score separation, candidate counts, low-cutoff side run).

Reads only files written by the census tools (bin/score_refs_decoys.py, bin/make_min_scores.py, bin/search_chunk.py).
Writes each figure as PNG (150 dpi, for the page) and PDF (vector, linked from the caption).
Run: /usr/bin/python3.12 fig_census_search.py --pathway carotenoid --step2-dir ... --search-dir ... --min-scores ... --out-dir ...
Colours: the first three slots of the validated categorical palette (blue, orange, aqua). Aqua is below 3:1 contrast
on the light surface, so every series is also direct-labelled and has its own marker shape.
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

BLUE, ORANGE, AQUA, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#8a8985", "#0b0b0b", "#52514e"
SURFACE = "#fcfcfb"


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GREY, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6, "axes.axisbelow": True,
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "pdf.fonttype": 42,
    })


def read_tsv(path):
    path = str(path)
    if path.endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    else:
        raw = Path(path).read_text()
    lines = [l for l in raw.splitlines() if not l.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def step_table(step_scores, min_scores, steps=None):
    """One dict per step: self_min, loo_min, max_decoy, rule (min_score of the rule), used (min_score used in the run)."""
    ss = {r["step"]: r for r in read_tsv(step_scores)}
    out = []
    for r in read_tsv(min_scores):
        if steps and r["step"] not in steps:
            continue
        s = ss[r["step"]]
        out.append(dict(step=r["step"], n_refs=int(float(s["n_refs"])), self_min=num(s["ref_min"]),
                        loo_min=num(s.get("loo_min")), max_decoy=num(s.get("max_decoy_score")),
                        rule=num(r.get("rule_min_score")) or num(r["min_score"]), used=num(r["min_score"])))
    return out


def count_candidates(search_dir):
    """step -> (rows, genomes) from candidates.chunk*.tsv.zst in search_dir."""
    rows, genomes = {}, {}
    for f in sorted(glob.glob(str(Path(search_dir) / "candidates.chunk*.tsv.zst"))):
        for r in read_tsv(f):
            s = r["step_id"]
            rows[s] = rows.get(s, 0) + 1
            genomes.setdefault(s, set()).add(r["protein_id"].split("_")[0])
    return {s: (rows[s], len(genomes[s])) for s in rows}


def scores_for(search_dir, step):
    out = []
    for f in sorted(glob.glob(str(Path(search_dir) / "candidates.chunk*.tsv.zst"))):
        out += [float(r["score"]) for r in read_tsv(f) if r["step_id"] == step]
    return out


def save(fig, out_dir, name):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_separation(tab, out_dir, name="fig1_step2_separation"):
    """Per step: highest decoy score, lowest held-out reference score, lowest own-reference score."""
    tab = sorted(tab, key=lambda r: (r["loo_min"] if r["loo_min"] is not None else r["self_min"]))
    fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(tab) + 1.4))
    for i, r in enumerate(tab):
        pts = [x for x in (r["max_decoy"], r["loo_min"], r["self_min"]) if x is not None]
        ax.plot([min(pts), max(pts)], [i, i], color="#cfcec9", lw=1.5, zorder=1)
        if r["max_decoy"] is not None:
            ax.scatter(r["max_decoy"], i, marker="s", s=34, color=GREY, zorder=3, edgecolor=SURFACE, linewidth=1)
        if r["loo_min"] is not None:
            ax.scatter(r["loo_min"], i, marker="D", s=34, color=ORANGE, zorder=3, edgecolor=SURFACE, linewidth=1)
        if r["self_min"] is not None:
            ax.scatter(r["self_min"], i, marker="o", s=38, color=BLUE, zorder=3, edgecolor=SURFACE, linewidth=1)
        ax.scatter(r["used"], i, marker="|", s=170, color=INK, zorder=4, linewidth=1.6)
    ax.set_yticks(range(len(tab)))
    ax.set_yticklabels([f'{r["step"]} (n={r["n_refs"]})' for r in tab])
    ax.set_xscale("log")
    ax.set_xlabel("HMM bit score (log scale)")
    ax.grid(axis="y", visible=False)
    handles = [plt.Line2D([], [], marker="s", ls="", color=GREY, label="highest decoy"),
               plt.Line2D([], [], marker="D", ls="", color=ORANGE, label="lowest held-out reference"),
               plt.Line2D([], [], marker="o", ls="", color=BLUE, label="lowest own reference (trains the HMM)"),
               plt.Line2D([], [], marker="|", ls="", color=INK, markersize=12, markeredgewidth=1.6, label="cutoff used in the search")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13 if len(tab) > 8 else -0.42), ncol=2, frameon=False, fontsize=8)
    ax.set_title("Step 2: references against decoys on the three mycosporine step HMMs")
    save(fig, out_dir, name)


def fig_counts(counts, out_dir, name="fig2_candidates_per_step"):
    items = sorted(counts.items(), key=lambda kv: kv[1][0])
    fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(items) + 1.2))
    ax.barh(range(len(items)), [v[1][0] for v in items], color=BLUE, height=0.62)
    for i, (s, (n, g)) in enumerate(items):
        ax.text(n * 1.08, i, f"{n:,} rows, {g:,} genomes", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([s for s, _ in items])
    ax.set_xscale("log")
    ax.set_xlim(left=1, right=max(v[0] for _, v in items) * 12)
    ax.set_xlabel("candidate rows (protein by step), log scale")
    ax.grid(axis="y", visible=False)
    ax.set_title("Step 4: candidates at the cutoffs used, all searched genomes")
    save(fig, out_dir, name)


def fig_side(side_scores, cuts, out_dir, name="fig3_low_cutoff_side_run"):
    """side_scores: step -> list of scores (floor run). cuts: step -> (rule, used). One panel per step."""
    steps = list(side_scores)
    fig, axes = plt.subplots(1, len(steps), figsize=(3.8 * len(steps), 3.8), squeeze=False)
    for ax, s in zip(axes[0], steps):
        v = side_scores[s]
        hi = max(v + [cuts[s][0]]) * 1.05
        ax.hist(v, bins=[i * 25 for i in range(int(hi // 25) + 2)], color=BLUE, edgecolor=SURFACE, linewidth=0.6)
        ax.set_yscale("log")
        rule, used = cuts[s]
        ax.axvline(used, color=ORANGE, lw=1.6, label=f"cutoff used in the full search ({used:g})")
        ax.axvline(rule, color=INK, lw=1.6, ls="--", label=f"cutoff of the rule ({rule:g})")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), frameon=False, fontsize=7.5)
        ax.set_title(f"{s}\n{len(v):,} proteins in 1,000 genomes", fontsize=9)
        ax.set_xlabel("HMM bit score (25-bit bins)")
        ax.set_ylabel("proteins (log scale)")
    save(fig, out_dir, name)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step2-dir", required=True, type=Path)
    ap.add_argument("--search-dir", required=True, type=Path)
    ap.add_argument("--min-scores", required=True, type=Path, help="min_scores TSV used in the run (has rule_min_score)")
    ap.add_argument("--side-dir", type=Path, help="side evaluation folder (candidates + its min_scores TSV)")
    ap.add_argument("--steps", nargs="*", default=None, help="keep only these steps (the search folder and cutoff table hold 18 steps of three pathways)")
    ap.add_argument("--tag", default="PROVISIONAL")
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args(argv)
    style()
    sfx = f".{a.tag}" if a.tag else ""
    tab = step_table(a.step2_dir / f"step_scores{sfx}.tsv", a.min_scores, a.steps)
    fig_separation(tab, a.out_dir)
    counts = count_candidates(a.search_dir)
    if a.steps:
        counts = {s: v for s, v in counts.items() if s in a.steps}
    fig_counts(counts, a.out_dir)
    if a.side_dir:
        cuts = {r["step"]: (r["rule"], r["used"]) for r in tab}
        side_steps = sorted({r["step_id"] for f in glob.glob(str(a.side_dir / "candidates.chunk*.tsv.zst")) for r in read_tsv(f)})
        side = {s: scores_for(a.side_dir, s) for s in side_steps}
        fig_side(side, cuts, a.out_dir)
    print("wrote figures to", a.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
