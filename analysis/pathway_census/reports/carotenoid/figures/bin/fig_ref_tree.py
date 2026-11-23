#!/usr/bin/env python3
"""Reference tree of a step 5 placement run: references, decoys and candidate representatives (PNG and PDF).

Reads the files written by bin/step5_place.py: ref.rooted.nwk, clade_members.tsv, clusters.tsv.zst.
A representative's marker area grows with the log of the number of candidate domains in its mmseqs cluster.
Blue circles are references, orange squares are decoys, grey circles are representatives. Labels name every reference and decoy
(direct labels, so colour is never the only code). The coloured bar at the right marks the `target` clade.
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import Phylo  # noqa: E402

BLUE, ORANGE, AQUA, GREY, INK, INK2, SURFACE = "#2a78d6", "#eb6834", "#1baf7a", "#8a8985", "#0b0b0b", "#52514e", "#fcfcfb"


def read_cluster_sizes(path):
    """rep id -> number of members, from a two-column TSV (rep, member), plain or .zst."""
    path = str(path)
    if path.endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    else:
        raw = Path(path).read_text()
    sizes = {}
    for line in raw.splitlines():
        if line.strip():
            rep = line.split("\t")[0]
            sizes[rep] = sizes.get(rep, 0) + 1
    return sizes


def read_roles(path):
    """leaf id -> (clade, role) from clade_members.tsv."""
    with open(path) as fh:
        return {r["label"]: (r["our_clade"], r["role"]) for r in csv.DictReader(fh, delimiter="\t")}


def read_subclades(path):
    """leaf id -> lineage subclade name from the column fig1_clade of clade_members.tsv (empty when none)."""
    with open(path) as fh:
        return {r["label"]: (r.get("fig1_clade") or "") for r in csv.DictReader(fh, delimiter="\t")}


def layout(tree):
    """Return (xy, edges, leaves): x = root-to-node branch length, y = leaf order (0 at the top), internal y = mean of children."""
    leaves = tree.get_terminals()
    ypos = {id(t): i for i, t in enumerate(leaves)}
    xy, edges = {}, []

    def walk(clade, depth):
        x = depth + (clade.branch_length or 0.0)
        if clade.is_terminal():
            y = ypos[id(clade)]
        else:
            ys = [walk(c, x) for c in clade.clades]
            y = sum(ys) / len(ys)
        xy[id(clade)] = (x, y)
        for c in clade.clades:
            edges.append((id(clade), id(c)))
        return y

    root_depth = 0.0
    tree.root.branch_length = 0.0
    walk(tree.root, root_depth)
    return xy, edges, leaves


def draw(tree_path, members_path, clusters_path, out_dir, name="fig4_reference_tree", title="Reference tree used for placement"):
    tree = Phylo.read(str(tree_path), "newick")
    roles = read_roles(members_path)
    subclades = read_subclades(members_path)
    sizes = read_cluster_sizes(clusters_path)
    xy, edges, leaves = layout(tree)
    fig, ax = plt.subplots(figsize=(7.2, 9.0))
    for a, b in edges:
        (xa, ya), (xb, yb) = xy[a], xy[b]
        ax.plot([xa, xa], [ya, yb], color="#b9b8b3", lw=0.6, solid_capstyle="butt")
        ax.plot([xa, xb], [yb, yb], color="#b9b8b3", lw=0.6, solid_capstyle="butt")
    xmax = max(x for x, _ in xy.values())
    decoy_ys, decoy_names = [], []
    outgroup_ys, outgroup_names = [], []
    for leaf in leaves:
        x, y = xy[id(leaf)]
        role = roles.get(leaf.name, ("", ""))[1]
        if role == "reference":
            ax.scatter(x, y, marker="o", s=46, color=BLUE, zorder=3, edgecolor=SURFACE, linewidth=0.8)
            ax.text(x + xmax * 0.012, y, leaf.name, va="center", fontsize=6.5, color=INK)
        elif role == "named_member":
            ax.scatter(x, y, marker="D", s=40, color=AQUA, zorder=3, edgecolor=SURFACE, linewidth=0.8)
            ax.text(x + xmax * 0.012, y, leaf.name, va="center", fontsize=6, color=INK)
        elif role in ("tree_outgroup", "outgroup"):
            ax.scatter(x, y, marker="^", s=42, color=INK, zorder=3, edgecolor=SURFACE, linewidth=0.8)
            outgroup_ys.append(y)
            outgroup_names.append(leaf.name)
        elif role == "decoy":
            ax.scatter(x, y, marker="s", s=40, color=ORANGE, zorder=3, edgecolor=SURFACE, linewidth=0.8)
            decoy_ys.append(y)
            decoy_names.append(leaf.name)
        else:
            ax.scatter(x, y, marker="o", s=6 + 5 * math.log10(max(sizes.get(leaf.name, 1), 1)) ** 2 * 4,
                       color=GREY, zorder=2, edgecolor=SURFACE, linewidth=0.4)
    block = []
    if decoy_names:  # one stacked label for decoys and outgroups, right of the tree, so the names do not overprint each other
        block.append("decoys:\n" + "\n".join(decoy_names))
    if outgroup_names:
        block.append("outgroups:\n" + "\n".join(outgroup_names))
    if block:
        ax.text(xmax * 1.03, -0.3, "\n".join(block), va="top", fontsize=5.6, color=INK, linespacing=1.1)
    # bars for lineage subclades (named in column fig1_clade), left of the target bar
    groups = {}
    for l in leaves:
        g = subclades.get(l.name, "")
        if g:
            groups.setdefault(g, []).append(xy[id(l)][1])
    for g, yy in sorted(groups.items(), key=lambda kv: min(kv[1])):
        ax.plot([xmax * 1.12, xmax * 1.12], [min(yy), max(yy)], color=INK2, lw=3, solid_capstyle="butt")
        ax.text(xmax * 1.135, (min(yy) + max(yy)) / 2, f"{g} ({len(yy)})", va="center", fontsize=7, color=INK2)
    # bar for the target clade (leaves whose clade is "target", contiguous in this drawing if the clade is monophyletic)
    ys = [xy[id(l)][1] for l in leaves if roles.get(l.name, ("", ""))[0] == "target"]
    if ys:
        ax.plot([xmax * 1.30, xmax * 1.30], [min(ys), max(ys)], color=BLUE, lw=3, solid_capstyle="butt")
        ax.text(xmax * 1.32, (min(ys) + max(ys)) / 2, f"target clade\n{len(ys)} leaves", va="center", fontsize=8, color=INK)
    ax.set_ylim(len(leaves) - 0.5, -0.5)
    ax.set_xlim(-xmax * 0.01, xmax * 1.65)
    ax.set_yticks([])
    ax.set_xlabel("substitutions per site (LG+G4); rooted on the " + ("outgroups" if outgroup_names else "decoys"))
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    handles = [plt.Line2D([], [], marker="o", ls="", color=BLUE, label="reference"),
               plt.Line2D([], [], marker="D", ls="", color=AQUA, label="named member (tree leaf only)"),
               plt.Line2D([], [], marker="^", ls="", color=INK, label="outgroup"),
               plt.Line2D([], [], marker="s", ls="", color=ORANGE, label="decoy"),
               plt.Line2D([], [], marker="o", ls="", color=GREY, label="candidate representative (area ~ cluster size)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=1, frameon=False, fontsize=7.5)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step5-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--name", default="fig4_reference_tree")
    ap.add_argument("--title", default="Reference tree used for placement")
    a = ap.parse_args(argv)
    d = a.step5_dir
    draw(d / "ref.rooted.nwk", d / "clade_members.tsv", d / "clusters.tsv.zst", a.out_dir, a.name, a.title)
    print("wrote tree figure to", a.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
