#!/usr/bin/env python3
"""Name the clades of the in-house reference tree: Zamocky 2022 Figure 1 clades, and our
backbone Clades I-V. Writes one row per leaf plus a per-clade report.

  name_clades.py --tree mfp.contree --anchors fig1_anchors.tsv --dedup dedup.tsv \
                 --out-members clade_members.tsv --out-report clade_report.tsv

Zamocky Figure 1 clades (`--anchors`): each is the most recent common ancestor (MRCA) of its
anchor leaves in the tree rooted on the `root` row's anchors, extended upward by the same rule as
our clades below (parent UFBoot >= --min-ufb, no other clade's anchors). Anchors come from the paper's text
(genera, phyla, named enzymes), not from this tree; the `evidence` column says which. Checks that
do not depend on the anchors are reported next to each clade: its UFBoot, its taxon composition,
whether it contains another Figure 1 clade's anchors (a conflict), and its mean length against
the average the paper states.

Our Clades I-V (from the `B__<name>__clade<N>` leaves): start at the MRCA of a clade's backbone
members (the member itself for a one-member clade), then move up while the parent has
UFBoot >= --min-ufb and contains no backbone member of another clade. A clade whose members' MRCA
already contains another clade's member is reported `not_monophyletic`.

A leaf can sit in no Figure 1 clade (deep backbone, unsupported) and in none of Clades I-V; the
columns are then empty. That is a result, not a gap to fill.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tree_outline import Node, parse_newick, reroot, ufb, group_of  # noqa: E402


def index(root: Node):
    """leaf name -> leaf node; node -> set of leaf names (by id)."""
    leaves, below = {}, {}

    def walk(n):
        if not n.children:
            leaves[n.name] = n
            below[id(n)] = frozenset([n.name])
            return below[id(n)]
        s = frozenset().union(*(walk(c) for c in n.children))
        below[id(n)] = s
        return s
    walk(root)
    return leaves, below


def mrca(leaves, names):
    names = list(names)
    missing = [x for x in names if x not in leaves]
    if missing:
        raise ValueError(f"anchors not in tree: {missing}")
    path = []
    n = leaves[names[0]]
    while n is not None:
        path.append(n)
        n = n.parent
    for n in path:
        if all(_under(leaves[x], n) for x in names):
            return n
    raise AssertionError("no common ancestor")


def _under(leaf, node):
    n = leaf
    while n is not None:
        if n is node:
            return True
        n = n.parent
    return False


def extendable(n: Node) -> bool:
    """True if n's parent may absorb n. Never into a child of the root: after rerooting, the
    whole ingroup carries the outgroup split's support, so extending there would swallow every
    non-outgroup leaf."""
    p = n.parent
    return p is not None and p.parent is not None and p.parent.parent is not None \
        and ufb(p) is not None


def load_anchors(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(
            (l for l in fh if not l.startswith("#")), delimiter="\t")]
    for r in rows:
        r["anchors"] = [a.strip() for a in r["anchors"].split(",") if a.strip()]
    return rows


def backbone_clades(leaves, below, min_ufb):
    members = {}
    for name in leaves:
        p = name.split("__")
        if p[0] == "B":
            members.setdefault(p[-1].replace("clade", ""), []).append(name)
    out = {}
    for clade, names in members.items():
        others = {n for c, ns in members.items() if c != clade for n in ns}
        node = mrca(leaves, names)
        if below[id(node)] & others:
            out[clade] = (None, "not_monophyletic", names)
            continue
        while (extendable(node) and ufb(node.parent) >= min_ufb
               and not (below[id(node.parent)] & others)):
            node = node.parent
        out[clade] = (node, "ok", names)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--dedup", required=True, help="dedup.tsv from build_refset.py (lengths)")
    ap.add_argument("--min-ufb", type=float, default=95)
    ap.add_argument("--out-members", required=True)
    ap.add_argument("--out-report", required=True)
    a = ap.parse_args(argv)

    anchors = load_anchors(a.anchors)
    root_row = next(r for r in anchors if r["clade"] == "root")
    root = reroot(parse_newick(open(a.tree).read()), set(root_row["anchors"]))
    leaves, below = index(root)
    with open(a.dedup, newline="") as fh:
        length = {r["label"]: int(r["length"]) for r in csv.DictReader(fh, delimiter="\t")}

    fig1 = [r for r in anchors if r["clade"] != "root"]
    anchor_owner = {x: r["clade"] for r in fig1 for x in r["anchors"]}
    nodes = {}
    for r in fig1:
        others = {x for x, c in anchor_owner.items() if c != r["clade"]}
        n = mrca(leaves, r["anchors"])
        while (extendable(n) and ufb(n.parent) >= a.min_ufb
               and not (below[id(n.parent)] & others)):
            n = n.parent
        nodes[r["clade"]] = n

    report = []
    for r in fig1:
        n = nodes[r["clade"]]
        lv = sorted(below[id(n)])
        conflict = sorted({anchor_owner[x] for x in lv
                           if x in anchor_owner and anchor_owner[x] != r["clade"]})
        mean = sum(length[x] for x in lv) / len(lv)
        report.append({
            "set": "zamocky_fig1", "clade": r["clade"], "n_leaves": len(lv),
            "ufboot": n.label or "root",
            "composition": ", ".join(f"{k} {v}" for k, v in
                                     Counter(group_of(x) for x in lv).most_common()),
            "mean_length": f"{mean:.0f}", "paper_mean_length": r.get("paper_mean_length", ""),
            "conflict_with": ";".join(conflict), "status": "conflict" if conflict else "ok",
            "evidence": r.get("evidence", "")})

    bb = backbone_clades(leaves, below, a.min_ufb)
    for clade in sorted(bb, key=lambda c: ["I", "II", "III", "IV", "V"].index(c)
                        if c in ["I", "II", "III", "IV", "V"] else 99):
        node, status, names = bb[clade]
        lv = sorted(below[id(node)]) if node is not None else []
        report.append({
            "set": "backbone", "clade": f"Clade {clade}", "n_leaves": len(lv),
            "ufboot": (node.label if node is not None else ""),
            "composition": ", ".join(f"{k} {v}" for k, v in
                                     Counter(group_of(x) for x in lv).most_common()),
            "mean_length": f"{sum(length[x] for x in lv) / len(lv):.0f}" if lv else "",
            "paper_mean_length": "", "conflict_with": "", "status": status,
            "evidence": f"{len(names)} backbone members"})

    def innermost(leaf, cands):
        """Of the clades containing leaf, the smallest (clades can nest, e.g. Clade V in Fam I)."""
        hit = [(len(below[id(n)]), name, n) for name, n in cands if n is not None
               and leaf in below[id(n)]]
        return min(hit, key=lambda x: x[0]) if hit else None

    members = []
    fig_c = list(nodes.items())
    bb_c = [(f"Clade {c}", v[0]) for c, v in bb.items()]
    for leaf in sorted(leaves):
        p = leaf.split("__")
        f = innermost(leaf, fig_c)
        b = innermost(leaf, bb_c)
        members.append({
            "label": leaf, "set": "backbone" if p[0] == "B" else "zamocky", "name": p[1],
            "group": group_of(leaf), "length": length.get(leaf, ""),
            "fig1_clade": f[1] if f else "", "fig1_ufboot": (f[2].label if f else ""),
            "our_clade": b[1] if b else "", "our_clade_ufboot": (b[2].label if b else "")})

    for path, rows in ((a.out_report, report), (a.out_members, members)):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    for r in report:
        print(f"{r['set']:13s} {r['clade']:28s} n={r['n_leaves']:<3} ufb={r['ufboot']:<9} "
              f"len={r['mean_length']:>4}/{r['paper_mean_length'] or '-':>3} {r['status']:16s} "
              f"{r['composition']}", file=sys.stderr)


if __name__ == "__main__":
    main()
