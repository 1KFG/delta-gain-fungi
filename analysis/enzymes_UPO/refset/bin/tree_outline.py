#!/usr/bin/env python3
"""Outline an IQ-TREE tree as nested, supported clades with their taxon composition.

Used to match clades of the in-house reference tree to Zamocky et al. 2022 Figure 1 by
membership. Standard library only.

Leaf labels follow build_refset.py:  B__<name>__clade<N>  or  Z__<abbrev>__<group>.
Internal labels are IQ-TREE's "SH-aLRT/UFBoot".

  tree_outline.py TREE --root-with Z__RdelHTP__Mucoro [...] [--min-ufb 95] [--min-size 3]

The tree is rooted on the branch above the smallest clade that contains every --root-with leaf
(searched over all rootings, so the outgroup need not be a clade in the input rooting). A clade is
printed when it has >= --min-size leaves and UFBoot >= --min-ufb; its unsupported descendants are
folded into it.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter


class Node:
    __slots__ = ("name", "length", "children", "parent", "label", "edge")

    def __init__(self, name="", length=0.0):
        self.name, self.length, self.children, self.parent, self.label = name, length, [], None, ""
        self.edge = None  # jplace edge number, from "{n}" after the branch length

    def leaves(self):
        if not self.children:
            return [self]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def parse_newick(s: str) -> Node:
    """Parse Newick; a name after ")" is the internal node's label (support). jplace edge
    numbers ("{n}" after a branch length) are kept on the node below the edge."""
    tokens = re.findall(r"\(|\)|,|:[^,();]+|[^,():;]+", s.strip().rstrip(";"))
    root = cur = Node()
    stack = []
    for t in tokens:
        if t == "(":
            child = Node()
            child.parent = cur
            cur.children.append(child)
            stack.append(cur)
            cur = child
        elif t == ",":
            parent = stack[-1]
            cur = Node()
            cur.parent = parent
            parent.children.append(cur)
        elif t == ")":
            cur = stack.pop()
        elif t.startswith(":"):
            m = re.match(r":([^{]+)(?:\{(\d+)\})?", t)
            cur.length = float(m.group(1))
            if m.group(2) is not None:
                cur.edge = int(m.group(2))
        elif cur.children:
            cur.label = t
        else:
            cur.name = t
    return root


def neighbours(n: Node):
    out = list(n.children)
    if n.parent is not None:
        out.append(n.parent)
    return out


def reroot(top: Node, target_leaves: set[str]) -> Node:
    """Return a new rooted tree whose root splits the smallest clade containing target_leaves
    (over all rootings) from the rest. Supports stay on the branch (node) they describe."""
    allnodes = []

    def walk(n):
        allnodes.append(n)
        for c in n.children:
            walk(c)
    walk(top)
    # Each undirected edge (child, parent) defines a bipartition; the child side is child.leaves().
    best = None
    total = {l.name for l in top.leaves()}
    for n in allnodes:
        if n.parent is None:
            continue
        side = {l.name for l in n.leaves()}
        for s, flip in ((side, False), (total - side, True)):
            if target_leaves <= s and (best is None or len(s) < best[0]):
                best = (len(s), n, flip)
    if best is None:
        raise ValueError("outgroup leaves not found in tree")
    _, edge_child, flip = best
    # Build the rerooted tree by traversing the undirected graph from a new root on that edge.
    a, b = edge_child, edge_child.parent

    def build(n, came_from, length, label):
        m = Node(n.name, length)
        m.label = label
        for nb in neighbours(n):
            if nb is came_from:
                continue
            if nb.parent is n:          # moving down an original edge
                child = build(nb, n, nb.length, nb.label)
            else:                       # moving up: the edge is n->parent; its support is n's
                child = build(nb, n, n.length, n.label)
            child.parent = m
            m.children.append(child)
        return m

    root = Node()
    half = edge_child.length / 2
    ra = build(a, b, half, a.label)
    rb = build(b, a, half, a.label)
    for c in (ra, rb):
        c.parent = root
        root.children.append(c)
    # Dissolve degree-2 nodes left where the old root was (n with one child, not the root).
    def clean(n):
        for i, c in enumerate(list(n.children)):
            while len(c.children) == 1:
                g = c.children[0]
                g.length += c.length
                g.parent = n
                c = g
            n.children[i] = c
            clean(c)
    clean(root)
    return root


def ufb(n: Node):
    if not n.label:
        return None
    try:
        return float(n.label.split("/")[-1])
    except ValueError:
        return None


def group_of(leaf: str) -> str:
    parts = leaf.split("__")
    if parts[0] == "B":
        return "backbone:" + parts[-1]
    return parts[-1]


def summarize(leaves):
    c = Counter(group_of(l.name) for l in leaves)
    return ", ".join(f"{k} {v}" for k, v in sorted(c.items(), key=lambda x: -x[1]))


def outline(n: Node, min_ufb: float, min_size: int, depth=0, out=None, leaves_depth=-1):
    out = [] if out is None else out
    for c in n.children:
        lv = c.leaves()
        sup = ufb(c)
        if c.children and len(lv) >= min_size and sup is not None and sup >= min_ufb:
            bb = [l.name.split("__")[1] for l in lv if l.name.startswith("B__")]
            out.append(f"{'  ' * depth}[{len(lv)} | {c.label}] {summarize(lv)}"
                       + (f"  BB: {' '.join(bb)}" if bb else ""))
            if depth <= leaves_depth:
                out.append(f"{'  ' * depth}  leaves: {' '.join(l.name for l in lv)}")
            outline(c, min_ufb, min_size, depth + 1, out, leaves_depth)
        elif c.children:
            outline(c, min_ufb, min_size, depth, out, leaves_depth)
    return out


def ladder(root: Node):
    """Walk from the root down the larger child; at each node print what splits off. This is the
    order in which clades branch from the backbone (Zamocky Fig 1 reads this way)."""
    out, n, step = [], root, 0
    while n.children:
        kids = sorted(n.children, key=lambda c: len(c.leaves()))
        big = kids[-1]
        for k in kids[:-1]:
            lv = k.leaves()
            out.append(f"{step:2d} [{n.label or 'root'}] {len(lv):3d} leaves: {summarize(lv)}")
        n, step = big, step + 1
    out.append(f"{step:2d} tip: {n.name}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("tree")
    ap.add_argument("--root-with", nargs="+", required=True)
    ap.add_argument("--min-ufb", type=float, default=95)
    ap.add_argument("--min-size", type=int, default=3)
    ap.add_argument("--leaves-depth", type=int, default=-1,
                    help="list the leaves of printed clades down to this nesting depth")
    ap.add_argument("--ladder", action="store_true", help="print the backbone ladder instead")
    a = ap.parse_args(argv)
    top = parse_newick(open(a.tree).read())
    root = reroot(top, set(a.root_with))
    print(f"root: {len(root.children[0].leaves())} | {len(root.children[1].leaves())} leaves")
    if a.ladder:
        print("\n".join(ladder(root)))
        return
    for line in outline(root, a.min_ufb, a.min_size, leaves_depth=a.leaves_depth):
        print(line)


if __name__ == "__main__":
    main()
