#!/usr/bin/env python3
"""Clade calls for PF01328 proteins by phylogenetic placement on the in-house reference tree.

Replaces nearest-neighbour identity to 19 enzymes (upo_sites `clade`) with likelihood-weighted
placement (epa-ng) on the 187-leaf tree of `data/refset_zamocky2022/tree/`. Standard library only;
the alignment and placement steps are external (hmmalign, epa-ng) and driven by
`../slurm/place_*.sh`.

Subcommands
  match-cols  hmmalign A2M/AFA output -> FASTA in reference-column space (inserts dropped)
  prune       drop one leaf from the reference tree + MSA (leave-one-out), emit it as the query
  assign      jplace + clade table -> one row per query with clade and its accumulated LWR
  loo         summarise leave-one-out assign outputs against the leaves' own clades

Clade of a placement. Each epa-ng placement is an edge with a like-weight ratio (LWR). An edge is
*inside* clade C when one side of the edge's bipartition is a proper subset of C's leaves; it is
the *stem* of C when one side equals C (the query then diverged before C's common ancestor, so it
is not counted as C). A query's support for C is the sum of LWR over its edges inside C. The call
is the clade with the highest support if that support is >= --min-lwr, otherwise `unassigned`.
Both clade systems are reported: our Clades I-V and Zamocky 2022 Figure 1 clades (they nest).
The rule is root-independent: it uses bipartitions, so the jplace tree's arbitrary root does not
matter.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "refset", "bin"))
from tree_outline import Node, parse_newick  # noqa: E402


# ---------- FASTA ----------

def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0] if line[1:].strip() else "", []
            elif line:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def write_fasta(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i + 60] + "\n")


def match_columns(aligned: str) -> str:
    """hmmalign A2M/AFA row -> match states only: drop lowercase (insert residues) and '.'
    (insert gaps). Uppercase and '-' are the HMM match columns."""
    return "".join(c for c in aligned if c == "-" or c.isupper())


# ---------- tree ----------

def to_newick(n: Node) -> str:
    def rec(x):
        s = "(" + ",".join(rec(c) for c in x.children) + ")" if x.children else x.name
        if x.children and x.label:
            s += x.label
        if x.parent is not None:
            s += f":{x.length:.10g}"
        return s
    return rec(n) + ";"


def prune(root: Node, leaf: str) -> Node:
    """Remove a leaf; dissolve the parent it leaves with one child (lengths summed)."""
    target = next(l for l in root.leaves() if l.name == leaf)
    p = target.parent
    p.children.remove(target)
    if len(p.children) == 1:
        only = p.children[0]
        if p.parent is None:                 # p was the root: promote its remaining child
            if only.children:
                only.parent = None
                only.length = 0.0
                return only
        else:
            g = p.parent
            only.length += p.length
            only.parent = g
            g.children[g.children.index(p)] = only
    elif p.parent is None and len(p.children) == 2:
        # the trifurcating root lost a child: unroot again (epa-ng wants an unrooted tree)
        c = next(x for x in p.children if x.children)
        o = next(x for x in p.children if x is not c)
        o.length += c.length
        o.parent = c
        c.children.append(o)
        c.parent, c.length = None, 0.0
        return c
    return root


# ---------- placement -> clade ----------

def load_clades(path):
    """clade_members.tsv (name_clades.py) -> {system: {clade: frozenset(leaves)}}."""
    sets = {"our": defaultdict(set), "fig1": defaultdict(set)}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["our_clade"]:
                sets["our"][r["our_clade"]].add(r["label"])
            if r["fig1_clade"]:
                sets["fig1"][r["fig1_clade"]].add(r["label"])
    return {k: {c: frozenset(v) for c, v in d.items()} for k, d in sets.items()}


def edge_bipartitions(root: Node):
    """edge number -> (side, other side) leaf sets."""
    allleaves = frozenset(l.name for l in root.leaves())
    out = {}

    def walk(n):
        s = frozenset([n.name]) if not n.children else frozenset().union(*(walk(c) for c in n.children))
        if n.edge is not None:
            out[n.edge] = (s, allleaves - s)
        return s
    walk(root)
    return out


def edge_clade_map(bip, clades):
    """edge -> (clade containing it or None, clade it is the stem of or None), per system."""
    out = {}
    for e, (a, b) in bip.items():
        inside, stem = None, None
        best = None
        for c, members in clades.items():
            if a == members or b == members:
                stem = c
            elif a < members or b < members:
                if best is None or len(members) < best:   # nested clades: innermost wins
                    inside, best = c, len(members)
        out[e] = (inside, stem)
    return out


def assign(jplace_path, clade_path, min_lwr):
    jp = json.load(open(jplace_path))
    root = parse_newick(jp["tree"])
    bip = edge_bipartitions(root)
    present = frozenset(l.name for l in root.leaves())
    # Restrict clades to leaves in this tree: under leave-one-out the held-out leaf is absent,
    # and its clade's stem must still be recognised as the stem.
    clades = {k: {c: m & present for c, m in v.items() if m & present}
              for k, v in load_clades(clade_path).items()}
    maps = {k: edge_clade_map(bip, v) for k, v in clades.items()}
    f = jp["fields"]
    ie, il, ip = f.index("edge_num"), f.index("like_weight_ratio"), f.index("pendant_length")
    rows = []
    for pl in jp["placements"]:
        names = pl.get("n") or [x[0] for x in pl.get("nm", [])]
        ps = sorted(pl["p"], key=lambda x: -x[il])
        row = {"best_edge_lwr": f"{ps[0][il]:.3f}", "pendant_length": f"{ps[0][ip]:.4f}",
               "n_placements": len(ps)}
        for system in ("our", "fig1"):
            sup, stem = defaultdict(float), defaultdict(float)
            for p in ps:
                inside, st = maps[system][p[ie]]
                if inside:
                    sup[inside] += p[il]
                if st:
                    stem[st] += p[il]
            if sup:
                c, v = max(sup.items(), key=lambda x: x[1])
            else:
                c, v = "", 0.0
            row[f"{system}_clade"] = c if v >= min_lwr else "unassigned"
            row[f"{system}_lwr"] = f"{v:.3f}"
            row[f"{system}_best"] = c
            row[f"{system}_stem_lwr"] = (f"{max(stem.values()):.3f}" if stem else "0.000")
        for n in names:
            rows.append({"query": n, **row})
    return rows


COLS = ["query", "our_clade", "our_lwr", "our_best", "our_stem_lwr", "fig1_clade", "fig1_lwr",
        "fig1_best", "fig1_stem_lwr", "best_edge_lwr", "pendant_length", "n_placements"]


def write_rows(path, rows, cols=COLS):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------- leave-one-out summary ----------

def loo(assign_files, clade_path):
    truth = {}
    with open(clade_path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            truth[r["label"]] = r
    out = []
    for p in assign_files:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                t = truth[r["query"]]
                for system, col in (("our", "our_clade"), ("fig1", "fig1_clade")):
                    exp = t[col] or "unassigned"
                    r[f"{system}_expected"] = exp
                    r[f"{system}_ok"] = "yes" if r[f"{system}_clade"] == exp else "no"
                out.append(r)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("match-cols")
    m.add_argument("aligned")
    m.add_argument("--length", type=int, help="expected column count (fail if different)")
    m.add_argument("--out", required=True)
    pr = sub.add_parser("prune")
    pr.add_argument("--tree", required=True)
    pr.add_argument("--msa", required=True)
    pr.add_argument("--leaf", required=True)
    pr.add_argument("--out-tree", required=True)
    pr.add_argument("--out-msa", required=True)
    pr.add_argument("--out-query", required=True)
    a_ = sub.add_parser("assign")
    a_.add_argument("jplace")
    a_.add_argument("--clades", required=True, help="clade_members.tsv from name_clades.py")
    a_.add_argument("--min-lwr", type=float, default=0.5)
    a_.add_argument("--out", required=True)
    lo = sub.add_parser("loo")
    lo.add_argument("assigned", nargs="+")
    lo.add_argument("--clades", required=True)
    lo.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    if a.cmd == "match-cols":
        n = 0
        with open(a.out, "w") as fh:
            for name, seq in read_fasta(a.aligned):
                s = match_columns(seq)
                if a.length and len(s) != a.length:
                    sys.exit(f"ERROR: {name} has {len(s)} match columns, expected {a.length}")
                write_fasta(fh, name, s)
                n += 1
        print(f"match-cols: {n} sequences", file=sys.stderr)
    elif a.cmd == "prune":
        root = prune(parse_newick(open(a.tree).read()), a.leaf)
        open(a.out_tree, "w").write(to_newick(root) + "\n")
        with open(a.out_msa, "w") as fm, open(a.out_query, "w") as fq:
            for name, seq in read_fasta(a.msa):
                write_fasta(fq if name == a.leaf else fm, name, seq)
    elif a.cmd == "assign":
        rows = assign(a.jplace, a.clades, a.min_lwr)
        write_rows(a.out, rows)
        print(f"assign: {len(rows)} queries", file=sys.stderr)
    else:
        rows = loo(a.assigned, a.clades)
        cols = COLS + ["our_expected", "our_ok", "fig1_expected", "fig1_ok"]
        write_rows(a.out, rows, cols)
        for system in ("our", "fig1"):
            ok = sum(r[f"{system}_ok"] == "yes" for r in rows)
            print(f"loo {system}: {ok}/{len(rows)} correct", file=sys.stderr)


if __name__ == "__main__":
    main()
