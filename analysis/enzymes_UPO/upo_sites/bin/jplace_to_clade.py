#!/usr/bin/env python3.12
"""
Turn EPA placements into clade assignments with a real uncertainty measure.

Rule: for each labelled clade, take the smallest subtree containing all its reference members
(its MRCA subtree). A query is assigned to clade C if its placement edges inside C's MRCA subtree
carry a majority of the query's total likelihood weight. The reported confidence is that summed
likelihood weight, which is a genuine placement uncertainty -- unlike percent identity to the
nearest reference, which is what the nearest-neighbour caller used.

Clades whose MRCA subtree contains reference members of ANOTHER clade are flagged non-monophyletic
and reported as such rather than silently assigned.
"""
import argparse, json, re, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--jplace", required=True)
ap.add_argument("--labels", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min-weight", type=float, default=0.5)
a = ap.parse_args()

J = json.load(open(a.jplace))
tree = J["tree"]

# ---- parse newick carrying {edge} tags -------------------------------------
class N:
    __slots__ = ("ch", "name", "edge")
    def __init__(self): self.ch = []; self.name = None; self.edge = None

def parse(s):
    pos = 0
    def node():
        nonlocal pos
        n = N()
        if s[pos] == "(":
            pos += 1
            while True:
                n.ch.append(node())
                if s[pos] == ",": pos += 1
                elif s[pos] == ")": pos += 1; break
        m = re.match(r"([^,():;{}]*)", s[pos:]); pos += m.end()
        nm = m.group(1).strip()
        if nm: n.name = nm
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":[^,(){};]*", s[pos:]); pos += m.end()
        if pos < len(s) and s[pos] == "{":
            m = re.match(r"\{(\d+)\}", s[pos:]); n.edge = int(m.group(1)); pos += m.end()
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":[^,(){};]*", s[pos:]); pos += m.end()
        return n
    return node()

root = parse(tree)

tips_of = {}
edges_of = {}
def collect(n):
    t, e = set(), set()
    if n.edge is not None: e.add(n.edge)
    if not n.ch and n.name: t.add(n.name)
    for c in n.ch:
        ct, ce = collect(c); t |= ct; e |= ce
    tips_of[id(n)] = t; edges_of[id(n)] = e
    return t, e
collect(root)

lab = {}
with open(a.labels) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2: lab[f[0]] = f[1]

nodes = []
def walk(n):
    nodes.append(n)
    for c in n.ch: walk(c)
walk(root)

clade_edges, mono = {}, {}
for cl in sorted(set(lab.values())):
    mem = [k for k, v in lab.items() if v == cl]
    best = None
    for nd in nodes:
        if all(m in tips_of[id(nd)] for m in mem):
            if best is None or len(tips_of[id(nd)]) < len(tips_of[id(best)]):
                best = nd
    foreign = [k for k in tips_of[id(best)] if k in lab and lab[k] != cl]
    mono[cl] = not foreign
    clade_edges[cl] = edges_of[id(best)]
    print(f"clade {cl}: {len(mem)} refs, MRCA subtree {len(tips_of[id(best)])} tips, "
          f"{len(clade_edges[cl])} edges, "
          f"{'monophyletic' if mono[cl] else 'NOT monophyletic, foreign: ' + ','.join(foreign)}",
          file=sys.stderr)

fields = J["fields"]
i_edge, i_lwr = fields.index("edge_num"), fields.index("like_weight_ratio")

out = open(a.out, "w")
out.write("query\tclade\tplacement_weight\tn_placements\tbest_edge\tnote\n")
counts = defaultdict(int)
for pq in J["placements"]:
    names = [n[0] if isinstance(n, list) else n for n in pq.get("nm", pq.get("n", []))]
    w = defaultdict(float); tot = 0.0; best = (None, -1)
    for p in pq["p"]:
        e, lwr = p[i_edge], p[i_lwr]; tot += lwr
        if lwr > best[1]: best = (e, lwr)
        for cl, es in clade_edges.items():
            if e in es: w[cl] += lwr
    for nm in names:
        if w:
            cl = max(w, key=w.get); wt = w[cl] / max(tot, 1e-12)
        else:
            cl, wt = "unassigned", 0.0
        if wt < a.min_weight:
            cl, note = "unassigned", f"best={cl}:{wt:.2f}<{a.min_weight}"
        else:
            note = "" if mono.get(cl, True) else "clade_not_monophyletic"
        counts[cl] += 1
        out.write(f"{nm}\t{cl}\t{wt:.3f}\t{len(pq['p'])}\t{best[0]}\t{note}\n")
out.close()
print("\nplacement-based clade counts:", file=sys.stderr)
for k, v in sorted(counts.items(), key=lambda x: -x[1]):
    print(f"  {k:12} {v}", file=sys.stderr)
