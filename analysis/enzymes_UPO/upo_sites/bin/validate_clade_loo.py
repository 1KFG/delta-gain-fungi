#!/usr/bin/env python3.12
"""
Leave-one-out validation of the upo_sites clade caller.

For each of the 19 calibration enzymes: remove it from the backbone, call its clade against
the remaining 18, and compare with the published assignment. This is the only honest estimate
of the caller's accuracy, because every other test would have the answer in the reference set.

Reports overall accuracy, per-clade accuracy, and the identity/margin distribution of correct
versus incorrect calls, so a confidence threshold can be chosen from evidence.
"""
from __future__ import annotations
import os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from upo_sites import (load_backbone, run_hmmalign, parse_stockholm,
                       match_column_indices, identity, build_position_map,
                       SITES, POCKET5, ALIPHATIC)

HMM = sys.argv[1]
BACKBONE = sys.argv[2]
ANCHOR = "HspUPO"

bb = load_backbone(BACKBONE)
print(f"backbone: {len(bb)} sequences")

workdir = tempfile.mkdtemp(prefix="loo.")
anchor_fa = os.path.join(workdir, "anchor.fa")
with open(anchor_fa, "w") as fh:
    for n, c, s in bb:
        if n == ANCHOR:
            fh.write(f">{n}\n{s}\n")
pmap = build_position_map(HMM, anchor_fa, ANCHOR, set(SITES))
pocket_use = [p for p in POCKET5 if p in pmap]

# one alignment of all 19 is enough: identities are computed over match columns
allfa = os.path.join(workdir, "all.fa")
with open(allfa, "w") as fh:
    for n, c, s in bb:
        fh.write(f">{n}\n{s}\n")
sto = os.path.join(workdir, "all.sto")
run_hmmalign(HMM, allfa, sto)
seqs, rf = parse_stockholm(sto)
mcols = match_column_indices(rf)
M = {}
for n, c, s in bb:
    k = next((k for k in seqs if k.split("|")[0] == n), None)
    M[n] = (c, "".join(seqs[k][i] for i in mcols).upper())

rows, correct = [], 0
for q, (qc, qm) in M.items():
    best = {}
    for r, (rc, rm) in M.items():
        if r == q:
            continue                      # leave one out
        v = identity(qm, rm)
        if rc not in best or v > best[rc][0]:
            best[rc] = (v, r)
    ranked = sorted(best.items(), key=lambda kv: -kv[1][0])
    call, (cid, cref) = ranked[0]
    margin = cid - ranked[1][1][0] if len(ranked) > 1 else cid
    ok = call == qc
    correct += ok
    rows.append((q, qc, call, cid, margin, cref, ok))

print(f"\n{'enzyme':14} {'true':5} {'call':5} {'ident':>6} {'margin':>7}  {'nearest':14} result")
for q, qc, call, cid, margin, cref, ok in sorted(rows, key=lambda r: (r[1], r[0])):
    print(f"{q:14} {qc:5} {call:5} {cid:6.1f} {margin:7.1f}  {cref:14} "
          f"{'OK' if ok else 'WRONG'}")
n = len(rows)
print(f"\nLeave-one-out accuracy: {correct}/{n} = {100*correct/n:.1f}%")

from collections import defaultdict
per = defaultdict(lambda: [0, 0])
for *_, qc, call, cid, margin, cref, ok in [(None,) + r for r in rows]:
    per[qc][1] += 1
    per[qc][0] += ok
print("\nper clade (correct/total):")
for c in sorted(per):
    a, b = per[c]
    print(f"  Clade {c:3} {a}/{b}")

cok = [r for r in rows if r[6]]
cbad = [r for r in rows if not r[6]]
def rng(v):
    return f"min={min(v):.1f} median={sorted(v)[len(v)//2]:.1f} max={max(v):.1f}"
if cok:
    print(f"\ncorrect calls   identity {rng([r[3] for r in cok])}")
    print(f"                margin   {rng([r[4] for r in cok])}")
if cbad:
    print(f"incorrect calls identity {rng([r[3] for r in cbad])}")
    print(f"                margin   {rng([r[4] for r in cbad])}")

import shutil
shutil.rmtree(workdir, ignore_errors=True)
