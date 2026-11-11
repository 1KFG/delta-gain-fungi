#!/usr/bin/env python3.12
"""
Apply the cross-validated T1 pocket score to every PF01328 protein in BFD.

THE RULE
--------
Count how many of the 5 alignment-reliable aliphatic pocket positions (HspUPO numbering
p73, p77, p81, p183, p187) match the consensus of the ACTIVE calibration enzymes: I L L A I.
All five are small hydrophobic residues.

On the 19 calibration enzymes the rule is monotone:

  matches   n   mean activity   active (>10 %)
  0/5       5        0.7          0 / 5
  1/5       4       12.8          1 / 4
  2/5       1        2.7          0 / 1
  3/5       2       29.2          2 / 2
  4/5       4       52.9          4 / 4
  5/5       3       38.6          3 / 3

A cut at >= 3 gives 9/9 precision and 9/10 recall on the calibration set.

WHAT THE CROSS-VALIDATION SAYS (score_t1.py)
--------------------------------------------
Leave-one-out, refitting the consensus on the 18 training enzymes every time:
  LOO R2 = +0.496, LOO AUC = 0.956, permutation p = 0.0005.
It beats a clade-only model out of sample (LOO R2 +0.434).

WHAT THE CROSS-VALIDATION ALSO SAYS, AND THIS IS THE LIMIT
-----------------------------------------------------------
9 of the 10 active calibration enzymes are clade V. Remove clade V and the AUC falls to 0.556,
which is chance, on ONE remaining active enzyme. So:

  This score detects clade-V-like aliphatic pockets. It is NOT demonstrated to find active
  enzymes outside clade V. The single out-of-clade active in the calibration set (GmaUPO-II,
  clade I, 49.3 % conversion) carries only 1 of 5 matches and the score misses it.

So a high score is a clade V hypothesis by another route, and a low score is NOT evidence of
inactivity. GmaUPO-II is the proof of that second point.

WHY EMIT IT ANYWAY
------------------
The placement caller declines 78.9 % of BFD. The pocket residues are measured for every protein
whether or not placement commits. For those 76,913 proteins the score is the only ranking signal
available. It is emitted with its measured limits so a user can rank candidates for assay, which
is a different and much weaker claim than predicting activity.

Substrate scope: this is 4-propylguaiacol peroxygenase conversion, one substrate, one assay, one
publication. It says NOTHING about which rings an enzyme cleaves.
"""
from __future__ import annotations
import argparse, subprocess, sys
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--table", required=True, help="upo_sites_bfd_final.tsv.zst")
ap.add_argument("--consensus", default="ILLAI", help="active consensus at p73,p77,p81,p183,p187")
ap.add_argument("--cut", type=int, default=3)
ap.add_argument("--out", required=True)
ap.add_argument("--report", default=None)
a = ap.parse_args()

CONS = a.consensus.upper()
assert len(CONS) == 5, "consensus must be 5 residues"

txt = subprocess.run(["zstd", "-dcq", a.table], capture_output=True, text=True).stdout
hdr = None
out = []
head_comments = []
tally = Counter()
by_trust = Counter()
by_clade = Counter()
unassigned_hits = Counter()

for line in txt.splitlines():
    if line.startswith("#"):
        head_comments.append(line)
        continue
    f = line.split("\t")
    if hdr is None:
        hdr = f
        ix = {k: j for j, k in enumerate(hdr)}
        out.append(f + ["t1_pocket_matches", "t1_pocket_flag"])
        continue
    pr = f[ix["pocket_residues"]].strip().upper()
    if len(pr) != 5:
        m, flag = -1, "no_pocket"
    else:
        m = sum(pr[k] == CONS[k] for k in range(5))
        flag = "clade_V_like" if m >= a.cut else "below_cut"
    tally[m] += 1
    if m >= a.cut:
        by_trust[f[ix["clade_trust"]]] += 1
        by_clade[f[ix["clade"]]] += 1
        if f[ix["clade"]] == "unassigned":
            unassigned_hits[f[ix["PHYLUM"]]] += 1
    out.append(f + [str(m), flag])

n = len(out) - 1
note = [
    "# t1_pocket_matches: matches (of 5) to the ACTIVE consensus " + CONS +
    " at HspUPO p73,p77,p81,p183,p187.",
    f"# t1_pocket_flag: 'clade_V_like' at >= {a.cut} matches.",
    "# CROSS-VALIDATED on 19 enzymes: LOO R2 +0.496, LOO AUC 0.956, permutation p 0.0005.",
    "# LIMIT: 9 of 10 active calibration enzymes are clade V. With clade V removed the AUC is",
    "#   0.556 (chance) on ONE active. This detects clade-V-like pockets. It is NOT shown to find",
    "#   actives outside clade V, and it MISSES the one that exists (GmaUPO-II, 49.3 % active,",
    "#   1/5 matches). A low score is therefore NOT evidence of inactivity.",
    "# Substrate: 4-propylguaiacol conversion only. Says nothing about ring cleavage.",
]
buf = head_comments + note + ["\t".join(r) for r in out]
p = subprocess.Popen(["zstd", "-qf", "-o", a.out], stdin=subprocess.PIPE, text=True)
p.communicate("\n".join(buf) + "\n")
print(f"wrote {n} rows to {a.out}", file=sys.stderr)

rep = ["# T1 pocket score across BFD\n",
       f"Consensus `{CONS}` at HspUPO p73, p77, p81, p183, p187. "
       f"Cut at >= {a.cut} matches.\n",
       "## Score distribution\n", "| matches | n | % |", "|---|---|---|"]
for k in sorted(tally):
    lab = "no pocket" if k < 0 else f"{k}/5"
    rep.append(f"| {lab} | {tally[k]} | {100*tally[k]/n:.1f} |")
hits = sum(v for k, v in tally.items() if k >= a.cut)
rep += [f"\n**{hits} proteins ({100*hits/n:.1f} %) score >= {a.cut}.**\n",
        "## Where those hits sit, by clade call\n", "| clade | n |", "|---|---|"]
for k, v in by_clade.most_common():
    rep.append(f"| {k} | {v} |")
rep += ["\n## By trust tier\n", "| tier | n |", "|---|---|"]
for k, v in by_trust.most_common():
    rep.append(f"| `{k}` | {v} |")
rep += ["\n## Hits the placement caller DECLINED, by phylum\n",
        "These are the proteins the score exists for: no clade call, but a measured pocket.\n",
        "| phylum | n |", "|---|---|"]
for k, v in unassigned_hits.most_common(15):
    rep.append(f"| {k or '(none)'} | {v} |")
if a.report:
    open(a.report, "w").write("\n".join(rep) + "\n")
    print(f"wrote {a.report}", file=sys.stderr)
print("\n".join(rep), file=sys.stderr)
