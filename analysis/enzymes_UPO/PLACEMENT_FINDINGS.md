# Phylogenetic placement: a clade caller with real uncertainty, and what it says about the clades

**Measured 2026-11-08.** Replaces nearest-neighbour clade assignment, which `PLAN.md` §3.1 flagged
as validated only at 52.7–96.2 % identity while BFD sits at 25–35 %.
Code: `upo_sites/bin/make_placement_aln.py`, `jplace_to_clade.py`, `upo_sites/slurm/{epa_place,refbs}.sh`.

## Method

PF01328 **match-state columns only** (fixed by the profile, so reference and query alignments are
directly comparable), 93 references, 183 columns. Reference tree by RAxML-NG (LG+G4) with **400
bootstrap replicates**. Queries placed by RAxML-8 EPA. A query is assigned to a clade when the
placement edges inside that clade's subtree carry the majority of its likelihood weight; the
reported confidence is that weight, which is a genuine placement uncertainty rather than a percent
identity.

### Two reference-set defects found and fixed

1. **13 `NEG_` negative controls align to 1 of 183 match columns.** DyP, LiP, MnP, CCP, HRP, KatG,
   several P450s and VCPO are not PF01328, so `hmmalign` places almost nothing and they become
   all-gap tips that attach anywhere at ~zero branch length. Excluded. All 19 clade-labelled
   references retained; 93 references kept.
2. **Columns that are gaps in more than half the references are dropped**, because EPA behaves
   badly on near-empty columns.

## The caller is substantially more conservative, and that is the point

On the 7,055 dereplicated catalytic-core representatives:

| | nearest neighbour | EPA placement |
|---|---|---|
| Unassigned | 19.5 % | **57 %** |
| Confidence measure | % identity to nearest reference | likelihood weight ratio |

Where both methods commit, they **agree on 2,825 of 3,051 calls (92.6 %)**. Placement declines to
call 4,003 of 7,055. At 25–35 % identity that is the correct behaviour: the nearest-neighbour
caller was forcing calls it could not support, and had no way to express doubt.

Cross-tabulation (rows = nearest neighbour, columns = placement):

| | I | II | III | IV | V | unassigned |
|---|---|---|---|---|---|---|
| **I** | 701 | 0 | 0 | 0 | 0 | 1150 |
| **II** | 0 | 252 | 0 | 0 | 0 | 333 |
| **III** | 0 | 26 | 486 | 0 | 1 | 566 |
| **IV** | 0 | 0 | 0 | **652** | 17 | 120 |
| **V** | 0 | 170 | 5 | 7 | 734 | 1550 |

## A methodological error I made, and the correction

My first monophyly test computed the **rooted** MRCA of each clade. **RAxML trees are unrooted**,
and the newick rooting is arbitrary. That test reported Clade III as spanning the whole tree — i.e.
polyphyletic — which was an artifact of MorUPO and MroUPO sitting on opposite sides of the
arbitrary root, not a biological result.

The correct test on an unrooted tree is **bipartition-based**: a clade is monophyletic if some
branch separates all its members from every other labelled reference. Redone that way:

| Clade | refs | monophyletic? | bootstrap support |
|---|---|---|---|
| **II** | 2 | yes | **100** |
| **III** | 2 | yes | **100** |
| I | 5 | yes | **53** |
| IV | 1 | trivially (single reference) | n/a |
| **V** | 9 | **no** — MphUPO falls inside it | — |

**Clade III is strongly supported after all.** The polyphyly was mine, not the data's.

## What the tree says about Clade IV — unresolved, not refuted

MphUPO, the sole Clade IV reference, falls **within** the Clade V radiation: no bipartition
separates the nine Clade V references from MphUPO except MphUPO's own terminal branch.

But the backbone nodes that place it there are **not supported**. Bootstrap values along the path
from the root to MphUPO:

```
91 tips  support = 99
90 tips  support = 100
39 tips  support = 28     <- clades II, IV, V
22 tips  support = 21     <- clades IV, V
14 tips  support = 88
13 tips  support = 20     <- the node that makes Clade V non-monophyletic
 8 tips  support = 27
```

**The nodes that matter carry 20–28 % support.** So the relationship between Clade IV and Clade V
is **unresolved in this reference framework**. The data support neither "Clade IV is a distinct
clade" nor "Clade IV is nested inside Clade V".

### What this does and does not change for G1 Finding 3

**Unchanged:** the group is real. 652 of 789 nearest-neighbour Clade IV representatives are
confirmed by placement onto MphUPO's branch. It remains the tightest group in BFD (66.5 % median
identity against 33–41 % elsewhere), the most catalytically intact, and anchored on a single
uncharacterized enzyme. **Characterizing a member is still the top experimental recommendation**,
and arguably a stronger one now: the group's phylogenetic position cannot be settled without
biochemistry or more references.

**Changed:** calling it "the best-supported clade" overstated it. It is best-*delimited* by
sequence identity, but its status as a separate clade is unresolved at 20–28 % bootstrap support.
Any text claiming Clade IV is a distinct subtype should say "unresolved relative to Clade V".

**A statistic not to quote:** Clade IV placements all score weight 1.00 because its subtree is a
single edge, so anything landing there trivially takes full weight. That is an artifact of having
one reference, not strength of evidence, and it is not comparable to Clade V's 27 edges.

## The real limitation, restated

Three of the five clades rest on one or two reference enzymes (IV: 1, II: 2, III: 2). Clade I has
5 and only 53 % support. **Every unresolved question here traces to 19 characterized enzymes.**
Expanding the reference set (Zamocky 2022 / RedoxiBase) would do more for clade annotation than any
further method work on our side.
