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

---

## Matched-identity validation: the caller is reliable above ~55 % identity, and BFD mostly is not

**Measured 2026-11-08**, closing the `PLAN.md` §3.1 item 3 that had been outstanding. Job 31868181,
22 minutes. Raw results: `data/placement/matched_identity_validation.tsv`.

Ordinary leave-one-out gave 18/19, but only tested the 52.7–96.2 % identity regime. This test
removes, for each labelled reference, **every other reference above an identity ceiling**, rebuilds
the tree, and places it — forcing the caller to work at BFD-like divergence. Accuracy is reported
against identity to the **nearest surviving reference**, which is directly comparable to a BFD
query.

| Identity to nearest surviving reference | n | correct | declined | **wrong** | precision when it commits |
|---|---|---|---|---|---|
| 70–100 % | 14 | 13 | 1 | 0 | **100 %** |
| 55–70 % | 18 | 15 | 1 | 2 | **88 %** |
| 45–55 % | 9 | 3 | 1 | **5** | **38 %** |
| 35–45 % | 2 | 1 | 0 | **1** | 50 % |

**Precision collapses below about 55 % identity.** At 45–55 % it is 38 %, barely above the 20 %
expected from guessing among five clades.

**And the caller does not protect itself by declining.** Only 3 of 43 scorable cases were declined.
Below 55 % it commits and is wrong rather than abstaining. That is the opposite of the behaviour
the confidence tiers were meant to provide, and it means **a placement weight above 0.5 is not
sufficient evidence at low identity**.

**71 of 114 test cases were impossible** — removing a reference's close relatives removed every
other member of its clade, because three of the five clades have only one or two references. That
is not a flaw in the test; it is the central limitation of the project stated quantitatively.

### What this means for the BFD results

Identity to the nearest reference across the 7,055 representatives: **median 43.8 %**.

| | representatives |
|---|---|
| ≥ 70 % identity | 823 (11.7 %) |
| **≥ 55 % identity (the reliable zone)** | **2,505 (35.5 %)** |
| ≥ 45 % | 3,365 (47.7 %) |

Of the 3,052 representatives placement assigns, **2,088 (68.4 %) sit above the 55 % threshold** and
can be trusted at 88–100 % precision. The remaining **964 assigned calls sit below it and should be
treated as unreliable**, with expected precision of 38–50 %.

**So the defensible output is 2,088 reliably-assigned representatives, not 7,055 and not 3,052.**

Reliable-zone calls by clade: I 679, IV 602, III 383, V 381, II 43.

### One piece of good news for Finding 3

**602 of the 652 Clade IV representatives fall in the reliable zone.** The Clade IV group is one of
the few parts of this dataset where clade assignment is actually trustworthy, because its members
sit at high identity to MphUPO. Its *rank* remains unresolved (20–28 % bootstrap), but *membership*
is well supported. The experimental recommendation stands on the firmer of the two claims.

### The honest headline

**A reliable identification model for fungal UPO subtypes does not currently exist, and cannot be
built from 19 reference enzymes.** What exists is a caller that is reliable above 55 % identity,
knows where that boundary is, and can therefore say which 2,088 of 97,520 proteins it is entitled
to an opinion about. Expanding the reference set is the only route to more, and no further method
work on our side substitutes for it.

---

## Does expanding the reference set with high-confidence BFD members help? No — it makes it worse

**Measured 2026-11-08.** Job 31868609. Raw results:
`data/placement/expanded/matched_identity_validation_expanded.tsv`.

Every limitation above traces to having 19 characterized enzymes. The cheap shortcut is to promote
confident BFD members to references. This tests whether that works, so the project does not rest on
an assumption either way.

**Construction.** 60 extended references drawn from the ≥70 % identity zone where matched-identity
validation measured 100 % precision, with placement weight ≥0.9, clustered at 60 % identity for
diversity and capped at 25 per clade. Reference set goes from 93 to 153, labelled from 19 to 79
(Clade I 30, IV 19, V 16, III 11, II 3).

**Test.** Queries are the **original 19 characterized enzymes only**. Extended references are
removed by the same identity ceiling as everything else, so an extended reference cannot leak the
answer for the enzyme it was derived from.

| Identity to nearest surviving reference | original (19 refs) | expanded (79 labelled) |
|---|---|---|
| 70–100 % | **100 %** (13 correct, 0 wrong) | **58 %** (11 correct, **8 wrong**) |
| 55–70 % | 88 % | 88 % |
| 45–55 % | 38 % | **29 %** |
| 35–45 % | 50 % (n=2) | 75 % (n=5) |
| **overall precision when committing** | **80 %** (32 correct, 8 wrong) | **66 %** (31 correct, **16 wrong**) |

**Expansion doubles the error rate, and the damage is worst exactly where the original set was
perfect.** At ≥70 % identity the caller went from flawless to 58 %.

### The mechanism, which is the useful part

The expanded run has more cases in the high-identity bin (19 versus 14) precisely **because the
extended references supply close neighbours**. Those new close neighbours carry **inferred** clade
labels, not measured ones. So a query that previously had to reach a genuine characterized enzyme
now finds a nearer BFD sequence instead — and inherits whatever that sequence's label got wrong.
Errors propagate into the tree topology as well as the labels.

This is the circularity risk made quantitative. It was flagged before the run; it is now measured.

### What this settles

**Semi-supervised reference expansion is not a substitute for characterized enzymes — it is worse
than doing nothing.** The 19-enzyme ceiling is real and cannot be worked around computationally.
Only genuinely characterized sequences (Zamocky 2022 / RedoxiBase, or new biochemistry) will move
clade annotation forward.

The production caller therefore keeps the **original 93-reference, 19-labelled** set. The expanded
set is retained only as the evidence for this negative result and must not be used for calling.
