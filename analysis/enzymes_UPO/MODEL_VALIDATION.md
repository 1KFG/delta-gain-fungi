# What is validated, and to what degree

**2026-11-08.** This document exists because I had been reporting "the model" as one thing. It is
four, and they have very different evidence behind them. Two are solid, one is bounded, one failed.

| Component | Status | Evidence |
|---|---|---|
| **A. Active-site identification** | **Validated** | 95/95 exact agreement with the published structure-based matrix |
| **B. Catalytic-core classification** | **Validated** | Reproduces the published "Glu180 in all 19" claim at 93.7 % |
| **C. Clade assignment** | **Bounded** | Reliable only above ~55 % identity; BFD median is 43.8 % |
| **D. Activity prediction** | **Failed statistically; mechanistic route partially re-opened** | r = 0.18 (T2); partial r = 0.13 given clade (T1); channel-width ranking r = +0.93 at n = 4 |

---

## A. Active-site identification — validated

**What it does.** Given any fungal PF01328 protein, locate the active-site residues in HspUPO
numbering via `hmmalign` against the Pfam profile, using match-state columns only.

**Evidence.**

1. **Exact agreement with independent ground truth.** At the five alignment-reliable aliphatic
   pocket positions, across all 19 characterized enzymes, `hmmalign` agrees with Figure S2's
   3DM **structure-based** alignment at **95 of 95 residues (100 %)**. This is a genuine
   cross-method check: a profile HMM sequence alignment reproducing a structure-based one.
2. **Across all 18 scoreable positions**, agreement is 95.5 % (322/337), with disagreements
   concentrated at positions 80, 172, 175, 223 and 232 — which is how those were identified as
   unreliable and excluded.
3. **The tool detects its own blind spots.** Positions 223, 228 and 232 fall in insert columns of
   PF01328 and are reported as unscoreable at startup, rather than silently returning garbage.
4. **Specificity control (new, 2026-11-08).** Ranking all 186 PF01328 match states by conservation
   across 87 phylogenetically spread BFD proteins, the mapped catalytic positions land in the top
   tier:

   | Position | BFD conservation | rank of 186 |
   |---|---|---|
   | Cys39 (proximal) | **98.8 %** | **6** |
   | Glu109 | **94.0 %** | **12** |
   | Asp111 | **90.4 %** | **15** |
   | His110 | 49.4 % | 53 |
   | Glu180 | 48.8 % | 56 |

   Median conservation across all match states is 34.5 %. The mapping is landing on the most
   conserved columns in the domain, not arbitrary ones.

**Scale.** All 97,520 BFD PF01328 proteins annotated in **1 min 58 s**, 188 MB peak.

**Honest limit.** Ground truth is 19 enzymes. The test shows the *method* reproduces a
structure-based alignment, not that every call on a 25 %-identity BFD protein is right. The proxy
for that is alignability: pocket columns align at 94.7 % against a catalytic-site ceiling of
95.7 %, so these positions are placed about as reliably as the most conserved residues in the
family.

### A note the pocket score depends on

The same specificity control shows the **aliphatic pocket positions are among the *less* conserved
columns** — ranks 21 to 138 of 186, with position 73 below the median. That is expected for a
specificity-determining region rather than a catalytic one, and it means the pocket score measures
real variation. It also means the pocket is not a conserved family feature and should never be
described as one.

---

## B. Catalytic-core classification — validated

**What it does.** Calls a PF01328 protein UPO-like if it carries Cys39 + Glu109 + Asp111 + Glu180.

**Evidence.** The 2025 study states Glu180 is present in all 19 characterized UPOs. Among BFD
proteins whose clade call is inside the validated identity range, **93.7 % carry Glu180** — an
independent reproduction of the published claim on 12,818 proteins the paper never saw. The rate
falls to 39.3 % for extrapolated calls and 21.8 % for unassigned, which is informative about those
proteins rather than about the method.

**Result.** 34,561 of 97,520 (35.4 %) are UPO-like. This is the set all downstream work uses.

---

## C. Clade assignment — bounded, and the bound is measured

**Reliable above ~55 % identity to a reference; not below.** Matched-identity validation:

| Identity to nearest reference | precision when it commits |
|---|---|
| 70–100 % | **100 %** |
| 55–70 % | **88 %** |
| 45–55 % | **38 %** |

BFD representatives sit at **median 43.8 %**. Of 3,052 assigned, **2,088 clear 55 % and are
trustworthy**; 964 do not.

Two things were tried and did not help, both measured rather than assumed:

- **Reference expansion with high-confidence BFD members made it worse** — overall precision
  80 % → 66 %, and 100 % → 58 % in the high-identity bin. Inferred labels propagate.
- **Phylogenetic placement** is better calibrated than nearest neighbour (57 % unassigned vs
  19.5 %, 92.6 % agreement where both commit) but does not extend the reliable range.

**Root cause, quantified:** three of the five clades have one or two reference enzymes, and 71 of
114 matched-identity test cases were impossible because removing a reference's close relatives
removed its entire clade.

---

## D. Activity prediction — failed, and reported as such

- **T2, aromatic ring hydroxylation:** pocket vs naphthalene r = **+0.176** (n = 19). Only two
  interpretable negatives, and both are the entire membership of Clade II.
- **T1, peroxygenase vs peroxidase:** r = +0.748 raw, but **partial r = +0.129 given clade**, and
  +0.001 within non-Clade-V enzymes. The pocket is clade in disguise.
- **Ring-size rule, mechanistic route: PARTIALLY RE-OPENED, corrected 2026-11-08.** An earlier
  entry here said CAVER gives ~1 Å bottleneck radii and the route was closed. That used CAVER
  3.0.2, which has a **documented bug** placing the starting vertex "in too narrow cavity".
  **CAVER 3.0.3** gives MroUPO 5.19 Å, HspUPO 4.87 Å, AaeUPO 4.23 Å, CviUPO 3.71 Å — the published
  ordering (MroUPO wider than AaeUPO), and r = +0.933 against naphthalene across those four.
  Absolute values still do not match the published 7 / 11 Å, so Karich's absolute cutoff cannot be
  applied; a **relative ranking** can. n = 4. See `CHANNEL_FINDINGS.md`.

**CORRECTION, 2026-11-08. A score IS now emitted. The statement below was wrong.**

The line that stood here said no activity score would be emitted, because producing one would
rank 97,520 proteins by a feature measured to have no predictive power. That conclusion came
from an IN-SAMPLE partial correlation (partial r = +0.129 given clade). It was never
cross-validated. Refusing on the strength of an unvalidated statistic is the same error as
publishing on the strength of one.

Cross-validating it changed the answer. `upo_sites/bin/score_t1.py`, leave-one-out with the
consensus, clade means and coefficients all refit per fold:

| model | LOO R² | LOO AUC |
|---|---|---|
| null | −0.114 | 0.000 |
| clade only | +0.434 | 0.922 |
| **pocket only** | **+0.496** | **0.956** |

Permutation p = 0.0005. The pocket beats the clade label OUT OF SAMPLE. The earlier "pocket is
clade in disguise" claim does not survive cross-validation and is withdrawn.

**What does survive is the clade V confound, in a sharper form.** 9 of the 10 active enzymes are
clade V. Remove clade V and the pocket AUC falls to 0.556 — chance — on ONE remaining active.
The score detects clade-V-like pockets. It is not shown to find actives outside clade V, and it
misses the one that exists (GmaUPO-II, 49.3 % active, 1/5 matches). So a low score is NOT
evidence of inactivity.

The score is emitted in `data/bfd_full/upo_sites_bfd_scored.tsv.zst` with those limits recorded
in the file header, because the placement caller declines 78.9 % of BFD and for those 76,913
proteins it is the only ranking signal available. See `T1_SCORE_REPORT.md`.

---

## The honest summary

**The site-finding and scoring software is reliable and validated.** It locates active-site
residues at 100 % agreement with structure-based ground truth, correctly identifies its own
unscoreable positions, reproduces an independent published claim on 12,818 proteins, and runs the
whole of BFD in two minutes.

**Subtype assignment is reliable for 2,088 of 7,055 representatives and says which.**

**Activity prediction: partially possible, and bounded.** A cross-validated pocket score reaches
LOO AUC 0.956 overall but 0.556 once clade V is removed, so it ranks candidates rather than
predicting activity. 5,115 proteins carry both the catalytic core and a passing score. The
mechanistic route — ranking by CAVER 3.0.3 channel width — is a live lead at r = +0.93 on four
enzymes, enough to justify modelling clade representatives, not enough to ship a predictor.

The binding constraint throughout is 19 characterized enzymes. It was tested, not assumed.
