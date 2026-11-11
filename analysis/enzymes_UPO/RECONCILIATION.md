# Reconciliation: two UPO workstreams, one main

**Drafted 2026-11-11** (laptop session `upo-iqtree-naming`) on branch `upo-reconcile`, which merges
the 40 commits of `upo-review-and-coordination` (the skylark session, `fungal-tf-clade-contrasts`,
work of 2026-11-08 after PR #1 merged) into `main` (which already carries PR #2, the Zamocky
reference set and the 187-leaf placement caller).

The two sessions worked in parallel on 2026-11-08/09 without seeing each other's later results.
Most of their work does not overlap. Where it does, both sessions reached the same diagnosis (19
characterized enzymes cannot support clade calls at BFD divergence) and took different routes out.
This document says which claims stand, which are superseded, and which wait on one measurement.

## What merges without contention

These skylark results do not depend on the clade caller and are carried onto main as written:

| Result | Document |
|---|---|
| CAVER 3.0.3 channel widths for 7 crystal structures; 3.0.2 bug and retraction; width–naphthalene lead r = +0.93 (n = 4) collapsing to +0.45 with TteUPO; Boltz-2 predicted structures fail channel validation | `CHANNEL_FINDINGS.md` |
| Novelty is deposition: GCF 2.0 % vs GCA-only 18.7 % novel, z = 7.4 | `G1_FINDINGS.md` §5b |
| Contamination check: UPO genes sit on more gene-rich contigs than random genes (lichens, Glomeromycota) | `G1_FINDINGS.md` §6 |
| G3 Phase A (Clade IV motifs incl. CCAAT box) and Phase B (RNA-seq, negative and interpretable) | `G3_*.md` |
| Site identification validated (95/95; conservation ranks 6/12/15 of 186) | `MODEL_VALIDATION.md` §A |
| T1 pocket ranking, clade-V-like (LOO AUC 0.956; out-of-clade 0.556); the reversal of the "no score" stance by cross-validation | `MODEL_CARD.md`, `T1_SCORE_REPORT.md` |
| External validation blocked by label incommensurability; 6-enzyme pre-registered experiment; 96-protein assay panel | `EXTERNAL_VALIDATION.md`, `ASSAY_PANEL.md` |
| `upo_scan.py` (proteome to scored table), coverage check | `upo_sites/README.md`, `data/bfd_full/COVERAGE_CHECK.md` |

Textual conflicts were only in `G1_FINDINGS.md` and `PLAN.md`. `G1_FINDINGS.md` keeps the
placement-based rewrite from main and adds skylark §5b and §6; `PLAN.md` keeps both sessions' risk
rows with statuses updated to what is now measured.

## 1. Two placement references: the one open measurement

| | skylark (`upo_sites`, 2026-11-08) | main (`upo_place`, 2026-11-09) |
|---|---|---|
| References | 93 tips, 19 labelled (characterized enzymes) | 187 tips, all labelled: 19 characterized + 168 Zamocky 2022 curated HTPs |
| Labels for non-characterized tips | none | from the in-house tree: clade = largest UFBoot >= 95 clade around the backbone members |
| Columns | 183 PF01328 match states | 224 MAFFT L-INS-i columns after trimAl -gt 0.75 |
| Tree / placement | RAxML-NG LG+G4, 400 bootstraps; RAxML-8 EPA | IQ-TREE ModelFinder topology, LG+G4 branch lengths; epa-ng |
| Plain leave-one-out | 18/19 | 184/187 (characterized 19: see below) |
| **Matched identity** (drop every reference above a ceiling) | **100 % precision at 70-100 %, 88 % at 55-70 %, 38 % at 45-55 %**; 71/114 cases impossible | **pending** (job `place_matched.sh`, same design) |
| BFD unassigned | 57 % of 7,055 representatives; 78.9 % of 97,520 in the final table | 42.0 % of 97,520 |

**Why these disagree on coverage, and why plain leave-one-out cannot settle it.** Main's
leave-one-out removed only the held-out leaf, so a query keeps its close relatives in the tree.
Its accuracy table was stratified by identity to the nearest of the **19 characterized** enzymes,
not to the nearest **surviving reference**. Skylark's matched-identity test is the harder and
correct comparison: it is what a BFD protein faces. So main's "32/32 below 30 % identity" is not
evidence against skylark's "38 % precision at 45-55 %". They measure different things.

**Why the 187-leaf reference might still do better at matched identity.** Skylark found 71 of 114
matched-identity cases impossible because removing a query's relatives removed its whole clade
(three clades had one or two references). The 187-leaf tree has 31/7/3/20/14 leaves in Clades I-V,
so far fewer cases lose their clade, and a query at 45 % identity to its nearest surviving
reference may still have several clade members at that distance. Whether that buys precision is an
empirical question.

**Why it might not.** Skylark showed that adding references with *inferred* labels degraded
precision from 80 % to 66 %, because wrong labels propagate. Zamocky's 168 sequences are
independently curated, but their *clade labels* here are inferred from the in-house tree. They are
not the same risk (labels come from a supported phylogeny, not from the caller being tested), but
they are not label-free either. The characterized-19 subset of the matched-identity test, whose
expected clades are published, is the independent check.

**Decision rule, written before the result:**

- If the 187-leaf reference, on the characterized 19 at matched identity, reaches **>= 80 %
  precision at 45-55 %** (skylark: 38 %) and does not lose precision at >= 55 %, main's placement
  calls become the production clade column, and the "unreliable below 55 %" statement is revised
  to the new measured boundary.
- If it does not, skylark's boundary stands for *both* references: `upo_place_bfd_all.tsv.zst`
  calls below 55 % identity to the nearest surviving reference are marked unverified, and the
  counts in `G1_BFD_SUMMARY.md` and `G1_FINDINGS.md` that rest on them get the same caveat.
- Either way the per-protein table keeps both callers' columns so the comparison stays visible.

## 2. Clade IV: membership agreed, rank now marginal rather than unresolved

| | skylark (93 refs) | main (187 refs) |
|---|---|---|
| Clade IV references | 1 (MphUPO) | 20 (MphUPO + 19 curated ascomycete HTPs) |
| IV vs V | MphUPO inside the Clade V radiation on nodes with **20-28 %** bootstrap | IV and V each recovered as clades |
| Support | standard bootstrap | Clade IV SH-aLRT/UFBoot **83.2/96** (LG+I+G4: 85.1/97); Clade V **93.8/86**; IV+V **98.5/99** |
| Membership | 652/789 NN Clade IV representatives confirmed | 8,437 BFD proteins; 5,316 were old IV calls, 2,853 old V |

**What stands:** both agree the Clade IV *group* is real and is anchored on one characterized
enzyme, so characterizing a member is the top experimental recommendation. **What changes:**
"unresolved" becomes "supported but marginal". Clade IV meets the usual joint threshold
(SH-aLRT >= 80 and UFBoot >= 95) only just, and Clade V misses it on UFBoot. The supports are not
directly comparable to skylark's standard bootstrap (UFBoot runs higher). Write "Clade IV, sister
to Clade V" and not "a distinct subtype"; `G1_FINDINGS.md` Finding 3 and `PLAN.md` say this.

The 2,853 proteins moving from old Clade V to IV depend on the IV/V boundary and are the calls
most sensitive to the matched-identity result in §1.

## 3. Clade I and Glu180: main's reading stands

Skylark's documents (`PLAN.md` row, earlier `G1_FINDINGS.md`) carried "the Clade I bin absorbs
non-UPOs" as **Confirmed**. Main overturned it on 2026-11-09, on two pieces of evidence: the
Glu180 deficit in placement Clade I is concentrated in ascomycetes (1.5 % vs 66.9 % in
basidiomycetes), and the 13 curated ascomycete Family II sequences, which are independent of BFD
and of either caller, all carry A/S/T at that position. Placement Clade I contains 99.8 % of the
old Clade I calls, so the phylum split very nearly describes the old bin as well; a direct count
on the old calls (one query on `upo_sites_bfd_all.tsv.zst`) would make that exact and is listed
in §5. The overturn does not depend on the §1 outcome, because Clade I precision is not in
question in either caller (31/31 leave-one-out; 100 % in skylark's >= 55 % zone). Open for both: whether ascomycete long UPOs lack the acid-base Glu or carry it at a
position that does not align to HspUPO 180.

## 4. Per-protein tables on the shared project

| File | Clade column from | Status |
|---|---|---|
| `data/bfd_full/upo_sites_bfd_final.tsv.zst`, `..._scored.tsv.zst` | skylark placement (`clade_trust`: reliable / extrapolated / unassigned) | carries the T1 pocket score; the skylark deliverable |
| `data/bfd_full/upo_place_bfd_all.tsv.zst` | main placement (`place_our_clade`, `place_fig1_clade`, LWR) | adds Zamocky Figure 1 lineages |

Next step after §1: one merged table with both clade columns plus the T1 score, and the decision
rule's outcome recorded in `STATUS.md` and `MODEL_CARD.md`.

## 5. Documents that need a follow-up edit once §1 is measured

- `STATUS.md`, `MODEL_CARD.md`: "No clade call for 78.9 % of BFD" and "precision 38 % below 55 %"
  are true of the skylark caller; add the 187-leaf result whichever way it lands.
- `PLACEMENT_FINDINGS.md`: pointer added; its own measurements stand.
- Glu180 by phylum on the *old* Clade I calls, to make §3 exact.
- `G1_BFD_SUMMARY.md`, `G1_FINDINGS.md`, `HANDOFF.md`: keep or caveat main's counts per the
  decision rule.
