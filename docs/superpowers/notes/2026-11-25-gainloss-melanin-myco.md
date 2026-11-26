---
title: "Gain and loss of the placed melanin and mycosporine steps on the genus tree"
subtitle: "DHN melanin, L-DOPA melanin and mycosporine: Mk models, 0.8 rule, tree-aware simulation null, Pagel"
date: 2026-11-25
---

::: {.callout-important}
## Status and limits

- Labels: MEASURED = read from a result file named here. PROPOSED = a choice of mine. UNVERIFIED = not checked. All calls are PROVISIONAL (the three pathway files have `verified: false`).
- Method as in `results/gainloss/carotenoid_noMicro_blasto_5steps/` (read only): FastTree support tree of the genus representatives, 30 Microsporidia tips dropped (tree identical to the carotenoid one, `cmp`), rooted on the 6 Blastocladiomycota tips, Mk models ER and ARD (castor), model with the lower AIC, a branch is an event when both ends have probability 0.8 or more, tree-aware simulation null (N = 2,000; as-specified and prevalence-matched, window 10%), Pagel (phytools). No BUSCO filter.
- **Only placed steps are modelled.** `hppd` is "family present" and is not modelled. `atp_grasp_ligase` is a search-hit step. `mg_glucosyltransferase` and `glutamicol_formation` have no placed call.
- **Several placed calls are close to search-hit calls** because almost every candidate falls in the target region (section 2). A "loss" is an absent call: a missing hit can be a real loss, a cutoff or placement effect, or a missing gene in one genome.
- One genome per genus: 1,609 tips with data (1,613 in the pruned tree; 4 representatives not searched).
- **The tree-aware null has a floor.** N = 2,000 gives p of at least 0.0005. With 144 tests per scope in the DHN set, BH q cannot go below 0.036 for two tests at that floor (section 4).
- Three sets, never mixed in one model: `melanin_dhn` (9 steps), `melanin_ldopa` (2 modelled steps) and `mycosporine` (2 steps). A fourth run, `melanin_dhn_anypks`, adds the union column "any DHN PKS" with the six non-PKS DHN steps (section 3).
- Nothing here has been tested at the bench.
:::

Paths are relative to `analysis/pathway_census/` unless they start with `docs/`. Run records with commands, versions and jobs: `results/gainloss/<set>/RUN.md`. Ready-to-paste report fragments were prepared for the two reports (see section 7, not saved as files).

## 1. Calls used

Source of each step: `assigned.PROVISIONAL.tsv` of the placement run (label `assigned` = present), through `bin/make_calls_table.py --prefix-aliases prefix_aliases.tsv`. Searched genomes: `results/search/melanin_myco_full/searched.chunk*.tsv` (melanin worktree). All 1,654 genus representatives; 4 are NA.

| Set | Step | Placement run used | Present / absent tips with data (tree) |
|---|---|---|---:|
| DHN | thn_synthase | `pks_joint_id0.8`, per-step rule | 348 / 1,261 |
| DHN | hexaketide_pks | `pks_joint_id0.8`, per-step rule | 85 / 1,524 |
| DHN | heptaketide_pks (the Alb1/AlbA/wA clade) | `pks_joint_id0.8`, per-step rule | 9 / 1,600 |
| DHN | ayg1 | `ayg1` | 577 / 1,032 |
| DHN | t4hn_reductase | `t4hn_reductase` | 612 / 997 |
| DHN | t3hn_reductase | `t3hn_reductase` | 141 / 1,468 |
| DHN | scytalone_dehydratase | `scytalone_dehydratase_v2_allplace` | 568 / 1,041 |
| DHN | laccase_abr2 | `laccase_abr2` | 98 / 1,511 |
| DHN | mco_abr1 | `mco_abr1` | 84 / 1,525 |
| L-DOPA | tyrosinase | `tyrosinase_v2_allplace` | 1,146 / 463 |
| L-DOPA | hmga | `hmga` | 11 / 1,598 |
| mycosporine | ddgs | `ddgs_allplace` | 664 / 945 |
| mycosporine | omt | `omt_v2` (no `_allplace` run exists) | 635 / 974 |

*Counts from `summary_steps.tsv` in each set folder (1,609 tips). The 1,654-representative counts, before pruning, are in `calls_all_steps.summary.tsv` (for example heptaketide_pks 10 present).*

Choices where the brief left room (PROPOSED, please check):

- `aspulvinone_e_synthase`: the base run, not `aspulvinone_e_synthase_id50_allplace`. The id50 run is a thinned diagnostic (119 of 168 proteins unassigned); the base run has 168 of 168 called. Either way the step is not modelled (below).
- `laccase_diphenol`: `laccase_diphenol_v2_allplace` (402 of 402 target; genomes 215).
- PKS: `pks_joint_id0.8`, not the `_allplace` folder. The all-placements note records 0 changed calls (assigned counts equal: 3,163; 609; 2,123 proteins).
- **Reading trap, the "present" counts differ between tables:** 1,654 representatives (before pruning) against 1,609 tips with data; genomes (for example tyrosinase 14,219 genomes) against genus tips (1,146).

### Steps not modelled (too few present tips)

PROPOSED rule: a step needs at least 5 present tips to fit a model. Among the 1,650 representatives with data (`calls_all_steps.summary.tsv`): `laccase_diphenol` 1 present, `asp_melanin_tyrosinase` 0, `aspulvinone_e_synthase` 0. They are listed in `calls_all_steps.tsv` but not analysed. 215, 42 and 49 genomes carry them (`n_call_genomes_not_in_reps` 214, 42, 49 outside the representatives). With 0 present tips a Mk model cannot be fitted. With 1 present tip any rate would be unidentifiable (UNVERIFIED, not tried). This is a negative result for the L-DOPA set: only tyrosinase and hmga are modelled.

### Calls that are close to search hits

Placed target against candidates (`other_steps_summary.tsv`): ayg1 11,835 of 11,835; mco_abr1 1,918 of 1,918; tyrosinase_v2_allplace 50,447 of 50,447; hmga 2,313 of 2,313; ddgs_allplace 8,864 of 8,865; t3hn_reductase 2,236 of 2,240. For these steps the placement removes almost nothing and the call is the HMM search hit. The figures mark them "target = all candidates".

## 2. Models and events (MEASURED, `summary_steps.tsv`)

ARD has the lower AIC for all 13 modelled steps. ER and ARD give the same loss counts for some steps and very different gain counts for others, so a gain count is only meaningful with its model.

| Set | Step | Delta AIC (ARD - ER) | Gains (tip / internal) | Losses (tip / internal) | Uncertain | Parsimony gains | Gains under ER |
|---|---|---:|---:|---:|---:|---:|---:|
| DHN | thn_synthase | -97.0 | 6 (6 / 0) | 36 (24 / 12) | 25 | 11 to 15 | 10 |
| DHN | hexaketide_pks | -63.7 | 0 (0 / 0) | 43 (34 / 9) | 143 | 24 to 26 | 24 |
| DHN | heptaketide_pks | -6.2 | 2 (0 / 2) | 0 | 3 | 2 to 3 | 2 |
| DHN | ayg1 | -71.6 | 14 (11 / 3) | 66 (53 / 13) | 27 | 36 to 37 | 28 |
| DHN | t4hn_reductase | -91.0 | 4 (3 / 1) | 76 (57 / 19) | 27 | 23 to 30 | 23 |
| DHN | t3hn_reductase | -91.6 | 0 | 18 (17 / 1) | 19 | 4 to 6 | 4 |
| DHN | scytalone_dehydratase | -68.8 | 16 (10 / 6) | 44 (34 / 10) | 30 | 20 to 26 | 17 |
| DHN | laccase_abr2 | -75.1 | 7 (7 / 0) | 39 (34 / 5) | 458 | 63 to 66 | 61 |
| DHN | mco_abr1 | -193.0 | 1 (1 / 0) | 44 (40 / 4) | 116 | 33 to 38 | 29 |
| L-DOPA | tyrosinase | -6.6 | 14 (12 / 2) | 47 (41 / 6) | 68 | 26 to 28 | 19 |
| L-DOPA | hmga | -4.8 | 3 (1 / 2) | 1 (1 / 0) | 0 | 3 | 3 |
| mycosporine | ddgs | -151.1 | 8 (7 / 1) | 139 (123 / 16) | 111 | 67 to 80 | 35 |
| mycosporine | omt | -247.9 | 3 (3 / 0) | 152 (137 / 15) | 93 | 59 to 73 | 21 |

Reading:

- **Weak steps.** heptaketide_pks has 9 present tips, hmga 11. Their ARD fits give a loss rate of 1.26 and 0.81 against gain rates of 0.006 and 0.009. With so few present tips the rates are poorly determined (UNVERIFIED, no confidence interval computed). heptaketide_pks: ARD is better than ER by only 6.2 AIC units.
- **Unstable ARD fits.** laccase_abr2 (458 uncertain branches, loss rate 6.27) and hexaketide_pks (143 uncertain, loss rate 5.02) have a high loss rate and many branches that stay below 0.8. For hexaketide_pks the 0 gains under ARD and 24 under ER show the dependence on the model.
- **Most events are on tip edges** (for example ddgs 123 of 139 losses; t4hn_reductase 57 of 76), so a missed gene in one genome would look the same.
- Gains of heptaketide_pks are on 2 internal edges; the larger has 9 tips (77.8% Eurotiomycetes) (`edge_events_primary.tsv`). The 9 present tips are the whole group: the call fits one gain of a Eurotiomycetes clade, not a repeated pattern (UNVERIFIED as biology).
- **Distribution by taxon (MEASURED from `calls.tsv`):** DHN steps are present almost only in Ascomycota (for example thn_synthase 348 of 1,038 Ascomycota tips and 0 of 429 Basidiomycota). Tyrosinase is present in Ascomycota (870 of 1,038), Basidiomycota (257 of 429) and a few Chytridiomycota (8 of 36) and Mucoromycota (1 of 40). ddgs and omt are in Ascomycota (579 and 614) and in part of Basidiomycota (88 and 24); both are absent from all 67 Pichiomycetes, 43 Saccharomycetes, 40 Mucoromycota and 36 Chytridiomycota tips.

## 3. "Any DHN PKS" (`melanin_dhn_anypks`)

`any_dhn_pks` = 1 if any of thn_synthase, hexaketide_pks, heptaketide_pks is present (`bin/add_union_step.py`). It is run with the six non-PKS DHN steps, not with its three members, because a union and its members always share events and would give a trivial pair. Present 402, absent 1,207 tips (ARD lower AIC 457.6 against ER 543.6).

**`any_dhn_pks` is present in 402 of 1,609 tips (25.0%).** ARD: 7 gains (5 tip, 2 internal), 41 losses (30 tip, 11 internal), 24 uncertain branches; ER gives 14 gains and 28 losses. Parsimony gains 16 to 19. Matched null (N = 2,000): present tips 399.4 (402 observed), gains 6.4, losses 68.0 (the null has more losses than observed). Co-events with the non-PKS steps (matched null, all edges; BH q over the 84 tests of this 7-step run, so not comparable with the 144 tests of the main run):

| Pair | Event | Observed (internal) | Null mean | p | q |
|---|---|---:|---:|---:|---:|
| any_dhn_pks x t4hn_reductase | co-loss | 21 (6) | 3.85 | 0.0005 | 0.014 |
| any_dhn_pks x ayg1 | co-loss | 17 (6) | 3.24 | 0.0015 | 0.025 |
| any_dhn_pks x scytalone_dehydratase | co-loss | 15 (4) | 3.27 | 0.014 | 0.098 |
| any_dhn_pks x laccase_abr2 | co-loss | 6 (0) | 0.60 | 0.0045 | 0.054 |
| any_dhn_pks x scytalone_dehydratase | co-gain | 2 (2) | 0.11 | 0.008 | 0.067 |

Co-loss with t4hn_reductase and ayg1 has q below 0.05. Co-loss with t3hn_reductase (4 against 1.46) and mco_abr1 (5 against 1.26) does not. Pagel (6 pairs, 792 to 2,294 s): all p below 0.001 (likelihood ratio 25.3 to 174.6); not read as shared events. Output: `results/gainloss/melanin_dhn_anypks/` and `melanin_dhn_anypks_null/`. The as-specified (stationary) null for this run was stopped after task 1 of 4 (500 replicates) to save compute; its partial files are in `melanin_dhn_anypks_null/stationary_partial_task1_only/` and are not used.

The per-step values of the six non-PKS steps are the same as in the main DHN run (gains, losses, uncertain are equal; the rates differ in the fourth digit because the optimiser restarts).

## 4. Shared events against the tree-aware null (MEASURED, `<set>_null/null_summary.tsv`, variant `stationary_match0.1`, all edges)

N = 2,000 replicates per variant (4 array tasks of 500). The as-specified null (stationary root, no matching) is misspecified for these steps: it simulates a rare trait. For DHN the mean present tips are 3.9 (hexaketide_pks, 85 observed), 28.8 (t4hn_reductase, 612), 52.6 (thn_synthase, 348) and 7.6 (mco_abr1, 84); for mycosporine 103 and 48 against 664 and 635. The matched null reproduces present tips (for example t4hn_reductase 610.0 against 612) and is the null read here. In the matched null the losses are above the observed counts for most steps (thn_synthase 60.9 against 36; ayg1 81.8 against 66; t4hn_reductase 93.4 against 76), which makes the co-loss test conservative; laccase_abr2 is the other way (gains 23.7 in the null against 7 observed; losses 25.4 against 39).

### DHN (36 pairs; 72 co-gain and co-loss tests; BH q over the 144 tests of the scope)

**Five tests have q below 0.05: co-loss of t4hn_reductase with thn_synthase (17 branches observed, null mean 3.10, p 0.0015, q 0.043), ayg1 (20 against 6.01, p 0.0010, q 0.043), laccase_abr2 (9 against 1.26, p 0.0015, q 0.043) and mco_abr1 (17 against 2.72, p 0.0005, q 0.036), and co-gain of t4hn_reductase with scytalone_dehydratase (4 against 0.11, p 0.0005, q 0.036).**

| Pair | Event | Observed (internal) | Null mean (95%) | p | q |
|---|---|---:|---:|---:|---:|
| thn_synthase x t4hn_reductase | co-loss | 17 (5) | 3.10 (0 to 11) | 0.0015 | 0.043 |
| ayg1 x t4hn_reductase | co-loss | 20 (4) | 6.01 (1 to 13) | 0.0010 | 0.043 |
| t4hn_reductase x laccase_abr2 | co-loss | 9 (1) | 1.26 (0 to 5) | 0.0015 | 0.043 |
| t4hn_reductase x mco_abr1 | co-loss | 17 (1) | 2.72 (0 to 8) | 0.0005 | 0.036 |
| t4hn_reductase x scytalone_dehydratase | co-gain | 4 (1) | 0.11 (0 to 1) | 0.0005 | 0.036 |

*Observed branch counts and null means from `summary_pairs.tsv`; internal-edge observed counts from the same table (null means for internal edges in `null_summary.tsv`: thn x t4hn 0.456, p 0.001; ayg1 x t4hn 0.89, p 0.019; t4hn x laccase_abr2 0.13, p 0.12; t4hn x mco_abr1 0.36, p 0.29).*

- **The tree-blind permutation test gives more pairs.** 16 of 72 co-gain and co-loss tests have q below 0.05 (`perm_q`), the null 5. The tree-aware null is the test used.
- **No pair involving heptaketide_pks has an observed co-event** (0 in all 8 pairs); it has only 9 present tips and 0 losses.
- **hexaketide_pks pairs are within chance after correction** (for example with ayg1: 8 against 2.99, p 0.050, q 0.45; with t4hn_reductase: 9 against 3.44, p 0.033, q 0.36). thn_synthase x hexaketide_pks co-loss is 5 against 1.10 (p 0.078, q 0.63).
- **q has a floor set by N.** Only the tests at p = 0.0005 reach q 0.036. A rank-2 test at the floor has q = 0.0005 x 144 / 2 = 0.036. Several tests with p between 0.003 and 0.017 (thn x ayg1 15 against 2.79, p 0.003, q 0.072; thn x scytalone 13 against 2.72, p 0.017) would need more replicates to settle (not done).
- Reading trap: **t4hn_reductase and t3hn_reductase are two reductases of one pathway but are not the same gene.** Their co-loss is 2 against 1.37 (p 0.33). I did not test why (UNVERIFIED).
- The null does not model missing data, a shared technical cause or a link between steps. Co-loss of steps with calls near a search hit (ayg1, t4hn_reductase) can have a shared cause such as a divergent family member in one lineage.

### L-DOPA

tyrosinase x hmga: no event shared in any of the 4 event types (observed 0; null means 0.027 co-gain, 0.038 co-loss). Not informative: hmga has 3 gains and 1 loss. **No test is possible beyond "no overlap".**

### Mycosporine

**ddgs x omt co-loss: 79 branches observed against a null mean of 16.68 (95% 7 to 29), p = 0.0005, q = 0.002.** Internal edges: 6 against 2.19 (p 0.031, q 0.124). Co-gain: 1 against 0.06 (p 0.059, q 0.118). 73 of the 79 are tip edges (Sordariomycetes 20, Lecanoromycetes 15, Tremellomycetes 11, Dothideomycetes 10, Leotiomycetes 10, Pezizomycetes 3, Eurotiomycetes 1, unlabelled 2); 6 are internal edges of 2 to 6 tips. The q is over 4 tests (one pair), not 40. omt took 1,472 draws per accepted replicate in the matched null (ddgs 189).

## 5. Pagel (MEASURED, `<set>/summary_pagel.tsv`)

| Set | Pairs run | Result |
|---|---:|---|
| DHN | 36 of 36 | 29 with p below 0.05 (26 below 0.001); 7 with p of 0.18 or more |
| L-DOPA | 1 | tyrosinase x hmga: LR 14.64, 4 df, p 0.0055; 1,537 s |
| mycosporine | 1 | ddgs x omt: LR 311.96, 4 df, p 2.8e-66; 803 s |

DHN pairs with p of 0.18 or more: heptaketide_pks with ayg1 (0.49), scytalone_dehydratase (0.18), t3hn_reductase (0.99); thn_synthase x heptaketide_pks (0.34); and three pairs with a negative likelihood ratio (hexaketide_pks x heptaketide_pks -0.92; hexaketide_pks x t3hn_reductase -2.06; t3hn_reductase x mco_abr1 -1.17). **A negative likelihood ratio means that the dependent model fitted worse than the independent model, which cannot happen for a converged fit; these three fits did not converge (UNVERIFIED reason).** Their p is reported as 1.0 by the script; I do not read it as a result. Run time 723 to 2,927 s per pair.

**Pagel is not read as evidence of shared events.** Presence is clustered by clade, so a small p follows from any two clustered traits (as in the carotenoid report). Large Pagel p values for pairs with a rare step (heptaketide_pks, 9 present tips) mean the test has little to fit.

## 6. Negative and unresolved results

- Three L-DOPA steps are not modelled (1, 0 and 0 present genus tips).
- heptaketide_pks (9 present tips) gives 2 gains and 0 losses and no co-event; no statement about its evolution is supported beyond the one clade.
- hexaketide_pks has 0 ARD gains and 143 uncertain branches; the ER model gives 24 gains. The model choice decides the answer.
- Co-gain is not shown for any pair except t4hn_reductase x scytalone_dehydratase (4 branches; internal 1).
- Co-loss within the PKS steps and between thn_synthase and ayg1 or scytalone_dehydratase is above the null mean but not at q below 0.05.
- The DHN co-loss signal centres on t4hn_reductase. I did not examine why.
- The q floor: 144 tests with N = 2,000 (section 4).
- No root comparison, ASTRAL comparison or BUSCO test was run. `gainloss_counts.tsv` has a second root (Chytridiomycota) for each step; I did not read it.

## 7. What was not done, and what would settle it

- **Report fragments as files.** The two fragment files requested (`results/gainloss/<set>/REPORT_FRAGMENT.md`) were refused by the write tool of this session ("subagents should return findings as text"). The fragment text was returned in the session result instead.
- A larger N (for example 20,000 replicates) for the DHN null so that q is not limited by 144 tests at p of 0.0005.
- A test of whether co-loss branches sit on low-completeness genomes: list the genomes at the tip co-loss edges (t4hn_reductase pairs; ddgs x omt 73 tips) and compare BUSCO completeness with the tips that keep both steps. A shift would mean a data effect.
- Placement quality: ayg1, mco_abr1, tyrosinase, hmga and ddgs are search-hit calls in effect. A tree with decoys that remove family members would change what a loss means. PKS: sampling label "not yet adequate", NNI-only bootstrap (owner decisions of 2026-11-25), wA unresolved against Alb1/AlbA, no outgroup. The effect on the gain/loss counts was not measured.
- 5 re-annotated genomes: the prefix alias file was used, but the placement did not re-run for them (BFD owner's decision pending, `2026-11-24-genome-prefix-mismatch.md`).
- Convergence of the three Pagel fits with a negative likelihood ratio.

## 8. Files

| What | Where |
|---|---|
| Calls, counts, events, co-events, Pagel, summaries, figures | `results/gainloss/{melanin_dhn,melanin_dhn_anypks,melanin_ldopa,mycosporine}/` |
| Simulation null | `results/gainloss/{melanin_dhn,melanin_dhn_anypks,melanin_ldopa,mycosporine}_null/` |
| Shared tree | `results/gainloss/tree.noMicrosporidia.nwk` (identical to the carotenoid file) |
| New code | `bin/add_union_step.py`, `bin/summarise_gainloss_set.py`, `reports/melanin_myco_gainloss/figures/bin/fig_gainloss_{tree,counts,coevents}.R` |
| Tests | `tests/test_add_union_step.py`, `tests/test_summarise_gainloss_set.py` |
| Existing code reused unchanged | `bin/make_calls_table.py`, `bin/drop_phylum_calls.py`, `bin/prune_tree_tips.R`, `bin/run_gainloss_roots.R`, `bin/run_pagel.R`, `bin/run_null_sim_models.R`, `bin/summarise_null_models.R`, `slurm/gainloss_*.sbatch` |
