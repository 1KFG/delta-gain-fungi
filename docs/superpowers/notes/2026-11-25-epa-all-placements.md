---
title: "epa-ng placement cap: survey and all-placements rerun"
subtitle: "Melanin and mycosporine steps and the PKS joint trees, default (7 placements) against all placements"
date: 2026-11-25
---

::: {.callout-important}
## Status and limits

- Labels: MEASURED = read from a result file or `sacct` named here. UNVERIFIED = not checked. All calls are PROVISIONAL (pathways `verified: false`).
- Only the placement stage was repeated. No tree, alignment or clade file was rebuilt. Any change below comes from the epa-ng filter alone.
- epa-ng default: `--filter-acc-lwr 0.99 --filter-max 7`. The rerun uses `--filter-acc-lwr 0.999 --filter-max 10000` (PROPOSED). `--filter-acc-lwr 1.0` crashes `upo_place assign`; `--filter-min-lwr` cannot be combined with `--filter-acc-lwr` (`other_steps_RUN.md`).
- "Distinct leaf sequences" = distinct sequences in the trimmed reference alignment `ref_set.trim.afa` (all leaves: references, decoys, named members, outgroups, representatives). It is not the number of representatives. Example: `laccase_diphenol_v2` has 87 leaves, 85 distinct, 67 representatives.
- "Cap hit" = a query with `n_placements >= 7` in `assign.tsv`. In the old runs this means the list was cut at 7. In the `_allplace` runs there is no cut, so the same column only counts queries with 7 or more placements.
- Not rerun (reasons in section 4): the other `_id50` diagnostic runs, `t4hn_reductase`, `t3hn_reductase`, `tyrosinase_par`, all carotenoid steps. 7 single-step runs and 2 PKS joint trees (6 step calls) were rerun.
- Fraction denominators are the queries of that step (proteins with a candidate row), not genomes.
:::

## 1. Survey of placed runs under the default (MEASURED)

Source: `results/step5/epa_cap_survey.tsv`, group `default_cap7` (script `bin/step5_epa_cap_report.py survey`, reads `placements_detail*.tsv`, `assign*.tsv`, `ref_set.trim.afa`). Calls from `placements_detail`: target, other side (decoy or outgroup clade), unassigned (best clade LWR below 0.5). `cap-hit unassigned` = unassigned queries that also hit the cap.

| Run | Leaves / distinct | Queries | Target | Other | Unassigned (fraction) | Cap-hit queries | Cap-hit and unassigned |
|---|---:|---:|---:|---:|---:|---:|---:|
| laccase_diphenol | 25 / 25 | 402 | 402 | 0 | 0 (0.0000) | 0 | 0 |
| laccase_diphenol_id50 | 21 / 21 | 402 | 402 | 0 | 0 | 0 | 0 |
| **laccase_diphenol_v2** | 87 / 85 | 402 | 198 | 0 | **204 (0.5075)** | 226 | **204** |
| asp_melanin_tyrosinase | 112 / 112 | 2,471 | 42 | 2,429 | 0 | 0 | 0 |
| asp_melanin_tyrosinase_id50 | 14 / 14 | 2,471 | 0 | 2,306 | 165 (0.0668) | 0 | 0 |
| aspulvinone_e_synthase | 55 / 54 | 168 | 83 | 85 | 0 | 0 | 0 |
| aspulvinone_e_synthase_id50 | 19 / 18 | 168 | 44 | 5 | 119 (0.7083) | 0 | 0 |
| hmga | 76 / 76 | 2,313 | 2,313 | 0 | 0 | 1 | 0 |
| hppd | 305 / 305 | 20,098 | 5,850 | 14,216 | 32 (0.0016) | 191 | 5 |
| hppd_v2 | 667 / 667 | 20,098 | 20,083 | 1 | 14 (0.0007) | 302 | 0 |
| tyrosinase (= tyrosinase_par) | 310 / 310 | 50,447 | 50,426 | 0 | 21 (0.0004) | 88 | 9 |
| tyrosinase_v2 | 1,010 / 1,010 | 50,447 | 50,433 | 0 | 14 (0.0003) | 62 | 4 |
| ayg1 | 283 / 283 | 11,835 | 11,835 | 0 | 0 | 1 | 0 |
| ayg1_id50 | 49 / 49 | 11,835 | 11,826 | 0 | 9 (0.0008) | 1 | 0 |
| laccase_abr2 | 267 / 267 | 9,631 | 3,809 | 5,822 | 0 | 0 | 0 |
| laccase_abr2_id50 | 41 / 41 | 9,631 | 1,045 | 7,778 | 808 (0.0839) | 0 | 0 |
| mco_abr1 | 124 / 124 | 1,918 | 1,918 | 0 | 0 | 0 | 0 |
| mco_abr1_id50 | 22 / 22 | 1,918 | 1,648 | 0 | 270 (0.1408) | 0 | 0 |
| scytalone_dehydratase_v2 | 218 / 218 | 8,604 | 8,050 | 533 | 21 (0.0024) | 255 | 4 |
| t4hn_reductase | 281 / 280 | 10,010 | 8,457 | 666 | 887 (0.0886) | 8 | 0 |
| t3hn_reductase | 281 / 280 | 2,240 | 2,236 | 0 | 4 (0.0018) | 0 | 0 |
| ddgs | 156 / 156 | 8,865 | 8,864 | 1 | 0 | 90 | 0 |
| ddgs_id50 | 33 / 33 | 8,865 | 8,861 | 3 | 1 | 0 | 0 |
| omt_v2 | 311 / 310 | 14,439 | 11,730 | 2,703 | 6 (0.0004) | 104 | 0 |
| PKS 0.7, per-step rule: thn_synthase | 264 / 264 | 7,286 | 3,163 | 3,677 | 446 (0.0612) | 0 | 0 |
| PKS 0.7, per-step: hexaketide_pks | 264 / 264 | 8,604 | 609 | 7,505 | 490 (0.0570) | 0 | 0 |
| PKS 0.7, per-step: heptaketide_pks | 264 / 264 | 2,200 | 2,123 | 0 | 77 (0.0350) | 0 | 0 |
| PKS 0.8, per-step: thn_synthase | 559 / 558 | 7,286 | 3,163 | 4,005 | 118 (0.0162) | 0 | 0 |
| PKS 0.8, per-step: hexaketide_pks | 559 / 558 | 8,604 | 609 | 7,838 | 157 (0.0182) | 0 | 0 |
| PKS 0.8, per-step: heptaketide_pks | 559 / 558 | 2,200 | 2,123 | 7 | 70 (0.0318) | 0 | 0 |

The PKS 0.7 run also has a union-decoy table (`pks_joint/`, same counts for thn_synthase and hexaketide_pks, no clade for heptaketide_pks); it is in `epa_cap_survey.tsv`. `tyrosinase_par` has the same counts as `tyrosinase`.

### Flags

- Flag A, cap is the direct cause of unassigned calls: `laccase_diphenol_v2` (all 204 unassigned queries hit the cap; 204 of 402 = 50.7%). Flag A, smaller: `hppd` (5 of 32), `tyrosinase` (9 of 21), `tyrosinase_v2` (4 of 14), `scytalone_dehydratase_v2` (4 of 21). A cap hit is necessary for the artefact, not sufficient (section 3).
- Flag B, few distinct sequences and high unassigned fraction, but no cap hit: `aspulvinone_e_synthase_id50` (18 distinct, 0.7083), `asp_melanin_tyrosinase_id50` (14, 0.0668), `laccase_abr2_id50` (41, 0.0839), `mco_abr1_id50` (22, 0.1408). Cap hits are 0, so the 7-placement cut is not the reason (a query with fewer than 7 placements was stopped by the 0.99 cumulative-weight limit, not by the cap; epa-ng filter semantics from its documentation, UNVERIFIED here). Not an epa-ng artefact by the data here.
- Flag C, unassigned fraction above 0.05 with many distinct sequences and few cap hits: `t4hn_reductase` (280 distinct, 0.0886, 8 cap hits, 0 of them unassigned). PKS 0.7 (0.0350 to 0.0612) and PKS 0.8 (0.0162 to 0.0318): 0 cap hits at any query.
- Few distinct sequences with no problem: `laccase_diphenol` v1 (25), `laccase_diphenol_id50` (21): 0 unassigned. No call changes expected; not rerun.

## 2. Rerun (MEASURED, `sacct`, partition `intel`, at most 14 CPUs of mine at once)

| Run folder (new, suffix `_allplace`) | Job | CPUs | Elapsed |
|---|---|---:|---|
| laccase_diphenol_v2_allplace | 32461561 | 2 | 00:01:03 |
| scytalone_dehydratase_v2_allplace | 32461562 | 2 | 00:01:18 |
| ddgs_allplace (control) | 32461563 | 2 | 00:02:55 |
| aspulvinone_e_synthase_id50_allplace | 32461564 | 2 | 00:00:12 |
| tyrosinase_v2_allplace | 32461558 | 4 | 01:00:17 |
| tyrosinase_allplace | 32461559 | 4 | 00:28:34 |
| hppd_allplace | 32461560 | 4 | 00:05:12 |
| pks_joint_allplace (0.7, per-step) | 32461566 | 4 | 00:43:13 |
| pks_joint_id0.8_allplace (0.8, per-step) | 32461567 | 4 | 00:48:15 |

All COMPLETED. Single-step runs: `submit_allplace.sh` (copies inputs, `step5_place_parallel.sbatch`). PKS: `submit_allplace_joint.sh` (copies tree and alignment, reruns `step5_joint.py clades --decoy-scope per-step --step-decoys`, then `step5_place_joint.sbatch`). `step5_place_joint.sbatch` now reads `PC_EPA_ARGS` (unset = unchanged). The clade files of both PKS reruns are byte-identical to the old per-step files (`cmp`). The dependency of jobs 32461559 and 32461560 on 32461558 was removed with `scontrol` to use free CPUs; this does not change results.

## 3. Old against new (MEASURED; `results/step5/epa_allplace_transitions.tsv`, `other_steps_summary.allplace.tsv`)

Proteins compared are the same set in each pair (only_old = only_new = 0 in all pairs). "Target/other flips" = proteins that went from target to other side or the reverse.

| Step | Target old to new | Other old to new | Unassigned old to new | Target/other flips | Genomes assigned old to new |
|---|---|---|---|---:|---|
| laccase_diphenol_v2 (402) | 198 to 402 | 0 to 0 | 204 to 0 | 0 | 135 to 215 |
| tyrosinase_v2 (50,447) | 50,433 to 50,447 | 0 to 0 | 14 to 0 | 0 | 14,216 to 14,219 |
| tyrosinase (50,447) | 50,426 to 50,446 | 0 to 0 | 21 to 1 | 0 | 14,212 to 14,219 |
| hppd (20,098) | 5,850 to 5,850 | 14,216 to 14,234 | 32 to 14 | 0 | 5,685 to 5,685 |
| scytalone_dehydratase_v2 (8,604) | 8,050 to 8,050 | 533 to 533 | 21 to 21 | 0 | 7,408 to 7,408 |
| ddgs, control (8,865) | 8,864 to 8,864 | 1 to 1 | 0 to 0 | 0 | 8,760 to 8,760 |
| aspulvinone_e_synthase_id50 (168) | 44 to 44 | 5 to 5 | 119 to 119 | 0 | 41 to 41 |
| PKS 0.7: thn_synthase, hexaketide_pks, heptaketide_pks | unchanged | unchanged | unchanged (446, 490, 77) | 0, 0, 0 | 3,124; 600; 2,091 unchanged |
| PKS 0.8: thn_synthase, hexaketide_pks, heptaketide_pks | unchanged | unchanged | unchanged (118, 157, 70) | 0, 0, 0 | 3,124; 600; 2,091 unchanged |

Changed calls are all `unassigned` to a call: laccase_diphenol_v2 204 (to target), tyrosinase_v2 14 (to target), tyrosinase 20 (to target), hppd 18 (to other side, decoy clade). No protein changed between target and other side in any run (0 of 7 single-step runs, 0 of 6 PKS step calls). Assigned proteins (target with support at least 0.5) equal the target count in each new row.

- Control: `ddgs` has 156 distinct leaf sequences, 0 unassigned under the default and 90 cap-hit queries (1.02% of 8,865). With all placements: 0 changed calls of 8,865 (MEASURED, `compare`). Cap-hit queries (n_placements 7 or more) rise to 161 in the new run because the list is no longer cut. The option does not change calls where the default already assigned every query.
- `laccase_diphenol_v2_allplace` reproduces the earlier `laccase_diphenol_v2_epafilter` row (402 / 0 / 0, genomes 215), the first measurement of this effect.
- Support minima of the new rows: laccase_diphenol_v2 0.999; tyrosinase_v2 0.808; tyrosinase 0.952; hppd 0.531.

## 4. Negative results and what was not shown

- PKS joint trees (0.7 and 0.8): 0 cap-hit queries even in the all-placements run (largest `n_placements` seen in `assign.thn_synthase.tsv` of the 0.8 run: 6). Their unassigned fractions (1.6% to 6.1%) are not an epa-ng cutoff artefact. The rerun changed 0 of 7,286 + 8,604 + 2,200 calls per tree.
- `scytalone_dehydratase_v2`: 4 cap-hit unassigned queries in the old run, 0 changed. A cap hit does not guarantee the call changes.
- `aspulvinone_e_synthase_id50`: 119 of 168 unassigned stays 119 (0 cap hits). Low identity thinning, not the filter, is the likely reason (UNVERIFIED; no test of that here).
- `hppd`: 14 unassigned remain (0.0007 of 20,098).
- Not rerun, UNVERIFIED by rerun: `asp_melanin_tyrosinase_id50`, `laccase_abr2_id50`, `mco_abr1_id50` (0 cap hits; no change expected); `t4hn_reductase` (887 unassigned, 8 cap-hit queries, none unassigned; the 8 could still change target or other, not tested); `t3hn_reductase` (0 cap hits). `omt_v2`, `hppd_v2`, `ayg1`, `ddgs_id50` etc. have 0 cap-hit unassigned queries.
- Not shown: whether a 0.999 cumulative weight (rather than all placements) is enough in every tree. Lower weights per placement were not examined. The effect on the gain/loss analysis was not run.
- The 5 re-annotated genomes are not in these tables (known, step 4 note).

## 5. Carotenoid steps (note only, not rerun)

Source: `results/step5/carotenoid/*` of the merged work, read only (`epa_cap_survey.tsv`, group `carotenoid_default_not_rerun`; 18 runs; `lycopene_cyclase` has no `assign.tsv`/`placements_detail.tsv` at that path and was not surveyed).

- Fewest distinct leaf sequences: `tdh` 6 (unassigned 1,389 of 6,142 = 0.2261), `phytoene_desaturase` 15 (24 of 9,784), `cart` 27 (230 of 5,533 = 0.0416), `tsp2` 63 (8 of 1,196), `lycopene_cyclase_v2` 66 (0 unassigned).
- Unassigned fraction above 0.05: `tdh` 0.2261, `cars` 0.1657 (132 leaves), `ggpp_synthase` 0.1103 (374 distinct; 5,787 of 52,474), `acaa` 0.0802, `cars_v2` 0.0676, `tsp3` 0.0648.
- Cap-hit and unassigned: only `ggpp_synthase` (23 of 5,787 unassigned; 348 cap-hit queries). 0 for `tdh`, `cars`, `cart`, `acaa`, `tsp3`. So the unassigned counts of the few-sequence carotenoid steps are not caused by the 7-placement cut. `ggpp_synthase` is the one carotenoid step worth a rerun (UNVERIFIED, not done: asked not to).

## 6. Not done

- No tree, bootstrap or rarefaction rerun. No gain/loss rerun. No change to any old result folder or to `other_steps_summary.tsv` rows (new rows appended, 7 rows `*_allplace`).
- `PC_EPA_ARGS` was added to `step5_place_joint.sbatch` only; the default behaviour is unchanged when it is unset.
- Tests: `tests/test_step5_epa_cap_report.py` (4 tests, new). Full suite: 974 passed, 1 skipped, 1 deselected (`test_step5_place.py::test_extract_end_to_end_with_failures`, a known failure since the carotenoid single-reference record).
- Open question for the owner: the unassigned call in `laccase_diphenol_v2` was an artefact of the default; `other_steps_summary.tsv` row `laccase_diphenol_v2` (198 target) is kept as a record of the artefact, the `_allplace` row (402 target) is the one to use. Which old `*_v2` values were quoted elsewhere was not checked.
