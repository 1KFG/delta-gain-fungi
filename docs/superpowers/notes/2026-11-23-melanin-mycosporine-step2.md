---
title: "Census steps 1 and 2: DHN melanin, L-DOPA melanin, mycosporine"
subtitle: "Step HMMs, reference and decoy scores, proposed min_score, sign-off table"
date: 2026-11-23
---

::: {.callout-important}
## Status and limits

- All results are PROVISIONAL. The three pathway files have `verified: false`. Every HMM and table carries the `.PROVISIONAL` tag.
- Nothing here searched BFD. The step 3 and step 4 searches were not run. The mycosporine pilot diamond search of 2026-11-22 (job 32359412, TIMEOUT) was not used.
- The scores are bit scores of the best domain of each target. They show how the HMMs separate the reference set and the decoy set. They do not show how they separate real BFD proteins. That needs step 5 placement.
- Leave-one-out (LOO) scores are not independent where two references are close relatives. The lowest LOO score is set by the most divergent reference.
- Decoys are few. This note was first written with four steps without a decoy. The owner decisions of 2026-11-23 (section 9) added decoys for five steps and changed three rows. Three steps still have decoys that score nothing on the step HMM (`scytalone_dehydratase`, `hmga`, `atp_grasp_ligase`). A step without a decoy score has no check of separation. Sections 1 to 8 describe the state at `e49c187`; where a row changed, it carries a marker "(changed 2026-11-23, section 9)".
- Nothing is decided here. The owner decisions are in section 7 with a recommendation each. Labels: MEASURED is read from a file named in the text. PROPOSED is a value or a choice. UNVERIFIED is not checked.
:::

## 1. What was run

| Item | Value |
|---|---|
| Branch | `pathway-census-melanin-myco` (from `origin/main` at `99b26a6`) |
| Tools | `bin/build_step_hmms.py --provisional`, `bin/score_refs_decoys.py` (via `slurm/score_refs_decoys.sbatch`), `bin/make_min_scores.py`, and the new `bin/target_scores_at_min.py` |
| Software | hmmer 3.4, mafft 7.505, Python 3.12 (`/usr/bin/python3.12`) |
| Step 1 | run on the login node (3 pathways, 18 HMMs built) |
| Step 2 jobs | 32407157 (melanin_dhn, 48 s), 32407158 (melanin_ldopa, 30 s), 32407159 (mycosporine, 18 s). All COMPLETED on the short queue (MEASURED, `sacct`) |
| Add-on jobs | 32407363 (melanin_dhn, 14 s), 32407364 (melanin_ldopa, 8 s), 32407365 (mycosporine, 8 s). All COMPLETED |
| HMMs | `results/hmms/<pathway>/<step>.PROVISIONAL.hmm` |
| Step 2 tables | `results/step2/<pathway>/` (`step_scores`, `leave_one_out`, `targets_at_or_above_min`, `min_scores`, `target_scores_at_min`, all `.PROVISIONAL.tsv`) |

Why a new script: `score_refs_decoys.py` counts targets at or above the lowest score of the own references on the full HMM. The owner rule uses a lower value (the lowest LOO score, or the highest decoy score). So the `targets_at_or_above_min` table of the tool is empty for all 21 steps (MEASURED). That empty table is not a clean result. At the rule value it is not empty. `bin/target_scores_at_min.py` re-runs each step HMM against the same target set and flags every decoy and other-step reference at or above the rule `min_score`. It has 3 unit tests (`tests/test_target_scores_at_min.py`).

## 2. Step 1: steps, references, files

Steps and reference counts agree between the pathway YAML, `refs/<pathway>/<step>.faa` and `provenance.tsv` for every step with a file (MEASURED, 22 + 15 + 8 = 45 references).

| Pathway | Steps in YAML | Steps with a reference file (HMM built) | References | Decoys (records across step files) |
|---|---|---|---|---|
| melanin_dhn | 9 | 9 | 22 | 31 (40 after section 9) |
| melanin_ldopa | 6 | 6 | 15 | 14 (33 after section 9) |
| mycosporine | 5 | 3 | 8 | 16 |

Reference file problems (MEASURED):

- **Steps without a reference.** `mycosporine`: `mg_glucosyltransferase` and `glutamicol_formation` have `refs: []` and no `.faa`. No HMM is possible. `build_step_hmms.py` stops with a mafft error on the first such step unless `--steps` is given. I built the three steps with `--steps ddgs omt atp_grasp_ligase`.
- **Non-amino-acid characters.** `melanin_dhn/hexaketide_pks.faa` (Q9Y7A7) has 37 `X` characters in 2,177 residues. Other files have none. hmmbuild accepts `X`. Effect on the HMM: not checked.
- **Duplicate ids, empty sequences, count mismatches.** None in the 45 reference files or in the decoy files.
- **Same decoy in several files.** A decoy that applies to several steps is written in each step's decoy file (for example Q03149 in three PKS decoy files, P17489 in two laccase files). The scorer treats them as separate targets. Counts in this note use unique accessions.
- **Same protein as reference and decoy.** Q55P57, P38993 and Q4WZB4 are references in one step and decoys in another (by design).
- **Decoy folders missing (at `e49c187`).** No decoy file for `ayg1`, `scytalone_dehydratase`, `hppd`, `hmga`, `aspulvinone_e_synthase`. Fixed 2026-11-23 (section 9).

Domain cuts: none. The build tool does not cut domains. It aligns the whole `.faa`. A cut would have to be done earlier, with `region:` in the YAML and a trimmed sequence in the `.faa`. The carotenoid pathway has 11 `region:` entries. These three pathways have none (MEASURED, grep). All steps here use full-length references. `fused_with` is empty for every step.

## 3. Sign-off table

Columns:

- **Ref range**: scores of the own references on the full HMM (the references trained the HMM, so these are high by construction).
- **LOO min**: lowest score of a held-out reference on an HMM built without it.
- **min_score**: the value from `make_min_scores.py`. Basis: LOO = lowest leave-one-out score; floor = highest decoy score (single reference).
- **Decoys >= min**: unique decoy accessions at or above `min_score`. For a floor step this includes the decoy that sets the floor, by construction.
- **Other refs >= min**: references of other steps of the same pathway at or above `min_score`.

### melanin_dhn

| Step | Refs | HMM length | Ref range | LOO min | min_score | Basis | Highest decoy | Decoys >= min | Other refs >= min | Flag |
|---|---|---|---|---|---|---|---|---|---|---|
| thn_synthase | 1 | 2187 | 5390.1 | single reference | 1908.6 | floor | 1908.6 (Q03149 wA) | 1 | 2 (Q9Y7A7 hexaketide 2008.0; A2QUI2 heptaketide 1933.9) | Floor is set by wA. Two PKS references of other steps score above it. |
| hexaketide_pks | 1 | 2177 | 5263.3 | single reference | 1839.0 | floor | 1839.0 (Q03149 wA) | 1 | 3 (P79068 2030.7; A2QUI2 1899.6; Q4WZA8 1872.0) | Same as above. |
| heptaketide_pks | 2 | 2159 | 4850.3 to 4853.9 | 3656.7 | 3656.7 | LOO | 3904.0 (Q03149 wA) | 1 | 0 | **Decoy above the lowest reference.** wA scores 3904.0. The LOO scores of the two references are 3656.7 and 3665.7. |
| ayg1 | 3 | 406 | 747.3 to 807.6 | 374.2 | 374.2 | LOO | 55.6 (Q93NG6) (changed 2026-11-23, section 9; was none) | 0 | 0 | LOO of Q6BCB5 (Exophiala) is 374.2; the other two are 615.9 and 632.6. Margin to the best decoy 318.6 bits. |
| t4hn_reductase | 5 | 275 | 462.2 to 523.6 | 356.0 | 356.0 | LOO | 342.9 (Q00791 stcU) | 0 | 1 (A7DWM6 t3hn 357.1) | **Margin of 13.1 bits to the best decoy.** A t3hn reference scores above `min_score`. |
| t3hn_reductase | 3 | 282 | 614.7 to 619.2 | 550.4 | 550.4 | LOO | 414.6 (Q00791 stcU) | 0 | 0 | Margin 135.8 bits. LOO range is narrow (550.4 to 570.7). Rule value is 90% of the lowest full-HMM score (candidate for lowering, section 7). |
| scytalone_dehydratase | 4 | 168 | 321.3 to 336.3 | 240.5 | 240.5 | LOO | none: 5 decoys added 2026-11-23, none scores (section 9) | 0 | 0 | The decoys give no check of separation, because none has a hit at E <= 1000. |
| laccase_abr2 | 1 | 587 | 1455.7 | single reference | 464.3 | floor | 464.3 (P17489 yA) | 1 | 0 | Floor set by one laccase decoy. Next decoy 98.2. Floor is 32% of the reference score. |
| mco_abr1 | 2 | 686 | 1213.3 to 1446.5 | 638.2 | 638.2 | LOO | 385.5 (P38993 FET3) | 0 | 0 | Margin 252.7 bits. |

### melanin_ldopa

| Step | Refs | HMM length | Ref range | LOO min | min_score | Basis | Highest decoy | Decoys >= min | Other refs >= min | Flag |
|---|---|---|---|---|---|---|---|---|---|---|
| tyrosinase | 7 | 563 | 394.7 to 839.8 | 63.6 | 63.6 | LOO | 36.2 (B1VTI5 grixazone synthase) | 0 | 0 | **LOO min is far below the others.** Q00234 (melO) gives 63.6. The next lowest LOO is 239.8 (A0A1S9DK56), then 377.5 (P00440). Margin to the best decoy is 27.4 bits. |
| laccase_diphenol | 3 | 624 | 1337.1 to 1437.6 | 1127.7 | 1127.7 | LOO | 180.5 (Q8LPL3 ascorbate oxidase) | 0 | 0 | Margin 947.2 bits. Rule value is 84% of the lowest full-HMM score (candidate for lowering). |
| hppd | 1 | 403 | 974.9 | single reference | 164.3 (changed 2026-11-23, section 9; was undefined) | floor | 164.3 (O52791 HmaS) | 1 | 0 | Floor set by one decoy, 4-hydroxymandelate synthase, which is the closest HPPD paralogue. The other 6 decoys do not score. Floor is 17% of the reference score. |
| hmga | 2 | 450 | 1069.3 to 1078.5 | 983.2 | 983.2 | LOO | none: 6 decoys added 2026-11-23, none scores (section 9) | 0 | 0 | Rule value is 92% of the lowest full-HMM score. The two references are close relatives (LOO 983.2 and 983.4). The decoys give no check of separation. |
| aspulvinone_e_synthase | 1 | 925 | 2248.1 | single reference | 1205.5 (changed 2026-11-23, section 9; was undefined) | floor | 1205.5 (Q0CT94 atrA) | 1 | 0 | Floor set by the atromentin synthetase of *A. terreus*. Two more NRPS-like decoys score 824.3 and 634.3. Floor is 54% of the reference score. The floor does not separate atromentin synthetases from the step. |
| asp_melanin_tyrosinase | 1 | 356 | 878.6 | single reference | 320.5 | floor | 320.5 (Q5BGU9 ivoB) | 1 | 0 | Floor set by ivoB, a tyrosinase-family decoy. Next decoy 283.8 (Q2UNF9). The tyrosinase references score 7.8 to 21.0 on this HMM. |

At `e49c187` the `make_min_scores.py` run on the full `step_scores` table of melanin_ldopa failed for the two undefined steps (MEASURED, `ValueError: step aspulvinone_e_synthase: no value for rule basis`), and the `min_scores` table was made from `step_scores.floor_defined.PROVISIONAL.tsv` (the table without those two rows). After section 9 both steps have a floor and the full table is used.

### mycosporine

| Step | Refs | HMM length | Ref range | LOO min | min_score | Basis | Highest decoy | Decoys >= min | Other refs >= min | Flag |
|---|---|---|---|---|---|---|---|---|---|---|
| ddgs | 3 | 412 | 779.4 to 827.6 | 463.5 | 463.5 | LOO | 315.5 (Q9ZAE9 acbC) | 0 | 0 | Margin 148.0 bits. The LOO min is the Phaffia reference (WEV88679.1). The two cyanobacterial references are 73% identical (yaml note), so their LOO scores (693.2, 676.6) are not independent. |
| omt | 2 | 335 | 430.0 to 481.0 | 8.3 | 100.6 (owner override 2026-11-23, section 9; rule value 8.3) | override | 100.6 (P21964 human COMT) | 1 (the decoy that sets the override) | 0 | **The rule value (8.3) is below all four decoys.** Each of the two references scores 8.3 or 9.7 when it is held out. The two references are 17% identical (yaml note), so the HMM of one does not detect the other. |
| atp_grasp_ligase | 3 | 500 | 688.7 to 1038.1 | 105.0 | 105.0 | LOO | none (6 decoys, none has a hit at E <= 10) | 0 | 0 | The LOO min is the cyanobacterial reference B2J6X7. The two Phaffia references are 98% identical (490 of 500 residues, `UNVERIFIED.md`), so their LOO scores (955.8, 961.2) are not independent. |

## 4. Steps without a decoy (floor undefined or untested)

This section is the state at `e49c187`. The decoys were added 2026-11-23; see section 9 for what was chosen and what it scores.

| Step | Refs | Consequence | Candidate decoy families (PROPOSED; none fetched, none checked against UniProt) |
|---|---|---|---|
| hppd | 1 | No `min_score`. | Other members of the vicinal oxygen chelate (VOC) superfamily: glyoxalase I, bleomycin resistance protein, extradiol dioxygenases (for example 2,3-dihydroxybiphenyl dioxygenase), methylmalonyl-CoA epimerase. Also HPPD-like proteins of other function (for example hydroxymandelate synthase). |
| aspulvinone_e_synthase | 1 | No `min_score`. | Other NRPS-like adenylation-thiolation-thioesterase enzymes: atromentin synthetase (basidiomycete), other aspulvinone or terphenylquinone synthetases that are not Asp-melanin (not apvA: the YAML says apvA cannot be a decoy), fungal alpha-aminoadipate reductase (Lys2), adenylate-forming acyl-CoA synthetases. |
| ayg1 | 3 | `min_score` defined (LOO). Separation untested. | Other members of the fold of Ayg1 and scytalone dehydratase (NTF2-like superfamily): ketosteroid isomerase, SnoaL-type polyketide cyclases, limonene-1,2-epoxide hydrolase. Which is the fold of Ayg1 is UNVERIFIED. |
| scytalone_dehydratase | 4 | Same. | NTF2-like superfamily members (above), plus the Ayg1 references of the same pathway as a cross-check. |
| hmga | 2 | Same. | Other cupin-fold dioxygenases and the homogentisate dioxygenase family members of other function. UNVERIFIED which family members exist in fungi. |
| atp_grasp_ligase | 3 | Decoys exist, none scores. No decoy score for the floor check. | Already in the file: MysD-like B2J6X6, fsqD, uncJ, aveB, DdlB, GshB. Add other ATP-grasp families: biotin carboxylase, carbamoyl-phosphate synthase large subunit, pyruvate carboxylase, glutathione synthetase-like proteins. |

## 5. Tyrosinase and laccase: do the HMMs score each other's sequences?

Inside `melanin_ldopa` (MEASURED, `target_scores_at_min.PROVISIONAL.tsv`):

- The `tyrosinase` HMM has 10 hit rows. They are 5 unique decoys (19.1 to 36.2; 4 of them appear in two decoy files) and the `asp_melanin_tyrosinase` reference Q0CRX0 (23.1). It has no hit on any laccase reference or laccase decoy.
- The `laccase_diphenol` HMM has a hit on 5 targets. All are laccase-type decoys (75.5 to 180.5). It has no hit on any tyrosinase sequence.

Across the two pathways (ad hoc `hmmsearch` of the three HMMs against all laccase, multicopper oxidase and tyrosinase references; not saved as a table; the command is in the README section):

| HMM | Highest laccase or MCO reference other than its own (score) | Tyrosinase reference with a hit |
|---|---|---|
| melanin_ldopa `laccase_diphenol` | 171.0 (A2R3T0, Abr1 of A. niger); E9RBR0 (Abr2) 108.9 | none |
| melanin_dhn `laccase_abr2` | 102.9 (J9VY90, Cryptococcus LAC) | not scored (not in that target set) |
| melanin_dhn `mco_abr1` | 223.5 (J9VQZ4, Cryptococcus LAC) | not scored |
| melanin_ldopa `tyrosinase` | none | own references only |

Result: the tyrosinase HMM and the laccase HMMs do not score each other's references. Each laccase HMM gives the other laccase-type references 80 to 224 bits. Each `min_score` is above those values.

Every decoy or other-step reference at or above the rule `min_score` of a step (MEASURED, unique accessions):

| Step (min_score) | Target | Score | Kind | Is it a problem for a search? |
|---|---|---|---|---|
| heptaketide_pks (3656.7) | Q03149 wA, *A. nidulans* | 3904.0 | decoy | Yes. wA makes YWA1, the same reaction as the heptaketide step (YAML note). A cutoff cannot separate it. Step 5 placement must. |
| hexaketide_pks (1839.0) | Q03149 wA | 1839.0 | decoy (sets the floor) | By construction. |
| hexaketide_pks | P79068 (thn_synthase), A2QUI2 and Q4WZA8 (heptaketide_pks) | 2030.7, 1899.6, 1872.0 | other-step references | Yes. The three PKS HMMs cannot separate the three PKS steps by score. Step 5 must. |
| thn_synthase (1908.6) | Q03149 wA | 1908.6 | decoy (sets the floor) | By construction. |
| thn_synthase | Q9Y7A7 (hexaketide), A2QUI2 (heptaketide) | 2008.0, 1933.9 | other-step references | Yes. Same reason. |
| laccase_abr2 (464.3) | P17489 yA, *A. nidulans* laccase | 464.3 | decoy (sets the floor) | By construction. The gap to the reference is large (1455.7). Next decoy 98.2. |
| t4hn_reductase (356.0) | A7DWM6 (t3hn_reductase) | 357.1 | other-step reference | Yes, by 1.1 bits. Placement must separate T4HN from T3HN reductases. |
| asp_melanin_tyrosinase (320.5) | Q5BGU9 ivoB | 320.5 | decoy (sets the floor) | By construction. |
| omt (8.3) | P21964 (human COMT) 100.6; O42898 (S. pombe) 96.5; A0A397HK53 (ankF) 57.9; B8NI24 (imqG) 44.8 | 44.8 to 100.6 | decoys | Yes, all four. The cutoff is below every decoy. It does not remove any of them. |

Decoys that score below the rule `min_score` but close to it (within 25% of the value): `t4hn_reductase` Q00791 (342.9 against 356.0). No other step has one.

Note for step 5: all other steps have no decoy or other-step reference at or above `min_score`. That is a statement about the reference and decoy set only. It says nothing about the real BFD proteins.

## 6. Mycosporine: reference set against the pathway steps

- The pathway YAML has five steps. Three have references: `ddgs`, `omt`, `atp_grasp_ligase`. These match the cyanobacterial MysA (DDGS, Npun_R5600), MysB (OMT, Npun_R5599) and MysC (ATP-grasp ligase, Npun_R5598) and the Phaffia rhodozyma DDGS, OMT and ATPG (WEV88679.1, WEV88680.1, WEV88681.1). The YAML uses the Phaffia names, not MysA to MysC (MEASURED, YAML and `provenance.tsv`).
- MysD (a D-Ala-D-Ala ligase-like protein, Npun_F5597, B2J6X6) is not a step. It is a decoy of `atp_grasp_ligase`. The YAML states that no paper that was read shows a second-amine step or an iminomycosporine in a fungus.
- `mg_glucosyltransferase` and `glutamicol_formation` have no reference (tier T4, gene not identified). They cannot be searched with an HMM. They are not in `results/hmms/mycosporine/`.
- Each of the three HMMs is built from one fungal reference (Phaffia) and one or two cyanobacterial references. The YAML calls the cyanobacterial proteins T2 at most and never fungal evidence. Whether a fungal search should use a fungal-only HMM is an owner decision (section 7). Five of the eight references are `verified: false` (`UNVERIFIED.md`).
- Fused or non-fused logic: all three steps have `fused_with: []`. The files that I read do not describe a fungal MysC-MysD fusion or a MysA-MysB fusion. I did not search the literature. Whether such fusions exist in fungi is UNVERIFIED. If one exists, the `atp_grasp_ligase` or `ddgs` HMM would score a fused protein only on the matching domain. The tool reports the best single domain per protein, so a fusion is not lost, but this step does not detect or label fusions. A fused-protein check would belong to step 6 (`assign_steps.py`).
- A related check on the existing decoys (MEASURED): the 1,583-residue pentafunctional AROM decoy (P07547) carries a 3-dehydroquinate synthase domain and scores 99.8 on the `ddgs` HMM. That is below `min_score` 463.5. So a domain-level hit of a fused AROM protein does not pass the `ddgs` cutoff.

## 7. Owner decisions, with a recommendation

Written at `e49c187`: nothing was applied then. No value was lowered. The `--scale` option of `make_min_scores.py` was not used. The owner answers of 2026-11-23 are in section 9.

1. **`omt` rule value (8.3).** The rule value is below all four decoys (44.8 to 100.6). It would pass nearly every O-methyltransferase in a search. The two references are 17% identical, so LOO says only that each HMM cannot detect the other reference.
   - Option A: keep 8.3 and rely on step 5.
   - Option B: use the highest decoy score, 100.6, as a floor.
   - Recommendation: B for the search (PROPOSED). It loses fungal OMTs that score below 100.6. How many there are is UNVERIFIED.
2. **`hppd` and `aspulvinone_e_synthase` (floor undefined).** Recommendation: add decoys (section 4) and re-run step 2 for `melanin_ldopa` before step 3. If you want to go on without decoys, the rule gives no value. I did not choose one.
3. **PKS steps (`thn_synthase`, `hexaketide_pks`, `heptaketide_pks`).** The three HMMs cannot separate the three PKS steps by score. wA (Q03149) is above the heptaketide `min_score`. Recommendation: use the rule values as search floors, and build one joint PKS tree with the references and decoys of all three steps in step 5, as the YAML note already says.
4. **`t4hn_reductase` against `t3hn_reductase`.** The t3hn reference A7DWM6 is 1.1 bits above the t4hn `min_score`. The t4hn decoy margin is 13.1 bits. Recommendation: keep the rule value. Place both steps on one SDR tree in step 5.
5. **Tyrosinase `min_score` (63.6).** One reference (Q00234 melO) sets it. The next LOO value is 239.8. The value is 27.4 bits above the best decoy. Recommendation: keep the rule value. It already separates the 5 decoys. Tyrosinase and laccase get separate trees with decoys inside (owner, 2026-11-22), so the cutoff has no cross-family job.
6. **Candidates for lowering (not applied).** Rule: `min_score` divided by the lowest full-HMM reference score is 0.80 or more, and the highest decoy leaves room. Reason to consider it: the rule value is high because the references are close relatives, so a more divergent ortholog may score below it (UNVERIFIED: no divergent ortholog was tested).

| Step | min_score / lowest full-HMM reference score | Highest decoy | Lowest factor allowed by the decoys (decoy / min_score) |
|---|---|---|---|
| hmga | 0.92 | none | no bound |
| t3hn_reductase | 0.90 | 414.6 | 0.75 |
| laccase_diphenol | 0.84 | 180.5 | 0.16 |

   Recommendation: no lowering now. Decide after step 5 shows whether hits sit just below the rule values. The same ratio for the others is 0.50 to 0.77 (LOO steps) and 0.32 to 0.36 (floor steps), except `tyrosinase` 0.16, `atp_grasp_ligase` 0.15, `omt` 0.02.
7. **Mycosporine HMM composition.** Keep the mixed fungal and cyanobacterial references or build fungal-only HMMs (one Phaffia reference per step, which turns each step into a single-reference floor step). Recommendation: keep as is for now. The pathway decisions are deferred (owner, 2026-11-22), and a fungal-only HMM has no LOO check.
8. **`Q9Y7A7` has 37 `X` characters.** Recommendation: leave it. Replace it with a clean sequence if the owner wants an HMM with no ambiguous residues.
9. **Decoys to add for `ayg1`, `scytalone_dehydratase`, `hmga`.** Recommendation: add them before step 5 (section 4). They do not block step 3 (a `min_score` exists).

## 8. Not done here

- No BFD search (steps 3 and 4). No genome chunking. No tree or placement (step 5).
- No domain cuts, no change to any reference or decoy file, no new decoy fetched.
- No lowered cutoff. No decision on `verified`.
- The cross-pathway scoring in section 5 is an ad hoc check, not a saved table.
- Whether fused MysA-MysB or MysC-MysD proteins exist in fungi is not checked.

## 9. Owner decisions and follow-up (2026-11-23)

Labels: MEASURED is read from a file or a run named here. PROPOSED is a choice or value. UNVERIFIED is not checked.

### 9.1 Owner answers

| Decision | Answer | Applied |
|---|---|---|
| 1 `omt` | Floor about 100 (the highest decoy score, 100.6). A second round may be needed. | Yes, as an explicit override (9.4) |
| 2 `hppd`, `aspulvinone_e_synthase` | Add decoys if they can be identified. | Yes (9.2, 9.3) |
| 3 PKS rule values | Keep. Joint PKS tree later. | No change |
| 4 t4hn / t3hn | Keep. One SDR tree later. | No change |
| 5 tyrosinase 63.6 | Keep. | No change |
| 6 lowering | None now. | No change |
| 7 mycosporine references | Keep the mixed set. | No change |
| 8 Q9Y7A7 | Leave. | No change |
| 9 decoys for `ayg1`, `scytalone_dehydratase`, `hmga` | Add before step 5. | Yes (9.2) |

### 9.2 Decoys added

All 28 accessions were fetched from `rest.uniprot.org` on 2026-11-23 (JSON and FASTA). The script `curation_scripts/melanin_myco_decoys/add_decoys.py` checks for each one: reviewed entry, accession, organism name in the FASTA header, sequence length equal in both fetches, and that the accession is not a reference of any of the four pathways. Organism, gene, name and length below are from the fetched record. A decoy file is `refs/<pathway>/decoys/<step>.faa`. Provenance rows are in `refs/<pathway>/decoys/provenance.tsv` (same columns as the carotenoid decoys). All rows have `verified: true` because all are reviewed entries whose fields were read. `pmid` is empty: no paper was read for these decoys. The decoys are also listed in `pathways/melanin_dhn.yaml` and `pathways/melanin_ldopa.yaml` (the pathway tests require YAML decoys and provenance rows to agree). Three entries have no gene name in UniProt. By the gene-name rule (`2026-11-22-owner-decisions.md`) they carry a protein abbreviation from the record: O52791 `HmaS` and P08659 `Luciferase` (short names in the UniProt record); Q93NG6 has neither, so it carries `Q93NG6_hydrolase_Pnico` (pattern `<accession>_<role>_<species code>`). The reason text says so. In the tables below "none in UniProt" means this case. FASTA headers use the format of the existing decoy files (`>ACC sp|ACC|ENTRY ...`). Score = best domain score on the step HMM (MEASURED, step 2 jobs 32410561 to 32410563; `target_scores_at_min.PROVISIONAL.tsv`). "none" = no hit at E <= 1000 in the add-on scoring.

**Correction of the owner brief for `ayg1`.** The brief names NTF2-like proteins as the Ayg1 family. The UniProt records of the references say otherwise (MEASURED): Q4WZB3 carries Pfam PF06500 (not PF02982, the scytalone dehydratase family), as does the *Exophiala* reference Q6BCB5. Reviewed PF06500 entries are esterases and hydrolases (alpha/beta-hydrolase type; the fold is not stated in the record, UNVERIFIED). So the `ayg1` decoys are PF06500 members of other function. The NTF2-like decoys went to `scytalone_dehydratase` (PF02982 in P56221 and Q00455).

#### hppd (pathway melanin_ldopa; step reference Q4WHU1; floor 164.3)

| Accession | Organism | Gene | Name | Length | Score |
|---|---|---|---|---|---|
| P50107 | *Saccharomyces cerevisiae* S288c | GLO1 | Glyoxalase I | 326 | none |
| Q56415 | *Serratia marcescens* | fosA | Glutathione transferase FosA (fosfomycin resistance) | 141 | none |
| Q99XE5 | *Staphylococcus aureus* Mu50 | ble | Bleomycin resistance protein | 132 | none |
| P06622 | *Pseudomonas putida* | xylE | Metapyrocatechase (catechol 2,3-dioxygenase) | 307 | none |
| A0A0U2WCB2 | *Mycosarcoma maydis* | RDO1 | Glyoxalase domain-containing protein RDO1 | 179 | none |
| Q0D1P8 | *Aspergillus terreus* NIH 2624 | terI | Lactoylglutathione lyase-like protein terB | 197 | none |
| O52791 | *Amycolatopsis orientalis* | none in UniProt | 4-hydroxymandelate synthase (HmaS) | 357 | **164.3** |

The first six decoys were the VOC-superfamily proteins named in the brief. None scores on the HMM. O52791 was added as a seventh, close paralogue with another reaction (EC 1.13.11.46; Pfam PF00903 and PF14696 in the record). Weak points: A0A0U2WCB2 and Q0D1P8 are fungal and listed because the brief asked for a fungal dioxygenase-like decoy, but the records give no function beyond "putative ring-cleaving dioxygenase" and "lactoylglutathione lyase-like" (UNVERIFIED function). Other fungal HPPDs (Q4WPV8, Q872T7, and others) are true HPPDs and were not used.

#### aspulvinone_e_synthase (melanin_ldopa; reference Q0CRX1; floor 1205.5)

| Accession | Organism | Gene | Name | Length | Score |
|---|---|---|---|---|---|
| Q0CT94 | *Aspergillus terreus* NIH 2624 | atrA | Nonribosomal peptide synthetase atrA (atromentin synthetase) | 920 | **1205.5** |
| A7XRY0 | *Aspergillus nidulans* FGSC A4 | tdiA | Didemethylasterriquinone D synthetase | 970 | 824.3 |
| B7STY1 | *Tapinella panuoides* | atrA | Atromentin synthetase atrA | 957 | 634.3 |
| P08659 | *Photinus pyralis* | none in UniProt | Luciferin 4-monooxygenase (luciferase) | 550 | 48.1 |
| Q01574 | *Saccharomyces cerevisiae* S288c | ACS1 | Acetyl-coenzyme A synthetase 1 | 713 | 13.3 |
| P07702 | *Saccharomyces cerevisiae* S288c | LYS2 | L-2-aminoadipate reductase | 1392 | 3.4 |

Weak points: Q0CT94, A7XRY0 and B7STY1 are NRPS-like quinone synthetases of the same architecture as the reference. They make other products (atromentin, terrequinone precursor). Whether the atromentin synthetases should count as "not the step" is a judgment: they are the nearest relatives, and they set a floor of 54% of the reference score. A cutoff at the floor cannot separate them from true aspulvinone E synthases by design; step 5 placement must. apvA (Q0CWD0) is not a decoy (YAML note).

#### ayg1 (melanin_dhn; references Q4WZB3, Q6BCB5, A2R3S8; rule value 374.2 LOO)

| Accession | Organism | Gene | Name | Length | Score |
|---|---|---|---|---|---|
| Q93NG6 | *Paenarthrobacter nicotinovorans* | none in UniProt | 2,6-dihydropseudooxynicotine hydrolase | 367 | 55.6 |
| S0EE80 | *Gibberella fujikuroi* (*Fusarium fujikuroi*) | FUS2 | 20-hydroxy-prefusarin hydrolase FUS2 | 419 | 40.1 |
| P04335 | *Escherichia coli* K-12 | frsA | Esterase FrsA | 414 | 39.3 |
| S8ASK9 | *Penicillium oxalicum* 114-2 | poxO | Hydrolyase poxO | 421 | 37.3 |

Weak points: only four decoys. S0EE80 and S8ASK9 are fungal PF06500 hydrolases of secondary-metabolite clusters. Their reaction differs from the Ayg1 reaction in the records, but they are close in family, so they are the strongest test. Whether they are orthologues of Ayg1 in another role is UNVERIFIED. The verticillium AYG1 (G2X3G9) is a true homologue with the Ayg1 role and is not a decoy.

#### scytalone_dehydratase (melanin_dhn; rule value 240.5 LOO)

| Accession | Organism | Gene | Name | Length | Score |
|---|---|---|---|---|---|
| P07445 | *Pseudomonas putida* | ksi | Steroid Delta-isomerase (ketosteroid isomerase) | 131 | none |
| Q9RN59 | *Streptomyces nogalater* | snoaL | Nogalonic acid methyl ester cyclase (SnoaL) | 134 | none |
| O52646 | *Streptomyces galilaeus* | acma | Aklanonic acid methyl ester cyclase AcmA | 144 | none |
| P61970 | *Homo sapiens* | NUTF2 | Nuclear transport factor 2 | 127 | none |
| Q9ZAG3 | *Rhodococcus erythropolis* | limA | Limonene-1,2-epoxide hydrolase | 149 | none |

Weak points: NTF2-like fold membership of these five is the brief's statement and the Pfam names in the records (PF12680, PF07366, PF02136, PF07858). I did not check a Pfam clan or a structure (UNVERIFIED). None scores, so these decoys test nothing at the HMM level. They are useful as outgroup members in the step 5 tree only.

#### hmga (melanin_ldopa; references Q4WHT9, Q00667; rule value 983.2 LOO)

| Accession | Organism | Gene | Name | Length | Score |
|---|---|---|---|---|---|
| O43123 | *Aspergillus nidulans* FGSC A4 | maiA | Maleylacetoacetate isomerase | 230 | none |
| O43708 | *Homo sapiens* | GSTZ1 | Maleylacetoacetate isomerase (GSTZ1) | 216 | none |
| P47096 | *Saccharomyces cerevisiae* S288c | BNA1 | 3-hydroxyanthranilate 3,4-dioxygenase | 177 | none |
| Q0QFQ2 | *Aquipseudomonas alcaligenes* | hbzE | Gentisate 1,2-dioxygenase 2 | 348 | none |
| Q16878 | *Homo sapiens* | CDO1 | Cysteine dioxygenase type 1 | 200 | none |
| P47125 | *Saccharomyces cerevisiae* S288c | BNA2 | Indoleamine 2,3-dioxygenase | 453 | none |

Weak points: none scores. The brief asked to check whether HmgC is a glutathione S-transferase. MEASURED from the records: O43708 (human, GSTZ1) is named "glutathione S-transferase zeta 1" with EC 2.5.1.18; the fungal O43123 has Pfam PF13410 and PF13409 (GST C-terminal and N-terminal-like domains). So HmgC-type isomerases are GST-family proteins. Whether HmgA (PF04209) shares a fold with the cupin dioxygenases (gentisate dioxygenase, CDO, HAAO) is UNVERIFIED, and nothing here depends on it. P47125 (heme dioxygenase) is probably unrelated to HmgA.

### 9.3 Step 2 re-run

Step 1 (HMMs) was rebuilt in the new worktree for all three pathways. The HMM of every step is unchanged by the decoys. Final step 2 jobs (after the decoy FASTA headers were put in the existing format): 32410561 (melanin_dhn, 1 min 49 s), 32410562 (melanin_ldopa, 1 min 40 s), 32410563 (mycosporine, 1 min 36 s). Add-on `target_scores_at_min` jobs: 32410599, 32410600, 32410601. All COMPLETED (MEASURED, `sacct`). Earlier runs of the same inputs (32410013, 32410014 before O52791 was added, 32410246, 32410015 and add-ons 32410247, 32410248, 32410324) gave the same `min_scores` (MEASURED, `min_scores` files compared by eye: all values equal). Scores and min_score of the other 12 steps agree with the first run (MEASURED: the section 3 tables and the new `min_scores` files).

| Step | min_score before | min_score now | Basis now |
|---|---|---|---|
| hppd | undefined | **164.3** | highest decoy score (single reference, floor), set by O52791 HmaS |
| aspulvinone_e_synthase | undefined | **1205.5** | highest decoy score (single reference, floor), set by Q0CT94 atrA |
| omt | 8.3 | **100.6** | owner override 2026-11-23 (9.4); rule value 8.3 |
| ayg1 | 374.2 | 374.2 | LOO. Highest decoy 55.6 |
| scytalone_dehydratase | 240.5 | 240.5 | LOO. No decoy scores |
| hmga | 983.2 | 983.2 | LOO. No decoy scores |

`hppd` has a floor only because of one added decoy (O52791). With the six VOC decoys of the brief alone, no decoy scored and the floor would be undefined (MEASURED, first run 32410014: `max_decoy_score` empty for `hppd`). A floor of 0 is not used anywhere. `aspulvinone_e_synthase` has a floor from the start of the run. No other-step reference or other decoy scores at or above the new `hppd` or `aspulvinone_e_synthase` min_score (MEASURED, `target_scores_at_min`).

### 9.4 omt override and tool changes

- `bin/make_min_scores.py` has a new option `--floor-override STEP=VALUE` (repeatable) with `--override-basis TEXT`. The rule value stays in `rule_min_score`. The command used for mycosporine: `--floor-override omt=100.6 --override-basis "owner override 2026-11-23: highest decoy score; revisit in second round"`. Result in `results/step2/mycosporine/min_scores.PROVISIONAL.tsv`: omt 100.6, basis as above, rule_min_score 8.3.
- New option `--allow-undefined`: a step with no rule value gets an EMPTY `min_score` and the basis "no floor defined: ...". The default is unchanged (stop with an error). `bin/target_scores_at_min.py` skips such a step and names it on stderr. `bin/search_chunk.py` stops with an error that names the step. This was not needed for the final tables (no step is undefined now). It makes the "no floor defined" case explicit if a later decoy set does not score.
- The committed review packet `docs/superpowers/plans/2026-11-21-stage1-review-packet.md` is generated from the YAML and was regenerated with `bin/make_review_packet.py` (its test `test_committed_packet_is_current` requires this). Its decoy counts changed.
- Tests: 11 new tests in `tests/test_make_min_scores.py`, `tests/test_target_scores_at_min.py` and `tests/test_search_chunk.py` (see the commit message for the full count).

### 9.5 omt: evidence for the second round (not acted on)

Data: `results/omt_check/` (git-ignored). Query: UniProt reviewed entries of fungi (taxonomy 4751), protein name containing "O-methyltransferase"; I removed names with tRNA, rRNA, fibrillarin, isoprenylcysteine, ubiquinone, regulator, wybutosine and methyltransferase-halogenase (the bifunctional halogenases). That left 113 entries (MEASURED); 108 of them have EC 2.1.1.* and were searched in the larger set below. Tool: `hmmsearch -E 1000 --domE 1000` with `results/hmms/mycosporine/omt.PROVISIONAL.hmm`. Score = best domain score.

| Set | Entries | With a hit at E <= 1000 | Score < 100.6 (no hit counts as below) | Score >= 100.6 |
|---|---|---|---|---|
| First 50 of the filtered list, in UniProt order (the sample asked for; not random) | 50 | 2 (O42898 96.5; B8NI24 44.8) | 50 | 0 |
| The 108 filtered entries with EC 2.1.1.* | 108 | 4 | 107 | 1 (Q8NKC1) |
| All 801 reviewed fungal EC 2.1.1.* entries | 801 | 5 | 800 | 1 (Q8NKC1) |

Hits in the 801-entry set (MEASURED): Q8NKC1 (*Schizosaccharomyces pombe* COMT2) 120.7; O42898 (*S. pombe* COMT1, already an omt decoy) 96.5; A0A397HK53 (ankF, already a decoy) 57.9; B8NI24 (imqG, already a decoy) 44.8; P0DUL6 (*Aspergillus hancockii*, O-methyltransferase hkm8) 42.6.

Reading, for the second round:

- 48 of the 50 sampled fungal O-methyltransferases have no hit at E <= 1000 on the omt HMM. The floor of 100.6 would remove both hits (96.5 and 44.8) of the sample. No sampled entry that is a true target is known, so this does not measure sensitivity for real mycosporine OMTs. It shows the HMM scores almost no fungal O-methyltransferase, so the HMM is specific to a narrow group (and the two references are only 17% identical).
- Q8NKC1 (120.7) would pass the 100.6 floor. The record names it "Probable catechol O-methyltransferase 2" (*S. pombe*, GN SPBPB21E7.04c, annotation only), so it is very probably a non-target that the floor does not remove. It is not a decoy now. Candidate for the decoy set in the second round (PROPOSED).
- The only known fungal target is the Phaffia reference OMT (WEV88680.1), which scores 430.0 to 481.0 on the full HMM; its leave-one-out score is 8.3 or 9.7 (section 3). Whether a fungal ortholog that is far from the Phaffia OMT would score above 100.6 is UNVERIFIED. No such ortholog was tested.

## 10. Not done in section 9

- No BFD search, no tree, no placement.
- No identity between the new decoys and the references was computed. The statement that a decoy is "close in family" rests on Pfam names in the records only.
- No domain cuts. No change to references. The pathway YAMLs got only the new `decoys:` blocks.
- The omt override is a PROPOSED value of the owner, not a result of an analysis of fungal OMTs.

## 11. Step 3 timing run (2026-11-23)

Labels: MEASURED is read from a file or `sacct` named in `analysis/pathway_census/results/step3/melanin_myco/RUN.md`. PROPOSED is an estimate or a choice. UNVERIFIED is not checked.

- One timing job on 200 genomes (chunk 0 of `chunk_manifest.size200`, the chunk of the carotenoid timing job 32375142), 18 HMMs, `omt` at 100.6. Job 32415026 (the first job, 32412179 on `short`, waited about 3 h on the 384-CPU user limit and was cancelled; 32415026 ran on `stajichlab`). Not the full search.
- MEASURED: job 4 min 21 s, script 209.3 s (staging 66.2 s, hmmsearch 143.1 s), MaxRSS 3.1 GB, 2,099,680 proteins in 200 genomes, 0 bad id prefixes, 1,367 candidate rows on 1,299 proteins. No tool change: a combined `min_scores` table and one HMM directory serve the three pathways.
- Candidate rows per step (MEASURED, chunk of 200): tyrosinase 439, hppd 154, omt 125, ayg1 87, atp_grasp_ligase 86, ddgs 81, t4hn_reductase 71, scytalone_dehydratase 64, hexaketide_pks 63, laccase_abr2 63, thn_synthase 53, t3hn_reductase 21, heptaketide_pks 15, hmga 15, asp_melanin_tyrosinase 13, mco_abr1 13, laccase_diphenol 4, aspulvinone_e_synthase 0.
- The `omt` rows all score 100.8 to 160.6 (median 111.4). Section 9.5 found almost no reviewed fungal O-methyltransferase above 100.6. The BFD rows are not reviewed entries. UNVERIFIED what they are. This is evidence for the second round of the `omt` floor.
- PROPOSED plan for the full search: 6 chunks of 4,000 genomes, 16 CPUs, `--mem=96G`, about 64 to 75 min per job, about 159,773 candidate rows (1,367 x 116.88). The default `--mem=32G` is too small (projected 47 to 62 GB, UNVERIFIED). Details in `RUN.md`.

## 12. omt decoy Q8NKC1 and new omt floor (2026-11-23)

Section 9.5 proposed Q8NKC1 as a decoy. It is added. Record fetched from `rest.uniprot.org` (MEASURED): reviewed entry COMT2_SCHPO, "Probable catechol O-methyltransferase 2", *Schizosaccharomyces pombe* (strain 972 / ATCC 24843), ORF SPBPB21E7.04c, 281 residues, EC 2.1.1.6, Pfam PF01596. Step 2 for mycosporine was re-run (job 32415136). Q8NKC1 scores 120.7 on the `omt` HMM, so the highest decoy score is 120.7 (before 100.6, P21964). Reference and leave-one-out scores are unchanged.

The `omt` `min_score` is now 120.7 (basis: owner rule floor = highest decoy score, through `--floor-override`; rule value 8.3). The step 3 timing job used 100.6, and its table is kept in `results/step3/melanin_myco/`. In the timing chunk 84 of 125 `omt` candidate rows lie between 100.6 and 120.7; 41 remain at or above 120.7. Whether the removed rows contain real orthologs is UNVERIFIED. Details in `results/step3/melanin_myco/RUN.md`.
