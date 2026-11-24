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
- Decoys are few. Four steps have no decoy and one has decoys that score nothing (section 4). A step without a decoy score has no check of separation.
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
| melanin_dhn | 9 | 9 | 22 | 31 |
| melanin_ldopa | 6 | 6 | 15 | 14 |
| mycosporine | 5 | 3 | 8 | 16 |

Reference file problems (MEASURED):

- **Steps without a reference.** `mycosporine`: `mg_glucosyltransferase` and `glutamicol_formation` have `refs: []` and no `.faa`. No HMM is possible. `build_step_hmms.py` stops with a mafft error on the first such step unless `--steps` is given. I built the three steps with `--steps ddgs omt atp_grasp_ligase`.
- **Non-amino-acid characters.** `melanin_dhn/hexaketide_pks.faa` (Q9Y7A7) has 37 `X` characters in 2,177 residues. Other files have none. hmmbuild accepts `X`. Effect on the HMM: not checked.
- **Duplicate ids, empty sequences, count mismatches.** None in the 45 reference files or in the decoy files.
- **Same decoy in several files.** A decoy that applies to several steps is written in each step's decoy file (for example Q03149 in three PKS decoy files, P17489 in two laccase files). The scorer treats them as separate targets. Counts in this note use unique accessions.
- **Same protein as reference and decoy.** Q55P57, P38993 and Q4WZB4 are references in one step and decoys in another (by design).
- **Decoy folders missing.** No decoy file for `ayg1`, `scytalone_dehydratase`, `hppd`, `hmga`, `aspulvinone_e_synthase` (section 6).

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
| ayg1 | 3 | 406 | 747.3 to 807.6 | 374.2 | 374.2 | LOO | none (no decoy file) | 0 | 0 | LOO of Q6BCB5 (Exophiala) is 374.2; the other two are 615.9 and 632.6. No decoy, so separation is untested. |
| t4hn_reductase | 5 | 275 | 462.2 to 523.6 | 356.0 | 356.0 | LOO | 342.9 (Q00791 stcU) | 0 | 1 (A7DWM6 t3hn 357.1) | **Margin of 13.1 bits to the best decoy.** A t3hn reference scores above `min_score`. |
| t3hn_reductase | 3 | 282 | 614.7 to 619.2 | 550.4 | 550.4 | LOO | 414.6 (Q00791 stcU) | 0 | 0 | Margin 135.8 bits. LOO range is narrow (550.4 to 570.7). Rule value is 90% of the lowest full-HMM score (candidate for lowering, section 7). |
| scytalone_dehydratase | 4 | 168 | 321.3 to 336.3 | 240.5 | 240.5 | LOO | none (no decoy file) | 0 | 0 | No decoy, so separation is untested. |
| laccase_abr2 | 1 | 587 | 1455.7 | single reference | 464.3 | floor | 464.3 (P17489 yA) | 1 | 0 | Floor set by one laccase decoy. Next decoy 98.2. Floor is 32% of the reference score. |
| mco_abr1 | 2 | 686 | 1213.3 to 1446.5 | 638.2 | 638.2 | LOO | 385.5 (P38993 FET3) | 0 | 0 | Margin 252.7 bits. |

### melanin_ldopa

| Step | Refs | HMM length | Ref range | LOO min | min_score | Basis | Highest decoy | Decoys >= min | Other refs >= min | Flag |
|---|---|---|---|---|---|---|---|---|---|---|
| tyrosinase | 7 | 563 | 394.7 to 839.8 | 63.6 | 63.6 | LOO | 36.2 (B1VTI5 grixazone synthase) | 0 | 0 | **LOO min is far below the others.** Q00234 (melO) gives 63.6. The next lowest LOO is 239.8 (A0A1S9DK56), then 377.5 (P00440). Margin to the best decoy is 27.4 bits. |
| laccase_diphenol | 3 | 624 | 1337.1 to 1437.6 | 1127.7 | 1127.7 | LOO | 180.5 (Q8LPL3 ascorbate oxidase) | 0 | 0 | Margin 947.2 bits. Rule value is 84% of the lowest full-HMM score (candidate for lowering). |
| hppd | 1 | 403 | 974.9 | single reference | **undefined** | floor needs a decoy | none (no decoy file) | n/a | n/a | `make_min_scores.py` stops with an error. No value was invented. |
| hmga | 2 | 450 | 1069.3 to 1078.5 | 983.2 | 983.2 | LOO | none (no decoy file) | 0 | 0 | Rule value is 92% of the lowest full-HMM score. The two references are close relatives (LOO 983.2 and 983.4). No decoy. |
| aspulvinone_e_synthase | 1 | 925 | 2248.1 | single reference | **undefined** | floor needs a decoy | none (no decoy file) | n/a | n/a | Same as `hppd`. |
| asp_melanin_tyrosinase | 1 | 356 | 878.6 | single reference | 320.5 | floor | 320.5 (Q5BGU9 ivoB) | 1 | 0 | Floor set by ivoB, a tyrosinase-family decoy. Next decoy 283.8 (Q2UNF9). The tyrosinase references score 7.8 to 21.0 on this HMM. |

For the two undefined steps the `make_min_scores.py` run on the full `step_scores` table fails (MEASURED, `ValueError: step aspulvinone_e_synthase: no value for rule basis`). The `min_scores` table for melanin_ldopa was made from `step_scores.floor_defined.PROVISIONAL.tsv`, which is `step_scores` without those two rows.

### mycosporine

| Step | Refs | HMM length | Ref range | LOO min | min_score | Basis | Highest decoy | Decoys >= min | Other refs >= min | Flag |
|---|---|---|---|---|---|---|---|---|---|---|
| ddgs | 3 | 412 | 779.4 to 827.6 | 463.5 | 463.5 | LOO | 315.5 (Q9ZAE9 acbC) | 0 | 0 | Margin 148.0 bits. The LOO min is the Phaffia reference (WEV88679.1). The two cyanobacterial references are 73% identical (yaml note), so their LOO scores (693.2, 676.6) are not independent. |
| omt | 2 | 335 | 430.0 to 481.0 | 8.3 | 8.3 | LOO | 100.6 (P21964 human COMT) | 4 | 0 | **Rule value is below all four decoys.** Each of the two references scores 8.3 or 9.7 when it is held out. The two references are 17% identical (yaml note), so the HMM of one does not detect the other. |
| atp_grasp_ligase | 3 | 500 | 688.7 to 1038.1 | 105.0 | 105.0 | LOO | none (6 decoys, none has a hit at E <= 10) | 0 | 0 | The LOO min is the cyanobacterial reference B2J6X7. The two Phaffia references are 98% identical (490 of 500 residues, `UNVERIFIED.md`), so their LOO scores (955.8, 961.2) are not independent. |

## 4. Steps without a decoy (floor undefined or untested)

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

Nothing is applied. No value was lowered. The `--scale` option of `make_min_scores.py` was not used.

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
