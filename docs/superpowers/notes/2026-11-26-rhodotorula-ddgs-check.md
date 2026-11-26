---
title: "Rhodotorula and ddgs: below the cutoff, or nothing close?"
subtitle: "Best ddgs HMM hit per genome with no score cutoff, mycosporine census"
date: 2026-11-26
---

::: {.callout-important}
## Status and limits

- Labels: MEASURED = read from a file named in the text. UNVERIFIED = not checked. All outputs are PROVISIONAL (mycosporine `verified: false`).
- This is a one-HMM search score check. A bit score is not a call. Census calls also need tree placement (step 5). Nothing here changes a census result.
- Search threshold was hmmsearch `-E 10`, not a score of 0. No genome in this check lacked a reported hit (0 of 491), so the threshold did not hide a genome. A weaker hit could still exist below E 10 (UNVERIFIED).
- No biology claim. The identity of the best Rhodotorula hit protein was NOT checked (UNVERIFIED).
- Pyricularia is the first 100 of 521 genomes with proteins by sorted locustag, not a random sample.
- Scores come from one HMM (`ddgs.PROVISIONAL.hmm`, 412 positions). Other ddgs-like HMMs, and a search of the DNA for unannotated genes, were not done.
:::

## 1. Answer

MEASURED: in all 215 Rhodotorula genomes the best ddgs domain score is 68.5 to 109.8 bits (median 102.2). The cutoff is 463.5 and the highest decoy is 315.5. The best Rhodotorula score is 353.7 bits below the cutoff. No genome is between 150 and 463.5. So the genomes do not score just below the cutoff. Nothing close exists in the predicted proteins, for this HMM.

## 2. Inputs (MEASURED)

| Item | Value | Source |
|---|---|---|
| Cutoff for `ddgs` | 463.5 (lowest leave-one-out score, 3 refs) | `results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv` (confirmed) |
| Max decoy score | 315.5 (reference range 779.4 to 827.6) | `results/step2/mycosporine/step_scores.PROVISIONAL.tsv` |
| Rhodotorula genomes | 220 in `genomes.corrected.tsv.zst`; 215 with protein files; 5 without | `results/genomes.corrected.tsv.zst` (tables worktree) |
| Cross-check | 215 = `n` for Rhodotorula in `reports/mycosporine/figures/out/genus_check.tsv`. That table shows ddgs_placed 0 and omt_hit 19 | same file |
| Prefix aliases | None of the 491 searched genomes is among the 5 re-annotated genomes (those are Sugiyamaella, Blastobotrys, Fusarium, Arxiozyma, Saccharomyces) | `selection_notes.txt` |

## 3. Distribution of best score (MEASURED, `bins.tsv`)

Denominator: genomes with protein files in the group. Bins in bit score.

| Group | n | <50 | 50-150 | 150-300 | 300-463.5 | >=463.5 | >=315.5 (max decoy) | median | max |
|---|---|---|---|---|---|---|---|---|---|
| Rhodotorula | 215 | 0 | 215 | 0 | 0 | 0 | 0 | 102.2 | 109.8 |
| Phaffia (positive control) | 12 | 0 | 5 | 0 | 0 | 7 | 7 | 743.2 | 776.3 |
| Cystobasidium | 16 | 0 | 0 | 0 | 0 | 16 | 16 | 604.2 | 608.4 |
| Cystofilobasidium | 8 | 0 | 8 | 0 | 0 | 0 | 0 | 104.6 | 106.2 |
| Rhodosporidiobolus | 14 | 0 | 14 | 0 | 0 | 0 | 0 | 106.9 | 109.1 |
| Sporobolomyces | 22 | 0 | 22 | 0 | 0 | 0 | 0 | 103.5 | 105.7 |
| Naganishia | 28 | 0 | 1 | 0 | 0 | 27 | 27 | 563.2 | 569.8 |
| Hortaea | 76 | 0 | 1 (+1 no hit) | 8 | 4 | 62 | 66 | 610.4 | 614.7 |
| Pyricularia (first 100) | 100 | 0 | 3 | 0 | 0 | 97 | 97 | 607.3 | 614.7 |

Notes on the table:

- Hortaea has 1 genome with no reported hit at E 10 (counted in no bin). The Hortaea row sums to 75 plus this 1 = 76.
- All 215 Rhodotorula scores are in 50 to 150. The lowest is 68.5 (*R. toruloides*).
- Coverage of the HMM by the best Rhodotorula hit: 0.568 to 0.636 (median 0.626) of 412 positions. Genomes with a ddgs call have median coverage of 0.98 to 0.99 (Cystobasidium, Phaffia, Naganishia, Pyricularia, Hortaea). MEASURED, `per_genome.tsv`.
- The same pattern occurs in Phaffia genomes without a ddgs call: 5 of 12 score 94.3 at coverage 0.580. The 7 others score 742.7 to 776.3 at coverage 0.983 to 0.995.

## 4. Cross-checks (MEASURED)

- The set of genomes at or above 463.5 here equals the set of genomes with a `ddgs` row in the step 3 candidate files (`results/search/melanin_myco_full/candidates.chunk*.zst`) in all 491 genomes (0 disagreements). The cutoff and HMM used here match the full search.
- The Phaffia count at or above 463.5 (7 of 12) and the Cystobasidium count (16 of 16) match `ddgs_placed` in `genus_check.tsv`.

## 5. Sensitivity table (not a result; MEASURED counts of genomes with best score at or above a lowered cutoff)

Counts genomes with a domain hit at that score. These are not census calls: step 5 placement was not run on them.

| Group | n | 463.5 | -10% (417.15) | -25% (347.62) | -50% (231.75) |
|---|---|---|---|---|---|
| Rhodotorula | 215 | 0 | 0 | 0 | 0 |
| Phaffia | 12 | 7 | 7 | 7 | 7 |
| Cystobasidium | 16 | 16 | 16 | 16 | 16 |
| Naganishia | 28 | 27 | 27 | 27 | 27 |
| Pyricularia (100) | 100 | 97 | 97 | 97 | 97 |
| Hortaea | 76 | 62 | 66 | 66 | 72 |

Lowering the cutoff by 10, 25 or 50 percent adds 0 Rhodotorula genomes. The best Rhodotorula score (109.8) is 23.7 percent of the cutoff, so a cutoff 76 percent lower would be needed to reach it. The full table for the other genera is `sensitivity.tsv`. Hortaea shows that a lowered cutoff can add genomes in a genus where scores lie below the cutoff (10 percent: +4; 50 percent: +10). That has not been examined here.

## 6. What the numbers show and do not show

Show:

- For this HMM, no Rhodotorula protein scores near the ddgs range. The best hits are at 63 percent HMM coverage and about 100 bits.
- The same low-score, partial-coverage hit appears in genera without a ddgs call (Cystofilobasidium, Rhodosporidiobolus, Sporobolomyces) and in 5 Phaffia genomes. This points to a background hit shared by genomes that lack a ddgs-like protein. The protein was not identified (UNVERIFIED).

Do not show:

- That Rhodotorula has no ddgs gene. Only predicted protein files (215 of 220 genomes) were searched. A gene missed by gene prediction, or on an assembly gap, would not appear. Assembly completeness (BUSCO) was not used.
- That Rhodotorula cannot make mycosporines. No biology was tested.
- The 5 Rhodotorula genomes without protein files were not searched (5 of 220).

## 7. Negative results

- No Rhodotorula genome is near the cutoff. The "just below" explanation is rejected for this HMM.
- Lowering the cutoff does not rescue any Rhodotorula call.

## 8. Not done here

- No tblastn or HMM search of the Rhodotorula assemblies for unpredicted genes.
- No identification of the 100-bit best hit (family, function).
- No search with a second ddgs HMM built from other reference sequences.
- No step 5 placement of any sub-cutoff hit.
- No change to any census table or cutoff.

## 9. Reproduction

`analysis/pathway_census/results/search/mycosporine_rhodotorula_ddgs/RUN.md`. Script `curation_scripts/mycosporine/rhodotorula_ddgs_check.py`. Tests `tests/test_rhodotorula_ddgs.py` (3 passed). Run time 30.8 s on 4 CPUs.
