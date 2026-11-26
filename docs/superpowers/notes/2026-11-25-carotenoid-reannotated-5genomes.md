---
title: "Carotenoid search repeated on the 5 re-annotated genomes"
subtitle: "Old annotation (merged full search) against the current protein files"
date: 2026-11-25
---

::: {.callout-important}
## Status and limits

- Labels: MEASURED is read from a file named in the text. UNVERIFIED is not checked. All candidate tables are PROVISIONAL (pathway `verified: false`).
- Scope: 5 genomes of 23,376. Only the search-hit level (rows at or above `min_score`) was compared. Tree placement and gain/loss were NOT redone (section 5).
- Cutoffs: the full search used `min_scores.exploratory25.PROVISIONAL.tsv` (x0.75 for `astaxanthin_synthase` and `ylo1`), not the rule table. This repeat uses the same file. MEASURED: the lowest candidate score per step in the full search is at or above the exploratory25 value for all 14 steps, and below the rule value for `astaxanthin_synthase` (1,216.4 against 1,332.5) and `ylo1` (480.5 against 640.7).
- Row matching: new protein ids are a new numbering, so ids cannot be compared. A row is matched by exact protein sequence (old sequence from the 2026-11-06 `seqdb` snapshot, new sequence from the current predict file). A protein with one changed residue counts as lost plus gained.
- `FF465DA3` (YJM1417) is a different assembly (`GCA_000977715.4` became `GCA_062075423.1`). Its comparison is not like for like.
- Which HMM files the full search used was not recorded. MEASURED: the 14 HMM files are dated 2026-11-22 11:03, before the full search (finished 12:55). UNVERIFIED byte identity.
:::

## 1. Run

MEASURED (`results/search/carotenoid_5genomes_reannot/`): 1 job, 32461202, partition `intel`, 16 CPUs, 32 GB, COMPLETED, 51 s elapsed, MaxRSS 282,896 K (`sacct`). Script time 10.4 s. 5 genomes searched, 39,700 proteins, 0 bad id prefixes, 37 candidate rows. Settings, paths and checksums are in the `RUN.md` of that folder.

## 2. Rows per genome (MEASURED, `compare/rows_per_genome_step.tsv`)

Rows are (protein, step) at or above `min_score`. Denominator: rows in the old full search for that genome.

| Genome (old / new prefix) | Proteins old / new | Identical sequences shared | Rows old | Rows new | Unchanged | Lost | Gained |
|---|---|---|---|---|---|---|---|
| F8B45790 / F7CB232C, *Blastobotrys sasensis* | 5,805 / 5,806 | 5,778 | 4 | 4 | 4 | 0 | 0 |
| FE750E90 / F317F9F5, *Arxiozyma pintolopesii* | 5,034 / 5,021 | 4,847 | 2 | 2 | 2 | 0 | 0 |
| FE41F6C9 / FB9275E9, *Fusarium oxysporum* Foc035 | 16,929 / 17,618 | 14,890 | 16 | 17 | 16 | 0 | 1+0... see below |
| FF465DA3 / F629B7B9, *S. cerevisiae* YJM1417 | 5,549 / 4,391 | 1,694 | 3 | 3 | 2 | 1 | 1 |
| F6D07108 / FAD1592B, *Sugiyamaella casensis* | 6,867 / 6,864 | 6,805 | 9 | 10 | 9 | 0 | 1 |
| Total | | | 34 | 37 | 33 | 1 | 4 |

Candidate rows by step (MEASURED, same file). Only pairs with a change are listed. All other (genome, step) pairs have equal row counts and identical sequences.

| Genome | Step | Rows old | Rows new | Change | Score of the changed row | Placement floor of the step |
|---|---|---|---|---|---|---|
| FE41F6C9 | acaa | 2 | 3 | 1 gained | 12.6 | 150 |
| FE41F6C9 | cars | 2 | 3 | 1 gained | 10.0 | 100 |
| F6D07108 | acaa | 1 | 2 | 1 gained | 13.8 | 150 |
| FF465DA3 | crtr | 2 | 2 | 1 lost, 1 gained | 451.3 (old), 451.1 (new) | crtr was not placed |

MEASURED: the 33 unchanged rows have identical scores in both searches (`compare/changed_rows.tsv` has no `score_differs` row). The floors for acaa (150) and cars (100) are in `results/step5/carotenoid/single_ref_steps_RUN.md` of the `pathway-census-realdata` worktree. The three gained acaa and cars rows score 10.0 to 13.8 bits, far under those floors. The README of that worktree states that carx and crtr were not placed. The `crtr` pair has nearly equal scores and 2 of 2 rows remain, so it is probably one locus with a changed sequence. UNVERIFIED: no alignment was run.

## 3. Genome-level presence at the search-hit level

Definition: a genome has the step if it has at least 1 candidate row for it.

MEASURED (`presence_changes` column; 24 genome-step pairs have a row in either search): 0 pairs change. 0 of 5 genomes change in any step. The 4 gained rows and the 1 lost row are in pairs that have another row in both searches. For YJM1417 the comparison is not like for like, and the 0 holds only for its 2 steps (`crtr`, `tdh`).

## 4. Effect on the census calls

- Placed steps of the census (`results/gainloss/carotenoid_noMicro_blasto_5steps/steps.config.tsv`): `phytoene_synthase`, `phytoene_desaturase`, `lycopene_cyclase`, `cart`, `ggpp_synthase`. MEASURED: none of the 5 has a changed row in any of the 5 genomes. FE41F6C9 has 1 row in each of the 5 steps. F6D07108 has 1 row in `ggpp_synthase`. The other 3 genomes have no row in these steps.
- Steps with a changed row: `acaa`, `cars`, `crtr`. The 3 new `acaa` and `cars` rows are under the placement floors, so they would not be placed. `crtr` is not a placed step.
- Because the sequences and scores of all rows in the 5 placed steps are identical, the placement of the new ids is expected to give the same call. UNVERIFIED: placement was not run.

## 5. Not done

- Tree placement and gain/loss were not redone. The number of new proteins that would need placing is 4 (3 under the floors, 1 `crtr`), fewer than 100. Placing them needs `bin/step5_place.py` and the step's tree, and the new proteins in `seqdb/proteins/<sha1[:2]>/<new prefix>.proteins.fa.gz`. Those files do not exist for the new prefixes (`2026-11-24-genome-prefix-mismatch.md`). The run is not simple, so it was not done.
- The 45 other genomes whose protein file changed under the same prefix (prefix-mismatch note, section 3) were not searched.
- Other pathways for these 5 genomes were not compared. The melanin and mycosporine search already used the current files.
- BUSCO and N50 of the new assemblies.

## 6. Negative results

- 0 of 5 genomes change step presence at the search-hit level.
- 0 changed rows in the 5 placed steps.
- 1 of 34 old rows is lost (and has a near-identical replacement). 33 of 34 are unchanged by sequence and score. 4 rows are gained, all in steps `acaa`, `cars`, `crtr`.
