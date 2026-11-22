---
title: "Census step list for review (DRAFT)"
subtitle: "BFD pathway census: carotenoid, DHN melanin, L-DOPA melanin, mycosporine"
date: 2026-11-22
---

::: {.callout-important}
## Status and limits

- DRAFT for the owner's review. Nothing here has been run on BFD.
- It follows the checklist in `2026-11-20-pilot-lessons.md` section 4 and the plan `2026-11-21-pathway-census-carotenoid-real-data.md` (tasks R1 to R10). It extends them to four pathways.
- Sizes and run times are UNVERIFIED until step 3 times one real job.
- Owner decisions of 2026-11-22 apply: strict reading of rule 1; gene-name rule accepted; names in `notes/2026-11-22-owner-decisions.md`. Other decisions are deferred and the current tables stand.
- Each output carries the `.PROVISIONAL` stamp until the owner sets `verified: true` per pathway (step 12).
:::

## Inputs (MEASURED from the repository, not re-run)

| Pathway | Reference sequences (`refs/<id>/*.faa` records) |
|---|---|
| carotenoid | 29 |
| melanin_dhn | 22 |
| melanin_ldopa | 15 |
| mycosporine | 8 |

- BFD genome table: 24,363 genomes. 23,379 have a protein file and a genome FASTA. 984 have neither (README, MEASURED 2026-11-20).
- Genomes without proteins are never read as absences.
- The carotenoid, DHN and L-DOPA tables use the strict reading. All routes end at T3 under it.

## Rules that apply to every step

1. A step call needs HMM evidence and a tree placement. A bare HMM or Pfam hit does not count.
2. NA is never "clean". Use `clean_mask` and `unknown_counts`.
3. A copy is a gene. Isoforms are collapsed.
4. Every number is labelled MEASURED, PROPOSED or UNVERIFIED.
5. No biological finding is written before a result file supports it.
6. Python: `/usr/bin/python3.12`. SLURM: no `BASH_SOURCE`. Temporary files in `${SCRATCH:?}`. Results to `/bigdata`.
7. Large tables are compressed with `zstd`.

## Steps

Each step lists what it does, what it needs and what you check. "Gate" marks a stop for your sign-off.

| # | Step | Needs | Output | You check |
|---|---|---|---|---|
| 0 | **Merge `pathway-census-realdata` to `main`, or run from the worktree.** | your decision | one code base | Gate. Merging is your call. |
| 1 | **Build per-step HMMs** (`build_step_hmms.py`) for each pathway, from the references. Fused proteins are cut per domain. | 0 | HMMs per step | HMM count per pathway; the carotenoid synthase and cyclase HMMs differ. |
| 2 | **Set `min_score` per step.** PROPOSED rule: lowest reference score. Run a leave-one-out check, because references train the HMM. Run each decoy against the HMMs. | 1 | table of `min_score` and decoy scores | Gate. Any decoy that scores above `min_score` is listed. |
| 3 | **Time one real search job** on a small genome chunk (batched driver, many genomes per job, all step HMMs). Then size the rest toward 1 to 1.5 h per job. | 1, 2 | timing table, chunk manifest | Gate. Job count and total hours. |
| 4 | **Full BFD search.** One searched-genome manifest lists every genome searched. Dereplicate at 90% identity only if a family has more than 50,000 proteins. | 3 | candidate tables (`.zst`), manifest | Searched genomes versus 23,379. |
| 5 | **Reference trees and placement.** One tree per step: references, decoys as outgroup, candidates. Tool: iqtree or raxml-ng for trees, `epa-ng` through pixi for placement. | 4 | placements TSV | Gate. Decoys must place outside the target clade. |
| 6 | **Assign steps** (`assign_steps.py`). | 4, 5 | step labels per protein | Counts of assigned, unassigned, fusion candidates. |
| 7 | **Quality flags.** Needs a fungal versus non-fungal reference set for contamination and an all-versus-all search for haplotype duplicates. Neither exists yet. | 6 | flags table | Gate. Choose the contamination reference set, or run without it (see Decisions). |
| 8 | **Copy matrix and phylum summary** (`build_matrix.py`). Raw and filtered counts side by side. | 6, 7 | genome by step matrix, phylum summary | Genomes with NA flags. |
| 9 | **Deposited arm.** Run steps 4 to 8 on deposited proteins of the 5,634 genomes with a deposited protein file. Collapse isoforms by gene from the GFF. | 8 | deposited matrix | ID mapping rule. |
| 10 | **Paired comparison** (`paired_compare.py`), deposited versus BFD on the same assemblies. | 8, 9 | paired table | Pair count. |
| 11 | **Absence check** (`absence_check.py`, tblastn on SLURM) with per-step `min_bitscore` calibrated on decoys. Merge with `merge_absence.py`. | 8 | absence summary | Outcome names (`tblastn_hit` rather than `present_unannotated`). |
| 12 | **Triage and report.** Run `triage.py`. Render chapters. Set `verified: true` per pathway only after your review. | 8, 10, 11 | chapters, tables, figures | Gate. `verified` per pathway. |

## Decisions I need from you

1. **Merge or work in the worktree** (step 0). Recommendation: work in the worktree until step 12.
2. **Contamination reference set** (step 7). The set does not exist. Options: (a) build it from the SwissProt table with a taxid-to-lineage map (UNVERIFIED whether the table is usable); (b) run first without it and report raw counts only. Recommendation: (b) first, so the matrix is not blocked.
3. **Order of pathways.** Option A: carotenoid first as the pilot, then the other three. Option B: all four through step 6 together. Recommendation: A, because the step list for the three other pathways reuses the same code and a first full run exposes defects at lower cost.
4. **Tyrosinase and laccase.** Run as separate trees with the decoys inside (see the owner note). Confirm.
5. **Ortholog naming** is not part of this list. It is a separate task after step 12.

## Known state of side work (MEASURED with `sacct`, 2026-11-22)

- Mycosporine diamond search, job 32359412: TIMEOUT after 02:00:18. The output file `diamond_hits.tsv` exists. I did not check whether it is complete.
- Mycosporine placement, job 32359421: CANCELLED after 0 s of run time.
- The mycosporine pilot should be re-run as part of steps 4 and 5, not repaired on its own.

## Step 3 result (carotenoid, 2026-11-22, MEASURED)

- Job 32375142, short queue, 16 CPUs, 1 chunk of 200 genomes (seeded shuffle), 14 step HMMs: job time 1 min 31 s, script time 78.9 s, peak memory 1.48 GB.
- Time split: staging and counting proteins 20.8 s; 14 hmmsearch runs 3.4 to 7.2 s each.
- 200 genomes, 2,099,659 proteins, 1,827 candidates (protein by step rows) at the accepted `min_score` values. All 200 genomes searched, 0 bad id prefixes.
- Genomes with protein files: 23,376 listed in the chunk manifest; 987 excluded (984 with no seqdb path, 3 with a seqdb path and a missing file). The README of 2026-11-20 gave 984 excluded.
- Tools: `bin/make_min_scores.py`, `bin/make_chunks.py`, `bin/search_chunk.py`, `slurm/search_chunk.sbatch`. No unit tests yet for these three scripts.
