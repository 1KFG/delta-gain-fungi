---
title: "Census step 4: full BFD search for DHN melanin, L-DOPA melanin, mycosporine"
subtitle: "Candidate counts per step, checks, and open decisions for step 5"
date: 2026-11-24
---

::: {.callout-important}
## Status and limits

- All results are PROVISIONAL. The three pathway files have `verified: false`.
- The rows are HMM hits at or above the step `min_score`. They are candidates, not orthologs. Placement (step 5) was not done. No function is claimed.
- One check failed: 39,700 proteins in 5 genomes have an id prefix that is not the locus tag of `genomes.tsv.zst` (section 2). The search of those genomes is complete. 13 rows cannot be mapped to a genome row.
- The jobs ran on the `intel` partition with `--time=04:00:00`, not on `stajichlab` (section 1). Job `sacct` records were not reachable after the run, so job State, Elapsed and MaxRSS from `sacct` are UNVERIFIED.
- The `omt` floor is 100.6 (owner decision 2026-11-24). Candidates between 100.6 and 120.7 stay in the table; they are a separate band in this note.
- Labels: MEASURED is read from a file named in the text. PROPOSED is a choice or an estimate. UNVERIFIED is not checked.
:::

Record with commands and job table: `analysis/pathway_census/results/search/melanin_myco_full/RUN.md`.

## 1. What was run

- 18 HMMs, one combined cutoff table (sha256 `9b939123...b8fec`, equal to the step 3 record), 6 chunks of 4,000 genomes (the last has 3,376), 16 CPUs, `--mem=96G`.
- First submission on `stajichlab` (jobs 32427757 to 32427762) stayed PENDING 21 min with reason Priority. The scheduler estimate was 2026-11-30. I cancelled them and resubmitted on `intel` (jobs 32427878 to 32427883, `--time=04:00:00`). All six started within 20 s.
- Script time per chunk: 5,224 to 6,438 s (MEASURED). The step 3 projection was 63 to 75 min (PROPOSED, from a ryzen node). The observed time is about 1.5 times that. No OOM, no TIMEOUT. The 4 h limit was not needed; 2 h would have been enough (longest 107 min).

## 2. Checks

| Check | Result (MEASURED) |
|---|---|
| Searched lists | 6 files; union 23,376 genomes; 0 duplicates; 0 missing; 0 extra against `has_proteins` in `genomes.tsv.zst` |
| Proteins searched | 244,451,875 (carotenoid full search: 244,451,335; difference +540) |
| Bad id prefixes | 39,700 in 5 genomes. Expected 0. FAIL. |

The 5 genomes: F8B45790, FE750E90, FE41F6C9, FF465DA3, F6D07108. Every protein in each file has a new prefix (F7CB232C, F317F9F5, FB9275E9, F629B7B9, FAD1592B). The carotenoid search had 0 bad prefixes for the same genomes. The two protein files I checked were rewritten on 2026-11-24 at 09:46 and 10:03, after the genome table was built. UNVERIFIED: the cause (a re-annotation run of the BFD pipeline is my reading; I did not check). The change in protein counts versus the carotenoid run is -1,158 to +689 per genome (5,805 to 5,806; 5,034 to 5,021; 16,929 to 17,618; 5,549 to 4,391; 6,867 to 6,864).

## 3. Candidates per step

Rows = (protein, step) pairs at or above `min_score`. Projection = the step 3 value (1,367 rows in 200 genomes x 116.88, PROPOSED). Ratio = observed / projection.

| Pathway | Step | min_score | Rows | Unique proteins | Genomes with a row | Max score | Projection | Ratio |
|---|---|---|---|---|---|---|---|---|
| melanin_dhn | ayg1 | 374.2 | 11,836 | 11,836 | 8,526 | 807.6 | 10,169 | 1.16 |
| melanin_dhn | t4hn_reductase | 356.0 | 10,010 | 10,010 | 8,092 | 523.6 | 8,298 | 1.21 |
| melanin_dhn | laccase_abr2 | 464.3 | 9,631 | 9,631 | 4,923 | 1,455.7 | 7,363 | 1.31 |
| melanin_dhn | hexaketide_pks | 1839.0 | 8,604 | 8,604 | 7,819 | 5,176.4 | 7,363 | 1.17 |
| melanin_dhn | scytalone_dehydratase | 240.5 | 8,604 | 8,604 | 7,511 | 336.3 | 7,480 | 1.15 |
| melanin_dhn | thn_synthase | 1908.6 | 7,286 | 7,286 | 6,535 | 5,334.6 | 6,195 | 1.18 |
| melanin_dhn | t3hn_reductase | 550.4 | 2,240 | 2,240 | 2,203 | 619.2 | 2,454 | 0.91 |
| melanin_dhn | heptaketide_pks | 3656.7 | 2,200 | 2,200 | 2,168 | 4,853.9 | 1,753 | 1.25 |
| melanin_dhn | mco_abr1 | 638.2 | 1,918 | 1,918 | 1,897 | 1,394.9 | 1,519 | 1.26 |
| melanin_ldopa | tyrosinase | 63.6 | 50,452 | 50,452 | 14,220 | 838.2 | 51,310 | 0.98 |
| melanin_ldopa | hppd | 164.3 | 20,103 | 20,103 | 16,847 | 974.9 | 18,000 | 1.12 |
| melanin_ldopa | asp_melanin_tyrosinase | 320.5 | 2,471 | 2,471 | 1,847 | 878.6 | 1,519 | 1.63 |
| melanin_ldopa | hmga | 983.2 | 2,313 | 2,313 | 2,280 | 1,078.5 | 1,753 | 1.32 |
| melanin_ldopa | laccase_diphenol | 1127.7 | 402 | 402 | 215 | 1,437.6 | 468 | 0.86 |
| melanin_ldopa | aspulvinone_e_synthase | 1205.5 | 168 | 168 | 127 | 2,240.8 | 0 | not defined |
| mycosporine | omt | 100.6 | 14,441 | 14,441 | 12,097 | 427.1 | 14,610 | 0.99 |
| mycosporine | atp_grasp_ligase | 105.0 | 9,167 | 9,167 | 8,859 | 1,038.1 | 10,052 | 0.91 |
| mycosporine | ddgs | 463.5 | 8,866 | 8,866 | 8,762 | 776.3 | 9,467 | 0.94 |
| All | | | 170,712 | 161,355 | 19,645 | | 159,773 | 1.07 |

Per pathway (MEASURED): melanin_dhn 62,329 rows, 53,024 unique proteins, 11,254 genomes. melanin_ldopa 75,909 rows, 75,909 unique proteins, 19,073 genomes. mycosporine 32,474 rows, 32,422 unique proteins, 12,739 genomes. The genome counts use the id prefix and include the 3 unmapped prefixes (section 2).
Reading: 15 of 18 steps are within 0.86 to 1.32 of the projection. `asp_melanin_tyrosinase` (1.63) is the largest gap. The scaling from 200 genomes has a limit: one 200-genome chunk gave 2 to 54 rows for the small steps.

## 4. What is notable (counts only)

- `aspulvinone_e_synthase`: 0 rows in the 200-genome chunk and 168 rows in the full search (127 genomes). Scores 1,205.7 to 2,240.8, median 1,298.1. The cutoff is 1,205.5. By taxonomy table: 167 rows in Eurotiomycetes (152 Aspergillus, 14 Penicillium, 1 *Aspergillaceae* sp.) and 1 in Sordariomycetes (*Bertia*). The cutoff is 54% of the reference score (step 3 record). Most rows are near the cutoff.
- `asp_melanin_tyrosinase` (1.63): 2,471 rows. 2,359 in Eurotiomycetes (1,528 Aspergillus, 662 Penicillium), 60 Sordariomycetes. 0 proteins hit both this HMM and `tyrosinase`.
- `laccase_diphenol`: 402 rows, all in Cryptococcus (class Tremellomycetes), 215 genomes. Max score 1,437.6.
- `tyrosinase`: 50,452 rows in 14,220 genomes. Rows per genome: 1 in 2,544 genomes; 2 in 2,455; 3 in 3,337; 4 in 2,049; 5 to 9 in 3,660 (the sum of 5, 6, 7, 8, 9: 1,364 + 1,034 + 703 + 352 + 207); 10 or more in 175 genomes; maximum 28. 3,475 rows (6.9%) score below 127.2 (2 x the cutoff). Median score 320.1. Top genera by rows: Fusarium 11,481; Aspergillus 2,951; Colletotrichum 2,730.
- `omt`: 14,441 rows (ratio 0.99 to the 100.6 projection). Split at the old floor: 10,097 rows in 100.6 <= score < 120.7 (8,829 genome prefixes) and 4,344 rows at 120.7 or above (4,282). Score median in the low band 109.3. In the high band 10 rows score 200 or more; maximum 427.1. The low band by phylum: Ascomycota 9,783; Basidiomycota 169; Mucoromycota 72; Glomeromycota 36. By class: Eurotiomycetes 3,345; Dothideomycetes 2,555; Sordariomycetes 2,122; Pichiomycetes 859; Lecanoromycetes 381. Top genera: Aspergillus 1,953; Fusarium 1,064; Penicillium 791; Macrophomina 491; Alternaria 275. The high band: Sordariomycetes 3,723 of 4,344; Fusarium 1,913; Pyricularia 446; Beauveria 369. The low band spreads over the genera that have many genomes in BFD (Sordariomycetes are 5,981 of the 23,376 genomes; Eurotiomycetes 3,201). I did not test whether the band follows genome count. What these proteins are is not known.
- PKS overlap (proteins that hit at least one of `heptaketide_pks`, `hexaketide_pks`, `thn_synthase`): 9,374 proteins. 994 hit one HMM; 8,044 hit two; 336 hit all three. Pairs: hexaketide and thn_synthase 6,316; hexaketide and heptaketide 1,648; heptaketide and thn_synthase 80. The PKS rows (18,090) lie on 9,374 proteins. The HMMs do not separate the three PKS steps (known from step 2).
- SDR: 583 proteins hit both `t3hn_reductase` and `t4hn_reductase` (2,240 and 10,010 rows). Laccase and MCO: 4 proteins hit two of `laccase_abr2`, `laccase_diphenol`, `mco_abr1`. Tyrosinase: the same 4 proteins are the only ones that hit `tyrosinase` plus a laccase or MCO step.
- Steps with 0 rows: none. The smallest are `aspulvinone_e_synthase` (168) and `laccase_diphenol` (402).
- Genomes with at least one candidate in any step: 19,645 of 23,376 (84%). Many HMMs in the set are broad families (tyrosinase, OMT, hppd), so this number is not a pathway presence number.

## 5. Not done here

- No placement, no tree, no biology claim, no rule on pathway completeness.
- No change to cutoffs. No second round on `omt`.
- The 5 genomes with changed locus tags were not re-mapped.
- The node speed difference (about 1.5 times slower than the projection) was not tested on a ryzen node.
- `sacct` records (State, Elapsed, MaxRSS, ExitCode) were not read (slurmdb timeout).

## 6. Open decisions for step 5 (owner)

All are PROPOSED. None has data beyond the counts above.

1. **Genome table for the 5 re-annotated genomes.** Rebuild `genomes.tsv.zst` rows for F8B45790, FE750E90, FE41F6C9, FF465DA3, F6D07108 with the new locus tags, or drop their 13 candidate rows. Needed before per-genome or taxon tables.
2. **Decoy or outgroup curation.** Steps with a very large count and a wide score range: `tyrosinase` (50,452 rows; cutoff set by one reference; 6.9% of rows below 2 x cutoff), `hppd` (20,103), `omt` (14,441; 10,097 in the band above the decoy floor), `ayg1`, `t4hn_reductase`, `laccase_abr2`, `atp_grasp_ligase`, `ddgs`. These need outgroup sequences (related families) in the trees. Steps with little data on separation: `scytalone_dehydratase`, `hmga`, `atp_grasp_ligase` have decoys that score nothing on the step HMM (step 2 note).
3. **PKS joint tree.** 8,380 of 9,374 PKS proteins hit two or three PKS HMMs. One joint tree of the three PKS steps (with outgroup PKS of known other function) is the way to assign steps. Per-step trees would repeat the same proteins.
4. **SDR joint tree.** 583 proteins hit both `t3hn_reductase` and `t4hn_reductase`. Place t4HN and t3HN reductases on one tree.
5. **Tyrosinase and laccase trees.** The 4 proteins that hit both `tyrosinase` and a laccase or MCO HMM are few. Separate trees for tyrosinase, for laccase (`laccase_abr2`, `laccase_diphenol`) and for `mco_abr1` are possible. Decide whether `asp_melanin_tyrosinase` (2,471 rows) goes into the tyrosinase tree as a clade or alone.
6. **Subsampling.** 50,452 tyrosinase proteins and 20,103 hppd proteins are too many for one tree. A clustering step (identity threshold) before the tree is needed. The threshold is not chosen.
7. **`omt` low band.** Keep the 10,097 rows (100.6 to 120.7) in step 5, or drop them (floor 120.7, as the earlier recommendation). Only a placement against outgroup COMTs can show what they are.
