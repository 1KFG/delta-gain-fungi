---
title: "Genome table prefixes: 5 genomes changed their locus tag"
subtitle: "Cause, how many genomes are affected, the alias fix, and what earlier results are affected"
date: 2026-11-24
---

::: {.callout-important}
## Status and limits

- Labels: MEASURED is read from a file or table named in the text. PROPOSED is a choice. UNVERIFIED is not checked.
- The cause is MEASURED to the level of `samples.csv`: the five genomes were re-annotated from new assembly versions and got new locus tags. Which BFD job or person did this is UNVERIFIED. No BFD file was changed.
- The prefix check reads the first protein id of each file. A file that mixes two prefixes would not be found. For the 5 genomes, all ids were counted (single prefix per file).
- The fix is in our copy only (alias file and corrected table). It maps the new proteins onto the OLD genome table row. The taxonomy is the same species. BUSCO and N50 in that row describe the OLD assembly. For `FF465DA3` the new assembly has a different accession (`GCA_062075423.1`), so the row is not the same genome. The BFD owner must decide the real fix.
- Only `bin/make_calls_table.py` was changed to read the alias file. Other tools that split the protein id at `_` were not changed (section 5).
- Not done: re-checking the 47 files that are newer than the table and have an unchanged prefix (section 3); any biology.
:::

## 1. The five genomes

Genome table: `results/genomes.tsv.zst`, mtime 2026-11-22 12:11 (MEASURED). Protein files are the paths in `seqdb/genomes.tsv` (`genome_annotation/<out>/predict_results/<out>.proteins.fa.gz`).

| Table prefix (old) | New prefix in file | Species | Protein file mtime | Proteins in carotenoid search (old file) | Proteins now (melanin search) |
|---|---|---|---|---|---|
| F8B45790 | F7CB232C | *Blastobotrys sasensis* | 2026-11-24 09:46 | 5,805 | 5,806 |
| FE750E90 | F317F9F5 | *Arxiozyma pintolopesii* RIGP | 2026-11-23 22:27 | 5,034 | 5,021 |
| FE41F6C9 | FB9275E9 | *Fusarium oxysporum* Foc035 | 2026-11-24 10:03 | 16,929 | 17,618 |
| FF465DA3 | F629B7B9 | *Saccharomyces cerevisiae* YJM1417 | 2026-11-24 09:39 | 5,549 | 4,391 |
| F6D07108 | FAD1592B | *Sugiyamaella casensis* | 2026-11-24 10:07 | 6,867 | 6,864 |

(MEASURED: first column pair from `prefix_aliases.tsv`; mtimes from `ls`; counts from `searched.chunk*.tsv` of `carotenoid_full` and `melanin_myco_full`. Three mtimes are later than the two the task named; all five files are newer than the table.) Each current file has exactly one prefix (MEASURED, `zcat | grep '>'`, counts equal the last column).

Where the prefixes are (MEASURED, DuckDB on `tables/*.parquet`, 2026-11-24):

- The old prefixes are in `species.parquet`, `samples.parquet` and `gene_info.parquet` (gene rows: 5,915; 5,246; 17,220; 5,837; 7,010 in the order of the table above). The tables are dated 2026-11-18 (`species`, `samples`) and 2026-11-19 (`gene_info`).
- The new prefixes are in none of these tables.
- The old prefixes are still in `seqdb/genomes.tsv` and in the snapshot `seqdb/proteins/<sha1[:2]>/<old>.proteins.fa.gz` (mtime 2026-11-06). Those old files hold 5,805; 5,034; 16,929; 5,549; 6,867 proteins, equal to the carotenoid search counts. There is no `seqdb/proteins` file for the new prefix `F7CB232C`.
- The current `Fungi_BFD/samples.csv` (mtime 2026-11-21 13:36) lists the new prefixes and not the old ones. Its ASMID for each new prefix differs from the ASMID of the old row: version `.1` became `.2` for four genomes, and `FF465DA3` (`GCA_000977715.4_Sc_YJM1417_v1`) became `GCA_062075423.1_ASM62075423v1` (column `asmid_samples_csv_for_alias` in `prefix_aliases.tsv`).

Cause (MEASURED up to this point, the trigger UNVERIFIED): BFD re-annotated these five genomes from new assembly versions. The run wrote new predict results into the same `genome_annotation/<out>/` folders. The new samples row has a new locus tag. `seqdb/genomes.tsv` was not rebuilt, so it maps the old tag to a file that now holds the new tag. The census read that file path.

## 2. Meaning of "bad id prefix"

`bin/search_chunk.py` counts a protein as bad when its FASTA header does not start with `>` + `<locustag>_` of the manifest row. The manifest row comes from `genomes.tsv.zst`. The search is not stopped by a bad prefix. Every protein is still searched. So the 39,700 proteins (5,806 + 5,021 + 17,618 + 4,391 + 6,864 = 39,700, MEASURED) were searched; only the prefix did not match. The carotenoid search ran when the files still had the old tags (0 bad).

## 3. All genomes

Script: `bin/check_genome_prefixes.py` (tests in `tests/test_check_genome_prefixes.py`). It reads the first protein id of each file for every genome with `has_proteins` true, 8 threads, about 5 min on a login node. Report: `results/genome_prefix_check.tsv`.

MEASURED, 23,376 genomes (all with `has_proteins`):

| Group | Files | Prefix ok | Prefix mismatch |
|---|---|---|---|
| mtime later than 2026-11-22 12:11 (table build) | 52 | 47 | 5 |
| mtime earlier | 23,324 | 23,324 | 0 |

So the 5 genomes are the only mismatches by first id, 5 of 23,376. All 5 are among the 52 files newer than the table.

Related MEASURED finding: the protein count differs between the carotenoid and melanin searches for 50 genomes (the 5 above and 45 more; searched-manifest comparison; totals 244,451,335 and 244,451,875, difference +540). In the 45 others the prefix is unchanged. Examples: `F75533D7` 5,752 to 5,952; `F13F4723` 4,250 to 4,328. The files were rewritten with the same locus tag. For those genomes the protein ids of the carotenoid candidates may name different proteins now. This was not checked (UNVERIFIED). Whether the 45 equal the 47 files "newer and ok" was not checked.

Also observed: two tags of `species.parquet` (`F3547434`, `F6B93D6D`) are not in the current `samples.csv`, but their protein files carry the same tag and are older than the table. Not investigated.

## 4. The fix in our copy

Files:

- `analysis/pathway_census/prefix_aliases.tsv` (tracked; copy in `results/`): 5 rows, columns `locustag` (table), `alias_prefix` (file), species, both ASMIDs, path, mtime, evidence.
- `results/genomes.corrected.tsv.zst` (git-ignored; made by `bin/apply_prefix_aliases.py`): the original table, 24,363 rows, plus the column `alias_prefix` (filled for 5 rows, empty for the rest). The original `results/genomes.tsv.zst` is unchanged. No other column was changed.
- `bin/make_calls_table.py --prefix-aliases prefix_aliases.tsv`: maps a protein id prefix to the table locustag before the genome sets are built. Without the option the behaviour is unchanged.

Rows now mapped (MEASURED, `melanin_myco_full`, 170,712 rows): 13 rows had an unmapped prefix before (FB9275E9 11, F7CB232C 1, FAD1592B 1). With the alias file, 13 of 13 map: 11 rows to `FE41F6C9` (steps `atp_grasp_ligase`, `ayg1`, `ddgs`, `hppd`, `omt`, `tyrosinase`), 1 to `F8B45790` (`hppd`), 1 to `F6D07108` (`hppd`). Rows with the OLD prefixes in the melanin run: 0. `F317F9F5` and `F629B7B9` have 0 rows.

Why `make_calls_table.py` needed it: its `searched` set holds the manifest locustag (old tag), but its hit sets are built from the protein id prefix (new tag). Without the alias, a genome with a hit would be counted as 0 (absent) and the hit prefix would count in `n_call_genomes_not_in_reps`. This is from code reading and a unit test with toy data; the real melanin calls table was not built here (UNVERIFIED on real data).

## 5. Which tools read the prefix

MEASURED by `grep` of `bin/` and `pathway_census/`.

- `bin/search_chunk.py`: only counts `n_bad_id_prefix`. No change needed. The counts stay in the record as measured.
- `bin/make_calls_table.py`: needs the alias file. Done.
- Not changed, same pattern `protein_id.split("_")[0]` (they need `--prefix-aliases` or a shared helper before they run on the melanin or mycosporine candidates): `bin/cutoff_sweep.py`, `bin/placed_vs_hits_by_class.py`, `bin/single_ref_nearest.py`, `bin/hits_by_class_floors.py`, `bin/step5_place.py`.
- `bin/step5_place.py` also reads protein sequences from `seqdb/proteins/<sha1>/<locustag>.proteins.fa.gz` by the prefix of the protein id. For a new prefix there is no such file, so placement of the 13 rows would fail or skip them. UNVERIFIED (not run).
- Gain and loss scripts (`run_gainloss*.R`) read `calls.tsv`; they do not read the prefix.

## 6. Earlier results

- Carotenoid (merged): the candidate rows of the 5 genomes carry the OLD prefixes (MEASURED: F8B45790 4 rows, FE750E90 2, FE41F6C9 16, FF465DA3 3, F6D07108 9; 0 rows with a new prefix; 0 unmapped rows of 238,898). They map through the old genome table without any alias. No carotenoid per-genome table mixes old and new prefixes. The old proteins exist in the 2026-11-06 snapshot `seqdb/proteins`, but not in the current predict files (new numbering), so the carotenoid results for these 5 genomes describe the previous annotation, not the current one.
- Melanin and mycosporine: counts per genome by prefix include 3 unmapped prefixes (as RUN.md states). With the alias file the 13 rows belong to 3 table genomes.
- If the owner rebuilds the BFD tables (new rows for the new prefixes), the carotenoid and melanin results of these 5 genomes will sit on different annotations of the same species. PROPOSED: treat the 5 genomes as one genome with a known annotation change; do not add both.

## 7. Message for the BFD owner

See the end of the task report. Short form: five genomes were re-annotated with new locus tags; `species.parquet`, `gene_info.parquet` and `seqdb/genomes.tsv` still list the old tags; `seqdb/proteins` holds the old proteins. Rebuild the tables and seqdb, or restore the old files. The census will not use the corrected files until you decide.
