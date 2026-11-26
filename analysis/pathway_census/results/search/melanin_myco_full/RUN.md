# Census step 4 for melanin_dhn, melanin_ldopa, mycosporine: full BFD search

Date: 2026-11-24. Branch `pathway-census-melanin-myco`. All outputs are PROVISIONAL (the three pathway files have `verified: false`).
Labels: MEASURED = read from a file named here. PROPOSED = a choice or estimate. UNVERIFIED = not checked.
Placement (step 5) was not done. No biology claim is made here.

## Inputs

| Item | Value |
|---|---|
| Cutoff table | `results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv`. sha256 `9b9391230fb5e4d28d2b7aa9ca49f7f4d81388ec8e506662e1d3fdfc3b5b8fec` (MEASURED, equals the value in `results/step3/melanin_myco/RUN.md`). `omt` = 100.6 by owner decision 2026-11-24. Step 2 tables were not changed. |
| HMMs | 18 files in `results/step3/melanin_myco/hmms_all/` (symlinks to `results/hmms/<pathway>/`). sha256 list: `hmms_all.sha256` (this folder). |
| HMM identity with the timing job | No checksum was recorded at the timing job. MEASURED: all 18 HMM files have mtime 2026-11-23 20:56:33 to 20:56:43. The timing job 32415026 started at 23:32:10 the same day. The Q8NKC1 re-run (job 32415136, 23:43) changed only the step 2 tables. So the HMM files are older than the timing job and were not rebuilt after it. This is the evidence; a byte comparison is not possible. |
| Manifest | `results/manifests/chunk_manifest.size4000.tsv.zst` (6 chunks), `results/manifests/genomes.tsv.zst`. Already in this worktree from the timing step. |

## Commands

```
cd analysis/pathway_census
export PC_ROOT=$PWD PC_PATHWAY=$PWD/pathways/melanin_dhn.yaml PC_HMM_DIR=$PWD/results/step3/melanin_myco/hmms_all \
  PC_CHUNKS=$PWD/results/manifests/chunk_manifest.size4000.tsv.zst \
  PC_MIN_SCORES=$PWD/results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv \
  PC_OUT_DIR=$PWD/results/search/melanin_myco_full PC_PROVISIONAL=1
cd $PC_OUT_DIR
for i in 0 1 2 3 4 5; do PC_CHUNK_ID=$i sbatch --export=ALL -p intel --mem=96G --time=04:00:00 --job-name=pc_mm_$i $PC_ROOT/slurm/search_chunk.sbatch; done
cd $PC_ROOT && /usr/bin/python3.12 curation_scripts/melanin_myco_step4_summary.py results/search/melanin_myco_full results/manifests results/search/melanin_myco_full/summary
```

## Jobs

First submission (partition `stajichlab`, `--mem=96G`, `--time=02:00:00`): jobs 32427757 to 32427762. All six stayed PENDING for 21 min with reason Priority. The scheduler estimate of the start was 2026-11-30 (MEASURED, `squeue`). The user held 488 CPUs in `stajichlab` (MEASURED, `squeue`), and 16-CPU jobs of the antismash run (job ids 32412xxx) were pending ahead of mine with estimated starts every 1 to 1.5 h. I did not count them. I cancelled the six jobs (`submit.txt`).
Second submission (partition `intel`): jobs 32427878 to 32427883. All six started within 20 s on nodes i01, i02, i18, i20, i20, i21. Changes from the instruction: partition `intel` instead of `stajichlab`, and `--time=04:00:00` instead of 02:00:00. Reason for the time: the speed of `intel` nodes was unknown. The same script, HMMs and table were used. UNVERIFIED: whether results are identical on a ryzen node (they should be; hmmsearch is deterministic).

| Chunk | Job | State | Node | Genomes | Proteins | Rows | Script time (s) | Job elapsed |
|---|---|---|---|---|---|---|---|---|
| 0 | 32427878 | COMPLETED | i01 | 4,000 | 41,749,741 | 29,525 | 6,270 | about 1 h 46 min |
| 1 | 32427879 | COMPLETED | i02 | 4,000 | 41,942,362 | 29,222 | 6,200 | about 1 h 44 min |
| 2 | 32427880 | COMPLETED | i18 | 4,000 | 41,990,403 | 29,132 | 6,438 | 1 h 46 min 32 s or more (last seen RUNNING at 1:46:32) |
| 3 | 32427881 | COMPLETED | i20 | 4,000 | 41,514,857 | 29,255 | 6,229 | about 1 h 44 min |
| 4 | 32427882 | COMPLETED | i20 | 4,000 | 41,930,219 | 28,755 | 6,226 | about 1 h 44 min |
| 5 | 32427883 | COMPLETED | i21 | 3,376 | 35,324,293 | 24,823 | 5,224 | about 1 h 27 min |

Script time, genomes, proteins and rows: MEASURED (`pc_chunk_*.log`, `timing.chunk*.tsv`). State: each job wrote its outputs and its log line; `sacct` was not reachable after the run (slurmdb timeout), so State, Elapsed, MaxRSS and ExitCode from `sacct` are UNVERIFIED. Elapsed is script time plus about 1 min (the timing job: 261 s elapsed, 209 s script). MaxRSS: `sstat` of chunk 0 at about 1 h 36 min showed 23,736,548 K (23.7 GB; a running value, not the final peak). The 96 GB request was enough (no OOM, no TIMEOUT).
Timing against the projection (PROPOSED 63 to 75 min, from ryzen): observed script time 5,200 to 6,440 s (87 to 107 min). That is about 1.5 times the projection. The hmmsearch time of the three PKS HMMs was 800 to 920 s each per chunk (ryzen projection: 20 x 20 to 22 s = 400 to 440 s). The `intel` nodes are slower than the ryzen node for this work, or the shared file system was busy. UNVERIFIED which.

## Checks (MEASURED, `summary.txt`)

| Check | Result |
|---|---|
| Every chunk wrote `searched.chunkNNNNN.tsv` | PASS (6 files, 4,000 x 5 + 3,376 rows) |
| Union of searched genomes | PASS: 23,376 unique. Equals the 23,376 genomes with `has_proteins` in `genomes.tsv.zst` (0 missing, 0 extra) |
| A genome searched twice | PASS: 0 duplicates |
| Status of genomes | PASS: all 23,376 `searched` |
| Proteins searched | 244,451,875. Carotenoid full search: 244,451,335 (difference +540) |
| Bad id prefixes | FAIL: 39,700 proteins in 5 genomes (all proteins of each genome). See below. |

Bad id prefixes. The 5 genomes: F8B45790 (*Blastobotrys sasensis*, 5,806 proteins), FE750E90 (*Arxiozyma pintolopesii*, 5,021), FE41F6C9 (*Fusarium oxysporum* Foc035, 17,618), FF465DA3 (*Saccharomyces cerevisiae* YJM1417, 4,391), F6D07108 (*Sugiyamaella casensis*, 6,864). The carotenoid search had 0 bad prefixes for the same genomes (5,805; 5,034; 16,929; 5,549; 6,867 proteins). MEASURED: the protein files of these genomes were rewritten on 2026-11-24 (mtime 09:46 and 10:03 for the two I checked; other files not checked), after `genomes.tsv.zst` was built. The new files carry a different locus tag: F7CB232C, F317F9F5, FB9275E9, F629B7B9, FAD1592B (MEASURED, one prefix per file). None of the five new prefixes is in `genomes.tsv.zst`. The genomes were searched in full (all proteins counted), so the search is complete. Effects: 13 candidate rows (prefixes FB9275E9 x 11, FAD1592B x 1, F7CB232C x 1) cannot be mapped to a genome row by prefix. The per-genome counts in this record use the prefix, so they include 3 unmapped prefixes. Step 5 needs the genome table rebuilt for these 5 genomes (owner decision).

## Results

Files in `summary/`: `candidates_per_step.tsv`, `tyrosinase_rows_per_genome.tsv`, `omt_low_band_genomes.tsv`, `taxon_breakdown.txt`. Also `summary.txt` (full script output). Candidate tables: `candidates.chunk0000N.PROVISIONAL.tsv.zst` (not committed).

Totals (MEASURED): 170,712 rows. 161,355 unique proteins. 19,645 genome prefixes with at least one candidate. Projection from RUN.md of step 3: 159,773 rows. Ratio 1.07.
By pathway: melanin_dhn 62,329 rows (53,024 unique proteins; 11,254 genomes); melanin_ldopa 75,909 rows (75,909; 19,073); mycosporine 32,474 rows (32,422; 12,739).

See the table in the note `docs/superpowers/notes/2026-11-24-melanin-myco-step4-results.md` for rows, unique proteins, genomes and the projection ratio per step.
