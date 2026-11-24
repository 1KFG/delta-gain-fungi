# Census step 3 for melanin_dhn, melanin_ldopa, mycosporine: timing run

Date: 2026-11-23. Branch `pathway-census-melanin-myco`. All outputs are PROVISIONAL (the three pathway files have `verified: false`).
Labels: MEASURED = read from a file or `sacct` output named here. PROPOSED = a choice or an estimate. UNVERIFIED = not checked.
This is a timing run on 200 genomes. It is not the full search. The full search was not run.

## Inputs

| Item | Value |
|---|---|
| HMMs | 18 step HMMs (9 DHN, 6 L-DOPA, 3 mycosporine), symlinked into `results/step3/melanin_myco/hmms_all/` from `results/hmms/<pathway>/` |
| Cutoff table | `results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv`: the three `results/step2/<pathway>/min_scores.PROVISIONAL.tsv` tables joined (18 rows). `omt` = 100.6 (owner override). sha256 `9b9391230fb5e4d28d2b7aa9ca49f7f4d81388ec8e506662e1d3fdfc3b5b8fec` |
| Chunk manifest | `results/manifests/chunk_manifest.size200.tsv.zst`, chunk 0. Copied (read only) from the carotenoid worktree. Same chunk 0 as the carotenoid timing job 32375142 |
| Full-search manifest (not used yet) | `results/manifests/chunk_manifest.size4000.tsv.zst` (6 chunks), `results/manifests/genomes.tsv.zst`, `results/manifests/chunk_excluded.tsv` |
| Tool change | None. `bin/search_chunk.py` takes one `--hmm-dir` and one `--min-scores` table. A pathway is not needed for the search. `--pathway` is used only for the `verified` gate. A combined table and one HMM directory with all 18 HMMs gives one search for the three pathways. Step names are unique across the three pathways (MEASURED, 18 distinct names). |

## Commands

```
cd analysis/pathway_census
export PC_ROOT=$PWD PC_PATHWAY=$PWD/pathways/melanin_dhn.yaml PC_HMM_DIR=$PWD/results/step3/melanin_myco/hmms_all \
  PC_CHUNKS=$PWD/results/manifests/chunk_manifest.size200.tsv.zst PC_CHUNK_ID=0 \
  PC_MIN_SCORES=$PWD/results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv \
  PC_OUT_DIR=$PWD/results/step3/melanin_myco/timing_size200 PC_PROVISIONAL=1
cd results/step3/melanin_myco/timing_size200
sbatch -p stajichlab --constraint="ryzen|cascade|broadwell" $PC_ROOT/slurm/search_chunk.sbatch
```

The `--pathway` file is melanin_dhn only to pass the gate. The other two pathway files also have `verified: false` (MEASURED, `grep` of the YAML).

## Jobs

| Job | Partition | State | Note |
|---|---|---|---|
| 32412179 | short | CANCELLED while PENDING | Pending about 3 h with reason AssocGrpCpuLimit. The user limit on `short` is 384 CPUs. The running antismash Nextflow jobs of the same user held all 384 for the whole wait (MEASURED, `squeue`). I cancelled it. |
| 32415026 | stajichlab | COMPLETED | Same script, same inputs, same 16 CPUs. Ran on node r11 (ryzen). The queue name is the only change. This is a deviation from "short queue". The partition `stajichlab` has ryzen (milan) nodes. The short partition also has ryzen/milan nodes. The timing may differ from a `short` node. UNVERIFIED how much. |

## Timing (MEASURED, job 32415026, `timing.chunk00000.tsv`, `sacct`)

| Item | This run (18 HMMs) | Carotenoid job 32375142 (14 HMMs) |
|---|---|---|
| Job elapsed | 4 min 21 s (261 s) | 1 min 31 s |
| Script time (`total`) | 209.3 s | 78.9 s |
| Stage and count proteins | 66.2 s | 20.8 s |
| hmmsearch, sum of 18 (or 14) steps | 143.1 s | 58.1 s |
| Peak memory (MaxRSS) | 3,108,144 K (3.1 GB) | 1.48 GB |
| Genomes searched | 200 of 200, all `searched` | 200 |
| Proteins | 2,099,680 | 2,099,659 |
| Bad id prefixes | 0 (all 200 genomes) | not recorded |
| Candidate rows | 1,367 (1,299 unique proteins) | 1,827 |

The protein count differs by 21 from the carotenoid job. One genome (FB53F6ED, *Batrachochytrium dendrobatidis* RTP6) has 8,532 proteins now and had 8,511 then. The protein file changed in between. UNVERIFIED why.
The staging time is 3.2 times the carotenoid value. The node and the file cache differ. UNVERIFIED which one. The stage time is not a function of the HMM count.

hmmsearch seconds per step (MEASURED): ayg1 5.3; heptaketide_pks 22.1; hexaketide_pks 20.3; laccase_abr2 5.0; mco_abr1 6.1; scytalone_dehydratase 4.1; t3hn_reductase 5.9; t4hn_reductase 6.3; thn_synthase 21.5; asp_melanin_tyrosinase 4.6; aspulvinone_e_synthase 8.8; hmga 4.6; hppd 4.6; laccase_diphenol 5.1; tyrosinase 4.8; atp_grasp_ligase 4.8; ddgs 4.6; omt 4.6.
The three PKS HMMs take 63.9 s (45% of the hmmsearch time). The other 15 HMMs take 4.1 to 8.8 s each.

## Candidates per step in the timing chunk (MEASURED, `timing_size200/candidates_per_step.tsv`)

Rows are (protein, step) pairs at or above `min_score`. `x116.88` = 23,376 / 200. This scaling is PROPOSED, an estimate from one seeded 200-genome chunk.
Check on the scaling (MEASURED from the carotenoid note): the carotenoid timing chunk gave 1,827 rows. 1,827 x 116.88 = 213,540. The full carotenoid run gave 238,898 rows (ratio 1.12). The two runs used different cutoffs for two steps, so the ratio is an approximate guide only.

| Step | min_score | Rows in chunk | Score range seen | Projected rows, 23,376 genomes (PROPOSED) |
|---|---|---|---|---|
| tyrosinase | 63.6 | 439 | 66.3 to 751.6 | 51,310 |
| hppd | 164.3 | 154 | 171.5 to 974.9 | 18,000 |
| omt | 100.6 | 125 | 100.8 to 160.6 | 14,610 |
| ayg1 | 374.2 | 87 | 395.6 to 807.6 | 10,169 |
| atp_grasp_ligase | 105.0 | 86 | 118.6 to 382.6 | 10,052 |
| ddgs | 463.5 | 81 | 477.9 to 613.7 | 9,467 |
| t4hn_reductase | 356.0 | 71 | 356.2 to 523.6 | 8,298 |
| scytalone_dehydratase | 240.5 | 64 | 241.4 to 335.4 | 7,480 |
| hexaketide_pks | 1839.0 | 63 | 1842.6 to 5166.2 | 7,363 |
| laccase_abr2 | 464.3 | 63 | 464.5 to 1455.7 | 7,363 |
| thn_synthase | 1908.6 | 53 | 1919.8 to 4940.6 | 6,195 |
| t3hn_reductase | 550.4 | 21 | 560.8 to 614.7 | 2,454 |
| heptaketide_pks | 3656.7 | 15 | 3734.6 to 4850.3 | 1,753 |
| hmga | 983.2 | 15 | 993.4 to 1078.5 | 1,753 |
| asp_melanin_tyrosinase | 320.5 | 13 | 325.3 to 415.6 | 1,519 |
| mco_abr1 | 638.2 | 13 | 650.1 to 1069.3 | 1,519 |
| laccase_diphenol | 1127.7 | 4 | 1272.3 to 1437.6 | 468 |
| aspulvinone_e_synthase | 1205.5 | 0 | none | 0 |
| Total | | 1,367 | | 159,773 |

Overlap between steps (MEASURED): 1,299 unique proteins in 1,367 rows. The three PKS steps have 131 rows on 67 proteins (62 proteins hit two steps, 1 hits three, 4 hit one). The two SDR steps have 92 rows on 88 proteins. The two laccase steps have 67 rows on 67 proteins.

Steps with a large count (reading of counts only, not findings):
- `tyrosinase`: 439 rows in 200 genomes (2.2 per genome). 33 rows are below 2 x the cutoff (127.2). The cutoff 63.6 is set by one reference (Q00234). Half of the rows score above 318. Needs step 5 placement.
- `hppd`: 154 rows. The floor 164.3 is the score of the decoy 4-hydroxymandelate synthase. 133 genomes have a hit. 26 rows are below 2 x the floor.
- `omt`: 125 rows on 125 proteins in 106 genomes (median score 111.4, mean 116.7). All rows are at 100.8 to 160.6, near the floor. This does not match the step 2 side check (section 9.5 of the step 2 note), where almost no reviewed fungal O-methyltransferase scored above 100.6. The BFD hits are not reviewed entries. UNVERIFIED what they are.
- PKS steps (`hexaketide_pks`, `thn_synthase`, `heptaketide_pks`): counts are low (63, 53, 15), but 62 of 67 proteins hit two HMMs. The HMMs do not separate the PKS steps by score (known from step 2). The 3,656.7 to 1,908.6 cutoffs are high, so these counts are small.
- `aspulvinone_e_synthase`: 0 rows in 200 genomes. The cutoff 1205.5 is 54% of the reference score.

## Projection for the full search (PROPOSED; UNVERIFIED until the first full chunk reports)

Basis: the carotenoid chunk-of-200 script time of 78.9 s gave 26.3 min predicted for 4,000 genomes. The observed times were 23.8 to 28.0 min for six jobs (ratio 0.90 to 1.06; carotenoid step 4 note). Scaling the same way (20 x):

| Item | Value |
|---|---|
| Script time, 4,000-genome chunk | 209.3 s x 20 = 69.8 min. With the carotenoid ratio 0.90 to 1.06: 63 to 74 min. |
| Job time | script time plus about 1 min start-up: about 64 to 75 min |
| Peak memory | 3.1 GB x 15.2 (the carotenoid growth from 200 to 4,000 genomes: 1.48 to 22.5 GB) = about 47 GB. At 20 x: about 62 GB. The default `--mem=32G` in `slurm/search_chunk.sbatch` is too small. |
| Candidate rows, 6 chunks | 159,773 (1,367 x 116.88) |

Recommendation (PROPOSED): keep 6 chunks of 4,000 genomes (`chunk_manifest.size4000.tsv.zst`), 16 CPUs, `--mem=96G`, `--time=02:00:00`. The projected job time (about 64 to 75 min) is inside the 1 to 1.5 h target and inside the 2 h limit of the short queue with a margin of 1.6 x. Six jobs use about 6 x 70 min = 7 h of job time. If all six run at once, the wall time is about 70 min. A chunk of 5,000 genomes (5 jobs) would project to about 87 min; the margin to 2 h falls to 1.4 x and the memory estimate rises to about 59 to 78 GB. Not recommended.
Memory reason (UNVERIFIED): `parse_domtbl` reads every hit at E <= 10 into Python lists before the score filter. A tool change to filter during the parse would cut memory. It is not made here.
Queue risk: job 32412179 waited about 3 h on the `short` queue while the antismash run held the 384-CPU user limit. Six 16-CPU jobs need 96 CPUs. Submit when that run is done, or use the `stajichlab` partition as in job 32415026.
