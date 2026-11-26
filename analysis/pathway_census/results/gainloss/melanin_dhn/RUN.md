# Gain and loss of the placed DHN melanin steps (9 steps)

Date: 2026-11-25 to 2026-11-26. MEASURED unless marked PROPOSED, UNVERIFIED. Calls PROVISIONAL. Branch `pathway-census-gainloss-melanin-myco`. Method, inputs, commands and versions are the same as in `../mycosporine/RUN.md` (read it for the exact command lines; only `--config`, `GL_OUT`, `GL_PRIMARY`, `NULL_OUT` and `GL_PAIRS` change). R 4.2.0; castor 1.7.10; phytools 1.2.0; ape 5.7.1; Python 3.12.15. Tree `../tree.noMicrosporidia.nwk` (identical to the carotenoid file). Partition `intel` with `--constraint=""` (the first submissions to `short` stayed pending, AssocGrpCpuLimit, and were cancelled: jobs 32464399, 32464401, 32464402, 32464404). At most 16 CPUs of mine at once (2 CPUs are allocated per job). No BUSCO filter. Null: 4 array tasks of 500 replicates per variant (N = 2,000), variants `stationary` (as specified) and `stationary_match0.1` (prevalence window 10%). Summaries: `bin/summarise_null_models.R`, `bin/summarise_gainloss_set.py`. Figures: `reports/melanin_myco_gainloss/figures/bin/fig_gainloss_{tree,counts,coevents}.R`.

## Steps (`steps.config.tsv`; absolute sources in `steps.config.all.tsv` is the same list)
thn_synthase, hexaketide_pks, heptaketide_pks: `DeltaGain_Fungi_wt_melanin/.../step5/melanin_dhn/pks_joint_id0.8/assigned.<step>.PROVISIONAL.tsv` (per-step decoy rule). ayg1, t4hn_reductase, t3hn_reductase, laccase_abr2, mco_abr1: `DeltaGain_Fungi_wt_step5b/.../step5/melanin_dhn/<step>/`. scytalone_dehydratase: `.../scytalone_dehydratase_v2_allplace/`. Tips with data 1,609; present counts in `summary_steps.tsv`.

## Jobs (all COMPLETED)
| job | what | elapsed |
|---|---|---|
| 32464408 | ASR, events, 10,000 permutations | 00:01:18 |
| 32464428 (4 tasks) | null, stationary | 47.6 to 49.4 min per task (4.64 s per replicate in task 1) |
| 32464429 (4 tasks) | null, prevalence-matched | 81.4 to 88.8 min per task (9.74 s per replicate in task 1; 6,053 draws per accepted hexaketide_pks replicate) |
| 32464510 (36 tasks) | Pagel, one pair per task | 12.1 to 48.8 min per task; script time 723 to 2,927 s |

Each unit is well under the 1 to 1.5 h guideline for the cheap jobs; the matched null tasks are about 1.4 h.

## Results in one place (from `summary_steps.tsv`, `summary_pairs.tsv`, `summary_pagel.tsv`)
ARD better than ER for all 9 steps. Matched null: 5 of 72 co-gain and co-loss tests with q below 0.05 (co-loss of t4hn_reductase with thn_synthase 17 vs 3.10, ayg1 20 vs 6.01, laccase_abr2 9 vs 1.26, mco_abr1 17 vs 2.72; co-gain with scytalone_dehydratase 4 vs 0.11). BH q is over 144 tests (floor 0.036 at p = 0.0005, N = 2,000). Pagel: 36 of 36 run; 3 with a negative likelihood ratio (non-converged fits).

## Failures and deviations
- Partition `short` blocked (see above).
- Null array task 2 of job 32464428 and task 4 of job 32464429 were held by the scheduler for a while (state PENDING, "None assigned"); they completed.
- `heptaketide_pks`: 9 present tips; the prevalence window of 10% allows only 9.
