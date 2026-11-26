# Gain and loss of the placed L-DOPA melanin steps (tyrosinase, hmga)

Date: 2026-11-25 to 2026-11-26. MEASURED unless marked PROPOSED, UNVERIFIED. Calls PROVISIONAL. Branch `pathway-census-gainloss-melanin-myco`. Method, inputs, commands and versions are the same as in `../mycosporine/RUN.md` (read it for the exact command lines; only `--config`, `GL_OUT`, `GL_PRIMARY`, `NULL_OUT` and `GL_PAIRS` change). R 4.2.0; castor 1.7.10; phytools 1.2.0; ape 5.7.1; Python 3.12.15. Tree `../tree.noMicrosporidia.nwk` (identical to the carotenoid file). Partition `intel` with `--constraint=""` (the first submissions to `short` stayed pending, AssocGrpCpuLimit, and were cancelled: jobs 32464399, 32464401, 32464402, 32464404). At most 16 CPUs of mine at once (2 CPUs are allocated per job). No BUSCO filter. Null: 4 array tasks of 500 replicates per variant (N = 2,000), variants `stationary` (as specified) and `stationary_match0.1` (prevalence window 10%). Summaries: `bin/summarise_null_models.R`, `bin/summarise_gainloss_set.py`. Figures: `reports/melanin_myco_gainloss/figures/bin/fig_gainloss_{tree,counts,coevents}.R`.

## Steps
tyrosinase: `step5b/.../melanin_ldopa/tyrosinase_v2_allplace`; hmga: `.../melanin_ldopa/hmga`. Not modelled (`steps.config.all.tsv`, `calls_all_steps.tsv`): laccase_diphenol (`laccase_diphenol_v2_allplace`; 1 present of 1,650 representatives with data), asp_melanin_tyrosinase (0), aspulvinone_e_synthase (base run; 0). `hppd` is "family present", not placed, not run.

## Jobs (all COMPLETED)
| job | what | elapsed |
|---|---|---|
| 32464411 | ASR, events, permutations | 00:00:19 |
| 32464417 (4 tasks) | null, stationary | 10.8 to 11.1 min per task |
| 32464418 (4 tasks) | null, prevalence-matched | 17.8 to 18.4 min per task |
| 32464419 | Pagel tyrosinase x hmga | 25.8 min (script 1,537 s) |

## Results (from `summary_*.tsv`)
tyrosinase 1,146 present / 463 absent: ARD (delta AIC -6.6), 14 gains, 47 losses, 68 uncertain (ER: 19 gains). hmga 11 / 1,598: ARD (-4.8), 3 gains, 1 loss. No shared event in any of the 4 event types. Pagel LR 14.64, 4 df, p 0.0055.
