# Any DHN PKS with the six non-PKS DHN steps

Date: 2026-11-25 to 2026-11-26. MEASURED unless marked PROPOSED, UNVERIFIED. Calls PROVISIONAL. Branch `pathway-census-gainloss-melanin-myco`. Method, inputs, commands and versions are the same as in `../mycosporine/RUN.md` (read it for the exact command lines; only `--config`, `GL_OUT`, `GL_PRIMARY`, `NULL_OUT` and `GL_PAIRS` change). R 4.2.0; castor 1.7.10; phytools 1.2.0; ape 5.7.1; Python 3.12.15. Tree `../tree.noMicrosporidia.nwk` (identical to the carotenoid file). Partition `intel` with `--constraint=""` (the first submissions to `short` stayed pending, AssocGrpCpuLimit, and were cancelled: jobs 32464399, 32464401, 32464402, 32464404). At most 16 CPUs of mine at once (2 CPUs are allocated per job). No BUSCO filter. Null: 4 array tasks of 500 replicates per variant (N = 2,000), variants `stationary` (as specified) and `stationary_match0.1` (prevalence window 10%). Summaries: `bin/summarise_null_models.R`, `bin/summarise_gainloss_set.py`. Figures: `reports/melanin_myco_gainloss/figures/bin/fig_gainloss_{tree,counts,coevents}.R`.

`calls.tsv` = the six non-PKS DHN steps plus `any_dhn_pks` (`bin/add_union_step.py --members thn_synthase hexaketide_pks heptaketide_pks --drop-members`), then Microsporidia dropped. Run separately from `../melanin_dhn/` so that no pair holds a union and one of its members. Present 402, absent 1,207 tips.

## Jobs (all COMPLETED except as noted)
| job | what | elapsed |
|---|---|---|
| 32464410 | ASR, events, permutations | 00:00:59 |
| 32465868 | null, stationary: task 1 completed (38.9 min); tasks 2 to 4 cancelled to save compute | partial, not used (`../melanin_dhn_anypks_null/stationary_partial_task1_only/`) |
| 32468344 (4 tasks) | null, prevalence-matched, N = 2,000 | 47.0 to 61.5 min for tasks 1 to 3 (task 4 ran 50.1 min or more when last seen) |
| 32470277 (6 tasks) | Pagel, `any_dhn_pks` against the 6 non-PKS steps | 13.3 to 38.3 min per task |

## Deviations
- Only the matched null has N = 2,000 here. The as-specified null was stopped after 500 replicates.
- BH q is over the 84 tests of this 7-step run.
