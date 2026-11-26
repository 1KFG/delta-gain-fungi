# Gain and loss of the placed mycosporine steps (ddgs, omt)

Date: 2026-11-25. MEASURED unless marked PROPOSED, UNVERIFIED or REPORTED. Calls PROVISIONAL (pathway `mycosporine` verified: false). Branch `pathway-census-gainloss-melanin-myco`, worktree `DeltaGain_Fungi_wt_gainloss_mm`. Same method as `results/gainloss/carotenoid_noMicro_blasto_5steps/` (read only). No BUSCO filter.

## Steps (`steps.config.tsv`)
| step | source (`assigned.PROVISIONAL.tsv`, label `assigned` = present) | present / absent / NA among 1,654 genus representatives (`calls_all_steps.summary.tsv`) |
|---|---|---|
| ddgs | `DeltaGain_Fungi_wt_step5b/.../results/step5/mycosporine/ddgs_allplace` (all-placements row) | 667 / 983 / 4 |
| omt | `.../mycosporine/omt_v2` (single-reference placement; no `_allplace` run exists) | 638 / 1,012 / 4 |

Steps of the pathway that are not placed (not analysed): atp_grasp_ligase (search-hit step), mg_glucosyltransferase, glutamicol_formation.

## Commands
```
# env: module unload R; module load R/4.2.0 ; Python /usr/bin/python3.12 ; cd analysis/pathway_census ; G=results/gainloss
REPS=BFD_dashboard/atlas_data/phylo/genus_representatives.tsv ; TAX=Fungi_BFD_runs/tables/species.parquet
SD=DeltaGain_Fungi_wt_melanin/analysis/pathway_census/results/search/melanin_myco_full   # searched.chunk0000{0..5}.tsv
TREE=BFD_dashboard/atlas_data/phylo/results/protein/buildtree/fungi_odb10/fasttree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fasttree.support.treefile
/usr/bin/python3.12 bin/make_calls_table.py --reps $REPS --taxonomy $TAX --search-dir $SD --config $G/mycosporine/steps.config.tsv --root $PWD --prefix-aliases prefix_aliases.tsv --out $G/mycosporine/calls_pre_drop.tsv
/usr/bin/python3.12 bin/drop_phylum_calls.py --calls $G/mycosporine/calls_pre_drop.tsv --phylum Microsporidia --out $G/mycosporine/calls.tsv --dropped-list $G/mycosporine/dropped_tips.txt
Rscript bin/prune_tree_tips.R $TREE $G/mycosporine/dropped_tips.txt $G/tree.noMicrosporidia.nwk     # 1,643 -> 1,613 tips
export PC_ROOT=$PWD GL_TREE=$PWD/$G/tree.noMicrosporidia.nwk GL_ROOTS=Blastocladiomycota,Chytridiomycota GL_NPERM=10000 GL_OUT=$PWD/$G/mycosporine
sbatch --export=ALL -p intel --constraint="" --time=01:00:00 slurm/gainloss_run_roots.sbatch                         # ASR, events, permutation test
export GL_PRIMARY=$PWD/$G/mycosporine NULL_OUT=$PWD/$G/mycosporine_null NULL_NREP=500 NULL_TOL=0                     # stationary root, as specified
sbatch --export=ALL -p intel --constraint="" --array=1-4%1 slurm/gainloss_null_models.sbatch                         # 4 x 500 = 2,000 replicates
NULL_TOL=0.1 sbatch --export=ALL -p intel --constraint="" --array=1-4%1 slurm/gainloss_null_models.sbatch            # prevalence-matched (window 10%)
export GL_OUT=$PWD/$G/mycosporine GL_PAIRS=$GL_OUT/pagel_pairs.tsv
sbatch --export=ALL -p intel --constraint="" --time=04:00:00 --array=1-1 slurm/gainloss_pagel.sbatch                 # Pagel, one pair
Rscript bin/summarise_null_models.R $G/mycosporine $G/mycosporine_null                                               # null_summary.tsv
/usr/bin/python3.12 bin/summarise_gainloss_set.py --set-dir $G/mycosporine --null-dir $G/mycosporine_null           # summary_*.tsv
Rscript reports/melanin_myco_gainloss/figures/bin/fig_gainloss_{tree,counts,coevents}.R ...                          # figures (see note)
```
Absolute paths of the inputs are in the note `docs/superpowers/notes/2026-11-25-gainloss-melanin-myco.md`. The tree is the same FastTree support tree as the carotenoid run: `cmp` of `tree.noMicrosporidia.nwk` with `carotenoid_noMicro_blasto_5steps/tree.noMicrosporidia.nwk` shows identical files; `dropped_tips.txt` is identical too (30 Microsporidia tips).

## Versions
R 4.2.0; castor 1.7.10; phytools 1.2.0; phangorn 2.10.0; ape 5.7.1; ggtree 3.4.0; ggplot2 3.4.2; Python 3.12.15.

## Jobs (partition `intel`, 2 CPUs allocated per job; at most 16 CPUs of mine at once)
| job | what | elapsed |
|---|---|---|
| 32464404 | first submission, partition short | cancelled while pending (AssocGrpCpuLimit), never ran |
| 32464412 | gainloss_run_roots (ASR, events, 10,000 permutations) | 00:00:20, COMPLETED |
| 32464420 (4 tasks) | null, stationary root | job elapsed 10.3 to 10.8 min per task of 500 replicates (script time 1.22 to 1.26 s per replicate) |
| 32464421 (4 tasks) | null, prevalence-matched | job elapsed 18.2 to 18.9 min per task (script time 2.16 to 2.24 s per replicate) |
| 32464422 | Pagel ddgs x omt | 803 s (script time), COMPLETED |

Each unit of work is small (minutes), well under the 1 to 1.5 h guideline. The null was split in 4 tasks of 500 replicates (carotenoid: 8 of 250) to reduce the number of jobs.

## Outputs
`calls.tsv`, `calls.summary.tsv`, `calls_all_steps.tsv` (all placed steps of the pathway, incl. unmodelled), `gainloss_counts.tsv`, `co_events.tsv`, `threshold_sensitivity.tsv`, `edge_events_primary.tsv`, `primary.rds`, `tree.rooted.primary.nwk`, `pagel/`, `summary_steps.tsv`, `summary_pairs.tsv`, `summary_pagel.tsv`, figures `gl_mycosporine_*.png|pdf`. Null: `../mycosporine_null/` (`null_summary.tsv`, `null_step_diagnostics.tsv`, `observed_co_counts.tsv`, `sim_*.tsv`, `sim_params.tsv`).

## Results in one place (all from `summary_steps.tsv`, `summary_pairs.tsv`, `null_summary.tsv`)
- 1,609 tips with data (1,613 tips in the pruned tree; 4 NA). ddgs 664 present / 945 absent; omt 635 / 974.
- ARD better than ER (delta AIC -151.1 and -247.9). ddgs: 8 gains, 139 losses, 111 uncertain; omt: 3 gains, 152 losses, 93 uncertain.
- Matched null, co-loss all edges: 79 observed, null mean 16.68 (7 to 29), p 0.0005, q 0.002. Internal edges: 6 observed, null mean 2.19, p 0.031, q 0.124. Co-gain: 1 observed, null 0.06, p 0.059, q 0.118.
- Pagel: LR 312.0, 4 df, p 2.8e-66.

## Failures and deviations
- The first submission to partition short stayed pending (AssocGrpCpuLimit) and was cancelled; all jobs ran on `intel` with `--constraint=""` (the scripts' constraint names ryzen, cascade and broadwell only).
- BH q in `null_summary.tsv` is over the 4 event types of one pair (4 tests per scope), not over 40 as in the carotenoid run, because there is only one pair.
- Unmodelled steps: none in this set (both placed steps are modelled).
