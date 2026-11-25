# Step 5 for the other melanin and mycosporine steps: run record (2026-11-24)

Status: PROVISIONAL. Label rule: MEASURED = read from a file named here or from `sacct`; PROPOSED = a choice made here; UNVERIFIED = not checked.
Branch `pathway-census-step5-other` (from `pathway-census-melanin-myco`). Worktree `DeltaGain_Fungi_wt_step5b`. No push, no merge. Nothing under `refs/`, `pathways/` or an existing tool was changed; the README got one appended section. Python `/usr/bin/python3.12`. Note: `docs/superpowers/notes/2026-11-24-melanin-myco-step5-other-steps.md`.

## Inputs (read-only copies or links in this worktree)
- `results/hmms`, `results/step2` (links to the PKS worktree `DeltaGain_Fungi_wt_melanin`), `results/genomes.tsv.zst` and `results/manifests/` (copies), `results/search/melanin_myco_full/candidates.chunk0000*.PROVISIONAL.tsv.zst` (copies; `zstd` refuses a symbolic link as input), `results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv`.
- The genome table is the one of the step 4 note (5 re-annotated genomes not corrected).

## Commands (from `analysis/pathway_census`; copies of the scripts are in `curation_scripts/melanin_myco_step5/`)
1. Family fetch: `curation_scripts/melanin_myco_step5/fetch_family.py --out results/step5/other_inputs/family PF13561 PF06500 PF02982 PF00394 PF00264 PF04209 PF00903 PF00501 PF01761 PF01596 PF02655` and, for atp_grasp_ligase, `--max-other 300 PF02222 PF02786 PF07478 PF08443 PF13535`. Log: `other_inputs/family/fetch_log.tsv` (32 sets, all `same_accessions` and `lengths_equal` true). Reference Pfam ids and EC: `other_inputs/ref_pfam.tsv` (fetched 2026-11-24).
2. Screens: `screen_family.sh` (`bin/screen_outgroups.py`, `bin/filter_decoys_for_tree.py --min-score 10 --coverage-mode best|union`).
3. Leaves: `make_leaves.sh` (`bin/select_family_leaves.py --mode best --max-score <floor> --exclude-ec ... --exclude-gene ... --n-named 10 --n-outgroup 6`), `assemble_inputs.sh` and `assemble_sdr_joint.sh` (`bin/cut_leaves_to_domain.py`, pad 20). Outputs: `other_inputs/leaves/`, `other_inputs/tree_inputs/<step>/`.
4. Cluster counts: `step5_place.py extract` of each step to `$SCRATCH`, then `mmseqs easy-cluster --min-seq-id <id> -c 0.8 --cov-mode 0` (table below).
5. Placement: `submit_step.sh <pathway> <step> <identity> intel` = `sbatch --parsable --export=ALL -p intel -c 15 -D <run folder> slurm/step5_prepare_ident.sbatch`, then `slurm/step5_multi_ref_clades.sbatch` (afterok), then `slurm/step5_place.sbatch` (afterok). Environment: `PC_ROOT=$PWD PC_PATHWAY PC_STEP PC_CAND_DIR=results/search/melanin_myco_full PC_OUT PC_DECOYS_FA PC_MEMBERS_FA PC_OUTGROUP_FA PC_OUTGROUP_ROLE=tree_outgroup PC_RAXML_THREADS=4 PC_RAXML_STARTS='pars{2}' PC_MIN_SUPPORT=0.5 PC_MIN_SEQ_ID=<id>`. `PC_CAND_DIR` is read by both parts and was exported for both. Lanes: 2 chains at a time (15 CPUs each, `-c 15`), by `--dependency=afterany` on the last place job of the lane, to stay at 30 CPUs.
6. SDR joint: `submit_sdr_joint.sh 0.8 intel`: prepare on `results/search/sdr_joint/` (rows of both steps relabelled `t4hn_reductase`; 11,667 proteins), then `derive_sdr_step.sh <step>` for each step (links to the tree files, domains of that step, `step5_multi_ref_clades.py` with the references of that step and the references of the other step as named members), then `step5_place.sbatch` per step with `PC_CAND_DIR=results/search/melanin_myco_full`.
7. Reports: `bin/step5_run_report.py --run <folder> --genomes results/genomes.tsv.zst`; `curation_scripts/melanin_myco_step5/collate.py` and `collate_tables.py`; `bin/target_overlap.py`.

## Floors (PROPOSED = owner floors; no placement floor; rows dropped 0)
Step, floor, rows, median score, share of rows below 1.1 x floor: `other_inputs/candidate_score_distribution.tsv` (force-added copy `candidate_score_distribution.tsv`). Floors: aspulvinone_e_synthase 1205.5, hmga 983.2, laccase_diphenol 1127.7, mco_abr1 638.2, laccase_abr2 464.3, asp_melanin_tyrosinase 320.5, ayg1 374.2, scytalone_dehydratase 240.5, t4hn_reductase 356.0, t3hn_reductase 550.4, ddgs 463.5, atp_grasp_ligase 105.0, omt 100.6, tyrosinase 63.6, hppd 164.3.

## Clusters at each identity (MEASURED; mmseqs2 17-b804f, `-c 0.8 --cov-mode 0`; domains extracted)
| step | domains | 0.5 | 0.7 | 0.8 | 0.9 | 0.95 | used |
|---|---:|---:|---:|---:|---:|---:|---|
| aspulvinone_e_synthase | 168 | 2 | 14 | 21 | 30 | 38 | 0.95 |
| laccase_diphenol | 402 | 1 | 1 | 2 | 5 | 5 | 0.95 |
| mco_abr1 | 1,918 | 2 | 35 | 104 | 214 | 337 | 0.8 |
| laccase_abr2 | 9,631 | 32 | 260 | 486 | 898 | 1,302 | 0.7 |
| asp_melanin_tyrosinase | 2,471 | 7 | 57 | 105 | 190 | 249 | 0.8 |
| ayg1 | 11,835 | 30 | 264 | 682 | 1,509 | 2,225 | 0.7 |
| ddgs | 8,865 | 8 | 43 | 131 | 576 | 1,027 | 0.8 |
| hmga | 2,313 | 1 | 1 | 1 | 13 | 73 | 0.95 |
| scytalone_dehydratase | 8,604 | 4 | 48 | 207 | 740 | 1,242 | 0.8 |
| omt | 14,439 | 66 | 420 | 839 | 1,537 | 2,143 | 0.7 (cap 300) |
| hppd | 20,098 | 78 | 662 | 1,395 | 2,445 | not run | 0.7 (cap 300) |
| tyrosinase | 50,447 | 2,460 | 5,820 | 7,928 | 10,938 | 13,352 | 0.5 (cap 300) |
| SDR union (t4hn + t3hn) | 11,667 | 3 | 86 | 250 | 773 | not run | 0.8 |
Rule (PROPOSED): the lowest identity with at least 100 clusters, 0.95 if none; cap 300 (stratified by phylum, seed 1). Round 1 ran every step at 0.5 first; 8 of those runs finished and are kept as `<step>_id50` (diagnostic).

## Jobs (MEASURED, `sacct`; partition intel; 16 CPUs allocated per prepare/place job)
| step | identity | prepare job (elapsed) | clades job | place job (elapsed) | State |
|---|---|---|---|---|---|
| ayg1 | 0.7 | 32435486 (00:15:54, COMPLETED) | 32435487 (COMPLETED) | 32435488 (00:05:58, COMPLETED) | |
| laccase_abr2 | 0.7 | 32435489 (00:15:51, COMPLETED) | 32435490 (COMPLETED) | 32435491 (00:04:41, COMPLETED) | |
| ddgs | 0.8 | 32435492 (00:05:20, COMPLETED) | 32435493 (COMPLETED) | 32435494 (00:03:10, COMPLETED) | |
| omt | 0.7 | 32435495 (00:13:22, FAILED) | 32435496 (CANCELLED by 1181) | 32435497 (00:00:00, CANCELLED by 1181) | |
| aspulvinone_e_synthase | 0.95 | 32435498 (00:01:10, COMPLETED) | 32435499 (COMPLETED) | 32435500 (00:00:29, COMPLETED) | |
| mco_abr1 | 0.8 | 32435501 (00:04:18, COMPLETED) | 32435502 (COMPLETED) | 32435503 (00:01:16, COMPLETED) | |
| laccase_diphenol | 0.95 | 32435504 (00:00:19, COMPLETED) | 32435505 (COMPLETED) | 32435506 (00:00:35, COMPLETED) | |
| asp_melanin_tyrosinase | 0.8 | 32435507 (00:01:53, COMPLETED) | 32435508 (COMPLETED) | 32435509 (00:00:34, COMPLETED) | |
| scytalone_dehydratase | 0.8 | 32435510 (00:05:07, COMPLETED) | 32435511 (FAILED) | 32435512 (00:00:00, CANCELLED by 1181) | |
| tyrosinase | 0.5 | 32435513 (00:48:48, COMPLETED) | 32435514 (COMPLETED) | 32435515 (01:17:18, COMPLETED) | |
| hmga | 0.95 | 32435516 (00:01:03, COMPLETED) | 32435517 (COMPLETED) | 32435518 (00:00:38, COMPLETED) | |
| hppd | 0.7 | 32435520 (00:17:50, COMPLETED) | 32435521 (COMPLETED) | 32435522 (00:09:25, COMPLETED) | |
| sdr_joint (t4hn_reductase, t3hn_reductase) | 0.8 | 32435712 (00:09:03) | derive 32435713, 32435715 | place 32435714 (00:02:00), 32435716 (00:00:40) | COMPLETED, COMPLETED, COMPLETED |

Round 1 (identity 0.5, `intel`): `other_steps_jobs.round1.txt`. Executed: aspulvinone_e_synthase 32434872-874, mco_abr1 32434875-877, ayg1 32434878-880, ddgs 32434881-883, laccase_diphenol 32434890-892, laccase_abr2 32434893-895, asp_melanin_tyrosinase 32434896-898, scytalone_dehydratase 32434899 (prepare; the clades job 32434900 failed, exit 3). Cancelled before they started: omt, hppd, tyrosinase, hmga, SDR (decoy-only leaf set had to be rebuilt after the gene-name exclusion and the identity choice).
Round 0 (`short`): `other_steps_jobs.short_cancelled.txt`: 38 jobs, all in `AssocGrpCpuLimit` (the user limit of 384 CPUs on `short`), cancelled by job id and resubmitted on `intel`.
Parallel hmmalign check for tyrosinase: job 32438154 (`slurm/step5_place_parallel.sbatch`, 10 min 7 s, COMPLETED), folder `tyrosinase_par`.
Timings (MEASURED): prepare 19 s to 48 min 48 s (tyrosinase), place 29 s to 1 h 17 min 18 s (tyrosinase). All within the 2 h limit. Tools: mafft 7.505, clipkit 1.3.0, RAxML-NG 2.0.4, EPA-ng 0.3.8, mmseqs2 17-b804f, HMMER 3.4.

## Failures and problems
- omt: prepare exit 3 (job 32435495): references not a clade free of decoys and outgroups. Not placed.
- scytalone_dehydratase: prepare ok (job 32435510), `step5_multi_ref_clades.py` exit 3 (job 32435511). Not placed. Also at 0.5 with the first leaf set.
- atp_grasp_ligase and the first leaf screens: no alignable leaf; no job.
- 12 candidate rows not extracted (ayg1 1, ddgs 1, hppd 5, tyrosinase 5): `protein_file_missing` (prefixes FB9275E9, F7CB232C, FAD1592B) and `bad_coordinates` (F75533D7, FADCCC61). `candidates_step.tsv` rows equal `placements.tsv` rows plus these failures for those 4 steps and are equal for the other steps (`other_steps_summary.tsv`, `candidates_equal_placements`).
- A dependency `afterany` on a job that stays `DependencyNeverSatisfied` does not start; cancelled the stuck job and removed the dependency (`scontrol update JobId=<id> Dependency=`) for 32435498 and 32435513.
- A first `scancel` pattern also matched jobs of another agent; it was blocked by the tool before it ran. After that all cancellations used job ids from the job files only.
- `module load ... | tail` runs in a subshell and loads nothing (as in CLAUDE.md); the loads were repeated without the pipe.

## Tests (MEASURED)
`/usr/bin/python3.12 -m pytest analysis/pathway_census/tests -q` (modules samtools, zstd, hmmer, mafft loaded): 924 passed, 1 failed, 1 skipped (926 collected). The failure is `test_step5_place.py::test_extract_end_to_end_with_failures` (`bad_coordinates ... length 0` instead of `id_not_in_protein_file`); it fails the same way in the carotenoid single-reference record (888 passed, 1 failed), and `bin/step5_place.py` is not changed (git diff empty). New tests: 23 in 4 files (`test_select_family_leaves.py` 10, `test_cut_leaves_to_domain.py` 4, `test_multi_ref_clades.py` 5, `test_step5_run_report.py` 4; 23 passed when run alone).

## Files
`other_steps_summary.tsv` (one row per run, 23 rows), `other_steps_assigned_by_class.tsv`, `other_steps_nearest_summary.tsv`, `other_steps_leaves.tsv`, `other_steps_score_by_call.tsv`, `tyrosinase_laccase_target_overlap.tsv`, `sdr_target_overlap.tsv`, `leaf_summary.tsv`, `decoy_alignability.tsv`, `leaves_selected.tsv`, `candidate_score_distribution.tsv` (copies of `other_inputs/` files, force-added). Per-run folders: `results/step5/<pathway>/<step>/` with `SUMMARY.tsv`, `report/`, `clade_report.json`, `placements_detail.tsv`.

## Follow-up runs v2 (2026-11-24 to 2026-11-25): owner decisions 1, 2, 4, 5
Label rule as above. Folders: `results/step5/<pathway>/<step>_v2/` (v1 untouched). New rows in `other_steps_summary.tsv` (rows 24 to 29, `*_v2`); the same rows alone are in `other_steps_summary.v2.tsv`. Job list: `other_steps_jobs.v2.txt`. Submit script: `curation_scripts/melanin_myco_step5/submit_v2.sh`. Inputs: `other_inputs/tree_inputs/{scytalone_dehydratase_v2,omt_v2}/` (new); tyrosinase, hppd and laccase_diphenol reuse the v1 input folders.

Changes to shared code (default behaviour unchanged when the new variables are unset): `slurm/step5_prepare_ident.sbatch` and `slurm/step5_multi_ref_clades.sbatch` read `PC_REFS_FA` (reference FASTA override; `refs/` is not changed) ; `step5_prepare_ident.sbatch` reads `PC_MIN_COV` (mmseqs `-c`, default 0.8); both place scripts read `PC_EPA_ARGS` (extra epa-ng arguments). New: `bin/step5_compare_calls.py`, `bin/step5_identity_distribution.py`, `tests/test_step5_v2_tools.py` (7 tests).

| run | prepare job (CPUs, elapsed) | clades job | place job (CPUs, elapsed) |
|---|---|---|---|
| scytalone_dehydratase_v2 (id 0.8, cap 300) | 32439796 (2, 00:25:20) | 32439797 COMPLETED | 32439798 (2, 00:02:01) |
| omt_v2 (id 0.7, cap 300, Phaffia alone) | 32439799 (2, 01:06:32) | 32439800 COMPLETED | 32439801 (2, 00:05:17) |
| tyrosinase_v2 (id 0.5, cap 1000) | 32439802 (8, 07:53:57) | 32439803 COMPLETED | 32439804 (8, 00:32:04; parallel hmmalign) |
| hppd_v2 (id 0.7, cap 1000 = all 662 clusters) | 32439805 (4, 02:12:10) | 32439806 COMPLETED | 32439807 (4, 00:13:42) |
| laccase_diphenol_v2 (id 1.0, cov 1.0) | 32441171 (2, 00:03:15) | 32441172 COMPLETED | 32441173 (2, 00:00:36) |
| laccase_diphenol_v2_epafilter (place only, epa-ng `--filter-max 10000 --filter-acc-lwr 0.999`) | copy of v2 | | 32442349 (2, 00:00:37) |

CPUs: at most 16 of mine at once (2 + 2 + 8 + 4). raxml-ng times (MEASURED, `ref.raxml.log`, `pars{2}`, two searches run one after the other): tyrosinase v1 310 taxa x 2,241 sites 2,716 s at 4 threads; tyrosinase v2 1,010 taxa x 3,839 sites 27,983 s (7 h 46 min) at 8 threads; hppd v1 305 x 800 822 s at 4 threads; hppd v2 667 x 936 7,191 s at 4 threads; omt v2 311 x 626 3,107 s at 2 threads. The estimate made before the run (taxa squared times sites, UNVERIFIED) was 8 to 13 h for tyrosinase at 1,000 representatives; the run took 7.8 h. A cap of 1,500 would be about 2.2 times longer by the same scaling (about 17 h, UNVERIFIED); it was not run.
Failed attempts (kept in the job file): epa-ng `--filter-min-lwr` with `--filter-acc-lwr` is refused (32441967); `--filter-acc-lwr 1.0` wrote an edge number that `upo_place assign` could not read (32442254, KeyError). `--filter-acc-lwr 0.999` worked (32442349). The coordinator message asked for `1.0`; it does not work here.
Tests (MEASURED): `/usr/bin/python3.12 -m pytest analysis/pathway_census/tests -q` (modules samtools, zstd, hmmer, mafft loaded): 931 passed, 1 failed, 1 skipped. The failure is the known `test_step5_place.py::test_extract_end_to_end_with_failures`.
