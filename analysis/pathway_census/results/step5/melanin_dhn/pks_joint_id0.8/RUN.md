# Census step 5, rerun of the PKS joint tree at mmseqs identity 0.8, with bootstrap

Date: 2026-11-24 to 2026-11-25. Branch `pathway-census-melanin-myco`. All results are PROVISIONAL (`pathways/melanin_dhn.yaml` has `verified: false`).
Labels: MEASURED = read from a file named here. PROPOSED = a choice. UNVERIFIED = not checked.
Previous run (untouched, except the added sub-folder `per_step_decoys/`): `../pks_joint/` (identity 0.7, commit 40ebc4e). Same inputs (copies in `inputs/`).
The target clades are tree regions bounded by labelled leaves. They are not ortholog sets. No outgroup is in the tree (see the note, section "Outgroup: return to this").

## Differences from the 0.7 run

| Item | 0.7 run | This run (MEASURED) |
|---|---|---|
| mmseqs identity, coverage | 0.7, 0.8 | 0.8, 0.8 (`--cov-mode 0`) |
| Clusters = representatives | 237 | 532 of 9,374 domains; `PC_MAX_REPS=100000`, so none was drawn (`reps_summary.json`: `capped` false) |
| Leaves | 264 | 559 (4 refs + 5 decoys + 18 named + 532 reps) |
| Alignment | 4,090 to 3,658 columns | mafft `--auto` 7.505 (under 5 min, so no `--retree`); clipkit smart-gap 4,747 to 4,118 columns (13.25% removed); gaps 49.86%; 3,888 patterns |
| Tree | 2 starts, 4 threads | raxml-ng 2.0.4, LG+G4, `--tree pars{2}`, seed 1, 2 workers x 2 threads; best log-likelihood -1,098,828.08 (start 1) and -1,098,842.13 (start 2); alpha 0.891; 1 near-zero branch; 25,358 s |
| Decoy scope | union (every decoy bounds every step) | per step (owner decision 2026-11-24); union kept in `union_scope/` |
| Bootstrap | none | 100 replicates, NNI-only search (see below) |

## Commands

```
cd analysis/pathway_census; source results/step5/melanin_dhn/pks_joint_id0.8/env.sh   # PC_* variables, PC_MIN_SEQ_ID=0.8
sbatch --export=ALL -D $PC_OUT slurm/step5_prepare_joint.sbatch      # job 32436765: extract, cluster, select-reps, mafft, clipkit; CANCELLED during raxml-ng (see Jobs)
sbatch --export=ALL -D $PC_OUT slurm/step5_tree_joint.sbatch         # job 32436974: raxml-ng ML tree (-c 4)
# clades, decoy scope per step (the tree job called clades from an older spooled script and stopped with rc 1; clades were run by hand):
bin/step5_joint.py clades --decoy-scope per-step --step-decoys $PC_STEP_DECOYS --tree $PC_OUT/ref.raxml.bestTree --refs $PC_REFS_STEPS \
   --decoys $PC_DECOYS_FA --members-fa $PC_MEMBERS_FA --reps $PC_OUT/rep_ids.txt --joint-id pks_joint --out-dir $PC_OUT
sbatch --export=ALL -D $PC_OUT -c 4 slurm/step5_place_joint.sbatch   # job 32453905
sbatch --export=ALL,PC_BS_PER_TASK=1,PC_BS_DIR=bootstrap_nni,PC_BS_EXTRA="--opt-topology nni" -D $PC_OUT --array=1-100%4 -c 2 slurm/step5_bootstrap_joint.sbatch
sbatch --export=ALL,PC_BS_DIR=bootstrap_nni -D $PC_OUT slurm/step5_support_joint.sbatch   # job 32453916: --bs-metric fbp,tbe on the best ML tree
bin/step5_bootstrap.py support|stability|boundary|compare|rep-coverage ...     # tables in bootstrap_report/ and sampling_check/
sbatch --export=ALL,PC_RARE_SEEDS="1 2 3" -D $PC_OUT slurm/step5_rarefaction_joint.sbatch  # job 32455520
PC_RUN=$PC_OUT PC_NEWDIR=union_scope PC_DECOY_SCOPE=union bash slurm/step5_reapply_decoy_scope.sh   # the union rule from the saved placement
```

## Jobs (partition `intel`; at most 16 CPUs of mine at once; `-c 1` tasks are given 2 CPUs on these nodes)

| Job | What | State | Elapsed |
|---|---|---|---|
| 32436765 | prepare (extract 9,374 of 9,374 cut, mmseqs, mafft, clipkit, raxml-ng) | CANCELLED by me at about 12 min, during raxml-ng (to set the tree and bootstrap CPUs apart) | |
| 32436974 | ML tree (`step5_tree_joint.sbatch`) | tree COMPLETED; job state FAILED (rc 1): the clades call of the older spooled script had no `--step-decoys` | 7 h 02 min 41 s |
| 32437015 | standard bootstrap array, 2 threads per task | tasks 1 and 2 COMPLETED (14,782 s and 16,964 s per replicate); tasks 3 to 6 cancelled after 1 h 45 min; 7 to 100 cancelled before start | |
| 32439665, 32439666 | pilots: `--opt-topology simplified --stop-rule kh-mult` (the `--fast` preset) and `--opt-topology nni`, 2 tasks each | fast pilot cancelled after about 58 min with no replicate finished (its data are deleted); nni pilot: replicates 1 and 2 done (514 s, 473 s) | |
| 32440121 | NNI bootstrap array, tasks 3 to 100 | COMPLETED | per replicate: mean 453 s, range 337 to 588 s (100 replicates, 2 threads; `bootstrap_nni/bs_*.time`) |
| 32453905 | place (`-c 4`): hmmbuild, hmmalign, epa-ng (145 s), assign, report | COMPLETED | 50 min 37 s |
| 32453906 | support | CANCELLED by me: it read the wrong bootstrap folder; the sbatch got `PC_BS_DIR` | |
| 32453916 | support, FBP and TBE on the 100 NNI trees | COMPLETED (the final `ls` of the first script version named a wrong file suffix; fixed) | |
| 32455492 | rarefaction | FAILED after 6 s: `epa-ng -w` folder missing; fixed (`mkdir -p`) | |
| 32455520 | rarefaction, 2 fractions x 3 seeds | COMPLETED | 14 min 35 s |

## Bootstrap method and its limit (read before using the support values)

- The standard raxml-ng bootstrap replicate (adaptive search, 2 threads) took 14,782 s and 16,964 s (4.1 h and 4.7 h; MEASURED, 2 replicates). 100 replicates with at most 12 CPUs would take about 70 to 80 h (arithmetic from the 2 measured times, UNVERIFIED as a total). This is more than the 6 h limit. Fewer standard replicates (50) would still take about 35 to 40 h.
- PROPOSED: the 100 replicates use an NNI-only topology search (`raxml-ng --bootstrap --opt-topology nni`). One replicate took 337 to 588 s.
- Check on the same two resampled alignments (seeds 1001 and 1002; MEASURED, `bootstrap/`, `bootstrap_nni/`): the standard and the NNI tree differ by relative Robinson-Foulds distance 0.113 and 0.129 (126 and 144 differing splits; 556 splits per tree). Standard tree against the best ML tree: 0.131 and 0.171. NNI trees against the best ML tree: mean 0.199 over 100. The NNI log-likelihoods are lower by 727 and 522 than the standard ones for the same alignments (-1,110,883 against -1,110,156; -1,092,005 against -1,091,483). So NNI trees are less optimised than standard ones. The effect on each support value is UNVERIFIED; two replicates cannot measure it. The NNI search is likely to give lower support, not higher (PROPOSED reading, UNVERIFIED).
- `--bs-trees autoMRE` was not run: it needs the standard search.
- The two standard replicates give, for the clade rule: thn_synthase, hexaketide_pks and heptaketide_pks (per-step) each have a target in 2 of 2; wA is inside the smallest side of Alb1 and AlbA in 2 of 2 (`bootstrap_report/stability_std2_per-step.tsv`).

## Checks (MEASURED)

- 9,374 domains cut, 0 failures (`extract_summary.json`). `assign.<S>.tsv` has 9,374 queries per step (checked by the sbatch). Candidate proteins equal placements per step: 7,286, 8,604, 2,200. Joint candidate rows 9,374 = domains placed.
- Proteins that are `target` for more than one step: 0 (per-step and union).
- Unit tests: `tests/test_step5_joint.py`, `tests/test_step5_bootstrap.py`; full suite in the note.

## Files (small summary files are in git; `results/` is ignored)

`SAMPLING_CRITERION.txt` (written before the call comparison), `sampling_check/` (`rep_coverage.json`, `calls_*.json/tsv`), `rarefaction/*/compare_vs_full.*`, `bootstrap_report/` (`support_nodes_*.tsv`, `stability_*.tsv`, `boundary.*.tsv`), `clade_report.json`, `joint_report/`, `union_scope/`.
