# Placement-only rerun of the PKS joint tree with all epa-ng placements: pks_joint_id0.8_allplace

Date: 2026-11-25. Old run (unchanged): `melanin_dhn/pks_joint_id0.8`. Label rule: MEASURED = read from a file or `sacct`. All results PROVISIONAL.
Tree, alignment, representatives and clusters were copied from the old run (no tree rebuilt; `domains.faa` was decompressed for the job and compressed again with zstd). `bin/step5_joint.py clades --decoy-scope per-step --step-decoys <refs/melanin_dhn/decoys/STEP.faa>` was run again; the three `clade_members.<step>.tsv` files are byte-identical to the old per-step files (`cmp`, MEASURED). Then `slurm/step5_place_joint.sbatch` with `PC_EPA_ARGS="--filter-acc-lwr 0.999 --filter-max 10000"` (new optional variable; unset = unchanged). epa-ng 0.3.8.
Submit script: `curation_scripts/melanin_myco_step5/submit_allplace_joint.sh`. Job 32461567, `intel`, 4 CPUs, elapsed 00:48:15 (MEASURED, `sacct`).
Result: `joint_report/calls_per_step.tsv` equals the old per-step table; 0 proteins changed call in all 3 steps (note `docs/superpowers/notes/2026-11-25-epa-all-placements.md`).
