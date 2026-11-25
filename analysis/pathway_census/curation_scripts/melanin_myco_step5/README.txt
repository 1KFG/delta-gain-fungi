One-off helpers of the step 5 run for the other melanin and mycosporine steps (2026-11-24, PROVISIONAL).
Order (from analysis/pathway_census; copies of the scripts also sit in results/step5/other_inputs/, where the run used them):
  fetch_family.py        UniProt reviewed entries of one Pfam family (network)
  screen:                bin/screen_outgroups.py of each family FASTA against the step HMM; bin/filter_decoys_for_tree.py for the decoys
  make_leaves.sh         bin/select_family_leaves.py per step
  assemble_inputs.sh     bin/cut_leaves_to_domain.py; decoys, members, outgroups per step
  assemble_sdr_joint.sh  joint t4hn_reductase + t3hn_reductase inputs
  submit_step.sh / submit_sdr_joint.sh / derive_sdr_step.sh   SLURM chains (prepare, multi_ref_clades, place)
  collate.py / collate_tables.py   tables from bin/step5_run_report.py
The shell scripts hold paths under results/step5/other_inputs/ and results/search/. They are a record, not a general tool.
