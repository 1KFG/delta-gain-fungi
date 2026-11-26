#!/bin/bash
# Regenerate tables and the generated blocks of REPORT.md. Run from analysis/pathway_census.
set -euo pipefail
M=${M:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_melanin/analysis/pathway_census/results}
S5=${S5:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_step5b/analysis/pathway_census/results}
SP=${SP:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/species.parquet}
OUT=${OUT:?set OUT to a scratch dir for genome_calls.tsv}
/usr/bin/python3.12 reports/melanin/bin/build_genome_calls.py --species-parquet "$SP" --steps reports/melanin/bin/steps.tsv \
  --root mel="$M" --root s5b="$S5" --search-dir "$M/search/melanin_myco_full" --prefix-aliases prefix_aliases.tsv --out "$OUT/genome_calls.tsv"
/usr/bin/python3.12 reports/melanin/bin/summarise_calls.py --calls "$OUT/genome_calls.tsv" --out-dir reports/melanin/tables
/usr/bin/python3.12 reports/melanin/bin/make_report_tables.py --report reports/melanin/REPORT.md --tables reports/melanin/tables \
  --pathways pathways --min-scores "$M/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv" \
  --search-summary "$M/search/melanin_myco_full/summary.txt" --summary-v1 "$M/step5/other_steps_summary.tsv" \
  --summary-v2 "$M/step5/other_steps_summary.v2.tsv" --summary-allplace "$M/step5/other_steps_summary.allplace.tsv" \
  --pks-calls "$M/step5/melanin_dhn/pks_joint_id0.8/joint_report/calls_per_step.tsv"
/usr/bin/python3.12 reports/melanin/figures/bin/fig_melanin.py --tables reports/melanin/tables --out-dir reports/melanin/figures/out
