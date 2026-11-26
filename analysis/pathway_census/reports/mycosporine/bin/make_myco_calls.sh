#!/usr/bin/env bash
# Make the two calls tables (hits and placed) for the mycosporine steps and the presence figure.
# Run from analysis/pathway_census. Required environment (no BASH_SOURCE use):
#   R1  results folder holding search/melanin_myco_full (searched.chunk*.tsv, candidates)  [melanin worktree]
#   R2  results folder holding genomes.tsv.zst and step5/mycosporine/*                      [step5b worktree]
#   SPECIES  species.parquet of the BFD tables
#   TMP  scratch folder for the intermediate tables
set -euo pipefail
: "${R1:?}" "${R2:?}" "${SPECIES:?}" "${TMP:?}"
mkdir -p "$TMP"
zstd -dc "$R2/genomes.tsv.zst" | /usr/bin/python3.12 -c "
import sys, csv
print('LOCUSTAG\tGENUS\tCLASS\tbusco_complete\tpep_path')
for x in csv.DictReader(sys.stdin, delimiter='\t'):
    if x['has_proteins'] == 'True':
        print('\t'.join([x['locustag'], x['genus'], x['class_'], x['busco_complete'], 'P']))
" > "$TMP/reps_all.tsv"
printf 'step\tsource\nddgs\thits\nomt\thits\natp_grasp_ligase\thits\n' > "$TMP/cfg_hits.tsv"
printf 'step\tsource\nddgs\t%s\nomt\t%s\n' \
  "$R2/step5/mycosporine/ddgs_allplace/assigned.PROVISIONAL.tsv" "$R2/step5/mycosporine/omt_v2/assigned.PROVISIONAL.tsv" > "$TMP/cfg_placed.tsv"
for m in hits placed; do
  /usr/bin/python3.12 bin/make_calls_table.py --reps "$TMP/reps_all.tsv" --taxonomy "$SPECIES" \
    --search-dir "$R1/search/melanin_myco_full" --config "$TMP/cfg_$m.tsv" \
    --prefix-aliases prefix_aliases.tsv --out "$TMP/$m/calls.tsv"
done
/usr/bin/python3.12 reports/mycosporine/figures/bin/fig_presence_myco.py \
  --calls-hits "$TMP/hits/calls.tsv" --calls-placed "$TMP/placed/calls.tsv" \
  --out-dir reports/mycosporine/figures/out
