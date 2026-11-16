#!/bin/bash
# Does the enriched-Pfam novelty survive counting Y1000+ as already public?
#
#   sbatch analysis/delta_gain_funnel/slurm/y1000_join.sh
#
# Runs the join over the enriched Pfam families (the 204 run in family_null_order_wide) and
# over the size-matched random-control families, which is the comparison that matters: a
# shrinkage that hits both lists equally is a property of the corpus, while one that hits
# only the enriched list is the deposition artefact the funnel's limitations warned about.
#
# Resources: the work is a per-family duckdb join over ~2.1 M regions, done one family at a
# time. Reading all 204 families at once OOM-killed a login node, which is why this is a
# batch job and why the script streams. 64 GB and 8 threads is comfortable.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -N 1 -n 1 -c 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -J y1000_join
#SBATCH -o analysis/delta_gain_funnel/logs/y1000_join_%j.out
#SBATCH -e analysis/delta_gain_funnel/logs/y1000_join_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

A=analysis/delta_gain_funnel
PY=/usr/bin/python3.12
R=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/results
SAMPLES=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/samples.parquet
OUT=$A/out
mkdir -p "$OUT"

# PRJNA736342 "Y1000+ Saccharomycotina yeasts genome sequencing and assembly" (UW-Madison,
# registered 2021) and PRJNA429441 "Sequencing of 196 yeast species from the subphylum
# Saccharomycotina" (UW-Madison, 2018), the precursor set. Both published proteomes through
# the project rather than through NCBI's annotation pipeline, so every genome in them is on
# the funnel's "new" side while its proteins have been public for years.
BP=(--bioproject PRJNA736342 --bioproject PRJNA429441)

run () {   # run <label> <null-dir> <regions-dir> <clusters-dir> <out-basename>
  echo "=== $1 $(date +%T)"
  $PY $A/bin/y1000_join.py \
      --null-dir "$R/$2" --regions-dir "$R/$3" --clusters-dir "$R/$4" \
      --samples "$SAMPLES" "${BP[@]}" \
      --threads 8 --memory 48GB --label "$1" --out "$OUT/$5.tsv"
}

run "enriched Pfam"        family_null_order_wide family_regions_wide family_clusters_wide y1000_enriched
run "random-control Pfam"  family_null_order_ctrl family_regions_ctrl family_clusters_ctrl y1000_control

echo "done $(date +%T)"
