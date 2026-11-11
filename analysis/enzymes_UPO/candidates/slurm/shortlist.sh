#!/bin/bash
# Expression shortlist of UPO candidates from under-characterized lineages.
#   sbatch analysis/enzymes_UPO/candidates/slurm/shortlist.sh
# Reads the shared BFD tables and UPO outputs read-only; writes candidates/out/ in this clone.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -J upo_shortlist
#SBATCH -o logs/shortlist_%j.out
#SBATCH -e logs/shortlist_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
TABLES=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
D=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full
OUT=$UPO/candidates/out; mkdir -p $OUT
W="${SCRATCH:?}/shortlist"; mkdir -p $W
module load duckdb/1.1.3 zstd 2>/dev/null || module load duckdb/1.1.3
$PY -W ignore -m unittest discover -s $UPO/candidates/test
sed -e "s#\${IDS}#$D/pf01328_ids.tsv#" -e "s#\${OUT}#$W#" $UPO/candidates/sql/export_qc.sql > $W/q.sql
( cd $TABLES && duckdb < $W/q.sql )
echo "geneqc rows: $(($(wc -l < $W/geneqc.tsv) - 1)); genomes: $(($(wc -l < $W/asm.tsv) - 1))"
$PY $UPO/candidates/bin/shortlist.py --placed $D/upo_place_bfd_all.tsv.zst \
    --scored $D/upo_sites_bfd_scored.tsv.zst --geneqc $W/geneqc.tsv --asm $W/asm.tsv \
    --clusters70 $D/novelty/clusters_070.tsv.zst --clusters90 $D/novelty/clusters_090.tsv.zst \
    --nr $D/novelty/reps090_vs_nr.tsv.zst --panel $D/assay_panel_96.tsv \
    --fasta $D/all_pf01328.faa --out $OUT
echo "done $(date +%T)"
