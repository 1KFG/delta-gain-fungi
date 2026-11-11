#!/bin/bash
# Tables for bfd_findings/REPORT.md (gene models, glomeromycete clade, ascomycetes in Basidiomycetous 1st).
#   sbatch analysis/enzymes_UPO/bfd_findings/slurm/report.sh
# Reads the shared BFD tables, input assemblies and UPO outputs read-only; writes bfd_findings/out/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 3:00:00
#SBATCH -J upo_bfdfind
#SBATCH -o logs/bfdfind_%j.out
#SBATCH -e logs/bfdfind_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
D=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full
OUT=$UPO/bfd_findings/out; mkdir -p $OUT
W="${SCRATCH:?}/bfdfind"; mkdir -p $W
module load duckdb/1.1.3 zstd 2>/dev/null || module load duckdb/1.1.3
$PY -W ignore -m unittest discover -s $UPO/bfd_findings/test
sed -e "s#\${IDS}#$D/pf01328_ids.tsv#" -e "s#\${OUT}#$W#" $UPO/bfd_findings/sql/export.sql > $W/q.sql
sed -i "1i SET memory_limit='40GB'; SET threads=8; SET temp_directory='$W/duck_tmp';" $W/q.sql
( cd $RUNS/tables && duckdb < $W/q.sql )
wc -l $W/*.tsv
$PY $UPO/bfd_findings/bin/report_tables.py --placed $D/upo_place_bfd_all.tsv.zst \
    --models $W/upo_models.tsv --genome-qc $W/genome_qc.tsv --tax $W/genome_tax.tsv \
    --fasta $D/all_pf01328.faa --genomes-dir $RUNS/input_clean_genomes --threads 8 --out $OUT
cp $W/genome_qc.tsv $OUT/bfd_genome_gene_model_qc.tsv
echo "done $(date +%T)"
