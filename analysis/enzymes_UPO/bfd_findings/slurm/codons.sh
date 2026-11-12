#!/bin/bash
# BFD start/stop flags vs genome sequence (bfd_findings/REPORT.md, finding 1).
#   sbatch analysis/enzymes_UPO/bfd_findings/slurm/codons.sh
# Reads the shared BFD tables and input assemblies read-only; writes bfd_findings/out/codon_check.tsv.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 3:00:00
#SBATCH -J upo_codons
#SBATCH -o logs/codons_%j.out
#SBATCH -e logs/codons_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
D=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full
OUT=$UPO/bfd_findings/out; mkdir -p $OUT
W="${SCRATCH:?}/codons"; mkdir -p $W
module load duckdb/1.1.3
$PY -W ignore -m unittest discover -s $UPO/bfd_findings/test
sed -e "s#\${IDS}#$D/pf01328_ids.tsv#" -e "s#\${OUT}#$W#" $UPO/bfd_findings/sql/codon_sets.sql > $W/q.sql
sed -i "1i SET memory_limit='40GB'; SET threads=8; SET temp_directory='$W/duck_tmp';" $W/q.sql
( cd $RUNS/tables && duckdb < $W/q.sql > /dev/null )
wc -l $W/codon_*.tsv
$PY $UPO/bfd_findings/bin/check_codons.py --cds $W/codon_cds.tsv --models $W/codon_models.tsv \
    --genomes-dir $RUNS/input_clean_genomes --proteins $D/all_pf01328.faa --threads 8 \
    --out $OUT/codon_check.tsv
echo "done $(date +%T)"
