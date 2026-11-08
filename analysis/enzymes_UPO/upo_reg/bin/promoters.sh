#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -J upo_prom
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/prom_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/prom_%j.err
set -euo pipefail
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
G3=$PROJ/data/g3
module load duckdb samtools/1.22.1
export TMPDIR="${SCRATCH:?SCRATCH not set}/upo_prom.$SLURM_JOB_ID"; mkdir -p "$TMPDIR"

/usr/bin/python3.12 "$PROJ/upo_reg/bin/make_promoter_coords.py" \
  --targets "$G3/target_genes.txt" \
  --gene-info /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/gene_info.parquet \
  --upstream 1000 --min-len 300 --bg-per-genome 5 --batch 300 --memory 48GB \
  --out "$G3/promoter_coords.tsv"

/usr/bin/python3.12 "$PROJ/upo_reg/bin/fetch_promoter_seqs.py" \
  --coords "$G3/promoter_coords.tsv" \
  --genomes-tsv /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb/genomes.tsv \
  --outdir "$G3/seqs" --samtools "$(which samtools)" --min-len 300

rm -rf "$TMPDIR"
echo "promoters done $(date +%T)"
