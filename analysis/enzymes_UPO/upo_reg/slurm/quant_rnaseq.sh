#!/bin/bash
#SBATCH -p epyc -c 16 --mem 48g --time 8:00:00
#SBATCH -J upo_rnaseq -o logs/quant_rnaseq.%j.log -e logs/quant_rnaseq.%j.log
# G3 Phase B: quantify the A. alternata stress series against the BFD Y9 annotation.
#
# Reference is Alternaria_alternata_Y9 (BFD LOCUSTAG F75F20B0), chosen because its gene ids are
# the ones carrying our clade calls. The RNA-seq is strain Z7, so this is CROSS-STRAIN
# quantification. The salmon mapping rate is recorded per run and is the QC that decides whether
# the comparison is usable at all.
#
# NEVER use BASH_SOURCE on SLURM. UPO_ROOT is passed in via --export.
set -euo pipefail
: "${UPO_ROOT:?set UPO_ROOT}"
: "${SCRATCH:?no SCRATCH}"

module load sratoolkit/3.2.0
module load parallel-fastq-dump/0.6.7
module load fastp/0.23.2
module load salmon/1.7.0

REF=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation/Alternaria_alternata_Y9/predict_results/Alternaria_alternata_Y9.cds-transcripts.fa.gz
MAN="$UPO_ROOT/data/g3/rnaseq/runs.tsv"
SRA="$UPO_ROOT/data/g3/rnaseq/sra"
QDIR="$UPO_ROOT/data/g3/rnaseq/quant"
mkdir -p "$QDIR"

WORK="$SCRATCH/upo_rnaseq"
mkdir -p "$WORK"
cp "$REF" "$WORK/ref.fa.gz"

echo "=== building salmon index ==="
salmon index -t "$WORK/ref.fa.gz" -i "$WORK/idx" -k 31 -p 16 --gencode 2>&1 | tail -3

while IFS=$'\t' read -r run cond nreads proj; do
  [ -z "$run" ] && continue
  if [ -s "$QDIR/$run/quant.sf" ]; then echo "have $run"; continue; fi
  f="$SRA/$run/$run.sra"
  if [ ! -s "$f" ]; then echo "MISSING sra $run"; continue; fi
  echo "=== $run ($cond) ==="
  rm -rf "$WORK/fq"; mkdir -p "$WORK/fq"
  parallel-fastq-dump --sra-id "$f" --threads 8 --outdir "$WORK/fq" \
      --split-files --gzip --tmpdir "$SCRATCH"
  r1="$WORK/fq/${run}_1.fastq.gz"; r2="$WORK/fq/${run}_2.fastq.gz"
  if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then echo "FAILED dump $run"; continue; fi
  fastp -i "$r1" -I "$r2" -o "$WORK/fq/t1.fq.gz" -O "$WORK/fq/t2.fq.gz" \
        -w 8 -j "$QDIR/$run.fastp.json" -h /dev/null 2>&1 | tail -2
  salmon quant -i "$WORK/idx" -l A -1 "$WORK/fq/t1.fq.gz" -2 "$WORK/fq/t2.fq.gz" \
        -p 16 --validateMappings -o "$QDIR/$run" 2>&1 | tail -3
  rate=$(grep -o '"percent_mapped": *[0-9.]*' "$QDIR/$run/aux_info/meta_info.json" | head -1 | tr -d ' ' | cut -d: -f2)
  echo "MAPRATE $run $cond $rate"
  rm -rf "$WORK/fq"
done < "$MAN"

echo "=== mapping rates ==="
for d in "$QDIR"/SRR*/; do
  r=$(basename "$d")
  [ -s "$d/aux_info/meta_info.json" ] || continue
  m=$(grep -o '"percent_mapped": *[0-9.]*' "$d/aux_info/meta_info.json" | head -1 | cut -d: -f2)
  echo -e "$r\t$m"
done
echo "=== done ==="
