#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -J upo_nr
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/nr_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/nr_%j.err
set -euo pipefail
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
OUT=$PROJ/data/bfd_full/novelty
NR=/srv/projects/db/ncbi/diamond/current/nr.dmnd
module load diamond/2.1.24
export TMPDIR="${SCRATCH:?SCRATCH not set}/upo_nr.$SLURM_JOB_ID"; mkdir -p "$TMPDIR"
echo "query: $(grep -c '^>' $OUT/reps_090.faa) representatives (90% derep of the catalytic core)"
diamond blastp --db "$NR" --query "$OUT/reps_090.faa" \
    --out "$TMPDIR/reps_vs_nr.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle \
    --max-target-seqs 5 --evalue 1e-5 \
    --threads "${SLURM_CPUS_PER_TASK:-32}" --tmpdir "$TMPDIR" --block-size 6 --index-chunks 1
zstd -qf "$TMPDIR/reps_vs_nr.tsv" -o "$OUT/reps090_vs_nr.tsv.zst"
rm -rf "$TMPDIR"
echo "done $(date +%T)"
