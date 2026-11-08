#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 48:00:00
#SBATCH -J upo_novelty
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/novelty_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/novelty_%j.err
set -euo pipefail

# G1 step 4: how novel are BFD UPOs against public protein space?
#
# Dereplicate first. Searching all 97,520 against a 523 GB nr is mostly redundant work:
# BFD contains many near-identical copies across closely related assemblies.
#
# KNOWN CONFOUND, stated in PLAN.md: absence from nr can mean "not deposited" rather than
# "novel". Many BFD genomes are themselves unpublished. A protein-only search therefore
# measures deposition status as much as novelty. The nucleotide/WGS comparison is a separate
# step and is required before any novelty claim is made.

PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
IN=$PROJ/data/bfd_full/all_pf01328.faa
OUT=$PROJ/data/bfd_full/novelty
NR=/srv/projects/db/ncbi/diamond/current/nr.dmnd
ID=${ID:-0.90}          # dereplication identity

module load mmseqs2/13-45111
module load diamond/2.1.24
mkdir -p "$OUT"
export TMPDIR="${SCRATCH:?SCRATCH not set}/upo_novelty.$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

if [ ! -s "$OUT/reps.faa" ]; then
  echo "=== dereplicate at $ID $(date +%T) ==="
  mmseqs easy-cluster "$IN" "$TMPDIR/clu" "$TMPDIR/tmp" \
      --min-seq-id "$ID" -c 0.8 --cov-mode 0 --threads "${SLURM_CPUS_PER_TASK:-32}"
  cp "$TMPDIR/clu_rep_seq.fasta" "$OUT/reps.faa"
  cp "$TMPDIR/clu_cluster.tsv"   "$OUT/clusters.tsv"
  echo "representatives: $(grep -c '^>' "$OUT/reps.faa")"
fi

echo "=== diamond vs nr $(date +%T) ==="
diamond blastp \
    --db "$NR" --query "$OUT/reps.faa" \
    --out "$OUT/reps_vs_nr.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore staxids sscinames stitle \
    --max-target-seqs 5 --evalue 1e-5 \
    --threads "${SLURM_CPUS_PER_TASK:-32}" \
    --tmpdir "$TMPDIR" --block-size 6 --index-chunks 1
zstd -qf --rm "$OUT/reps_vs_nr.tsv"
zstd -qf --rm "$OUT/clusters.tsv" 2>/dev/null || true
rm -rf "$TMPDIR"
echo "done $(date +%T)"
