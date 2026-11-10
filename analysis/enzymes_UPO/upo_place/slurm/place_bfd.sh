#!/bin/bash
# Placement clade calls for the BFD PF01328 set: one array task per chunk.
#   sbatch --array=1-20 analysis/enzymes_UPO/upo_place/slurm/place_bfd.sh   (after prepare_ref.sh)
# Reads the BFD chunks (read-only) from the shared project; writes only to this clone:
#   analysis/enzymes_UPO/data/bfd_place/chunk_NNN.tsv  (merge with merge_bfd.sh / see README)
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 6:00:00
#SBATCH -J upo_place_bfd
#SBATCH -o logs/place_bfd_%A_%a.out
#SBATCH -e logs/place_bfd_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source analysis/enzymes_UPO/upo_place/slurm/env.sh
CHUNKDIR=${CHUNKDIR:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full/chunks}
OUT=$UPO/data/bfd_place; mkdir -p $OUT
mapfile -t ALL < <(ls $CHUNKDIR/chunk_*.fa | sort)
F=${ALL[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
B=$(basename $F .fa)
W="${SCRATCH:?}/place_$B"; mkdir -p $W
hmmalign --outformat afa $PREF/ref.hmm $F > $W/aln.afa
$PY $TOOL match-cols $W/aln.afa --length $NCOL --out $W/q.fa
$PX epa-ng --ref-msa $TREEDIR/aln.trim.afa --tree $PREF/ref.raxml.bestTree --query $W/q.fa \
    --model $PREF/ref.raxml.bestModel -w $W --redo -T ${SLURM_CPUS_PER_TASK:-4}
$PY $TOOL assign $W/epa_result.jplace --clades $TREEDIR/mfp.clade_members.tsv --out $OUT/$B.tsv
gzip -c $W/epa_result.jplace > $OUT/$B.jplace.gz
echo "done $B $(date +%T)"
