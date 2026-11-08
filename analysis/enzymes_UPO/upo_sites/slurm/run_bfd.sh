#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 4:00:00
#SBATCH -J upo_sites
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/upo_sites_%A_%a.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/upo_sites_%A_%a.err
set -euo pipefail

# Absolute paths only. Never $(dirname "${BASH_SOURCE[0]}") -- it resolves wrong under SLURM.
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
CHUNKDIR=$PROJ/data/bfd_full/chunks
OUTDIR=$PROJ/data/bfd_full/out
HMM=$PROJ/upo_id/models/PF01328.hmm
BACKBONE=$PROJ/data/clade_backbone/raw.fa
PY=/usr/bin/python3.12

# Chunks per array task. Sized so one task runs roughly 1-1.5 h; see PLAN.md section 5 step 7.
PER_TASK=${PER_TASK:-5}

module load hmmer/3.4
module load zstd 2>/dev/null || true
mkdir -p "$OUTDIR"

TASK=${SLURM_ARRAY_TASK_ID:-1}
mapfile -t ALL < <(ls "$CHUNKDIR"/chunk_*.fa | sort)
START=$(( (TASK - 1) * PER_TASK ))
if [ "$START" -ge "${#ALL[@]}" ]; then
  echo "task $TASK: nothing to do (${#ALL[@]} chunks total)"; exit 0
fi
MINE=( "${ALL[@]:$START:$PER_TASK}" )
echo "task $TASK handling ${#MINE[@]} chunks: ${MINE[*]##*/}"

# node-local scratch for the hmmalign temporaries
export TMPDIR="${SCRATCH:?SCRATCH not set}/upo_sites.$TASK"
mkdir -p "$TMPDIR"

for f in "${MINE[@]}"; do
  base=$(basename "$f" .fa)
  out="$OUTDIR/${base}.tsv.zst"
  if [ -s "$out" ]; then echo "skip $base (exists)"; continue; fi
  echo "=== $base $(date +%T) ==="
  $PY "$PROJ/upo_sites/bin/upo_sites.py" \
      --proteome "$f" \
      --hmm "$HMM" \
      --backbone "$BACKBONE" \
      --pre-selected \
      --chunk 4000 \
      --cpu "${SLURM_CPUS_PER_TASK:-8}" \
      --out "$out"
  echo "=== $base done $(date +%T) ==="
done

rm -rf "$TMPDIR"
echo "task $TASK complete"
