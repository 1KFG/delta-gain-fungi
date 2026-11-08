#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 2:00:00
#SBATCH -J upo_boltz
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/boltz_%A_%a.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/boltz_%A_%a.err
set -euo pipefail

# Model all 19 calibration enzymes with the heme in place.
#
# WHY: only 5 of the 19 have crystal structures, so the channel-width versus naphthalene test
# currently rests on n=5 (and on n=4 once TteUPO is set aside). All 19 have naphthalene
# measurements from Figure S4. Modelling the other 14 takes the test to n=19.
#
# The 5 with crystal structures are modelled TOO, deliberately: comparing predicted against
# experimental channel width on those 5 is the validation that decides whether a predicted
# structure can be measured this way at all. If it fails, the other 14 are not used.
#
# PYTHONNOUSERSITE=1 is required: ~/.local numpy 2.5.2 otherwise shadows the environment's
# own numpy 1.26.4 and boltz fails at import with a message blaming numba.

PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
BOLTZ=/opt/linux/rocky/8.x/x86_64/pkgs/boltz/2.2.1/env/bin/boltz
OUT=$PROJ/data/boltz/out
export BOLTZ_CACHE=/srv/projects/db/boltz
export PYTHONNOUSERSITE=1
mkdir -p "$OUT"

mapfile -t ALL < <(ls "$PROJ"/data/boltz/yaml/*.yaml | sort)
TASK=${SLURM_ARRAY_TASK_ID:-1}
PER=${PER_TASK:-3}
START=$(( (TASK-1) * PER ))
if [ "$START" -ge "${#ALL[@]}" ]; then echo "task $TASK: nothing to do"; exit 0; fi
MINE=( "${ALL[@]:$START:$PER}" )
echo "task $TASK: ${#MINE[@]} structures"

for y in "${MINE[@]}"; do
  b=$(basename "$y" .yaml)
  if [ -d "$OUT/boltz_results_$b" ]; then echo "skip $b"; continue; fi
  echo "=== $b $(date +%T) ==="
  $BOLTZ predict "$y" --out_dir "$OUT" --cache "$BOLTZ_CACHE" \
      --accelerator gpu --devices 1 --output_format pdb \
      --diffusion_samples 1 --override || echo "FAILED $b"
  echo "=== $b done $(date +%T) ==="
done
echo "task $TASK complete"
