#!/bin/bash
# Placement-only rerun of step 5 runs with epa-ng keeping all placements (PC_EPA_ARGS), into <run>_allplace folders.
# Trees, alignments, clusters and clade files are copied from the old run folder; nothing in the old folder is changed.
# Usage (paths from environment, no BASH_SOURCE): PC_ROOT=<worktree>/analysis/pathway_census R5=<results/step5 of the old runs>
#   submit_allplace.sh single <pathway> <step_folder> <cpus> [dependency_job_id]
# Prints the job id. The joint PKS trees are handled by submit_allplace_joint.sh.
set -euo pipefail
: "${PC_ROOT:?}" "${R5:?}"
kind=$1; pathway=$2; folder=$3; cpus=$4; dep=${5:-}
[ "$kind" = single ] || { echo "kind must be single" >&2; exit 2; }
OLD=$R5/$pathway/$folder; NEW=${OLD}_allplace
[ -d "$OLD" ] || { echo "no $OLD" >&2; exit 1; }
[ ! -e "$NEW" ] || { echo "$NEW exists; not overwritten" >&2; exit 1; }
mkdir -p "$NEW"
# inputs only (no placement-stage output)
rsync -a --exclude='assign*' --exclude='epa*' --exclude='placements*' --exclude='assigned*' --exclude='candidates_step*' \
  --exclude='ref.hmm' --exclude='hmm*' --exclude='SUMMARY.tsv' --exclude='report' --exclude='step5_*.out' --exclude='calls_*' "$OLD"/ "$NEW"/
# the step name is the step_id of the old run (folder names carry suffixes such as _v2, _id50)
STEP=$(awk -F'\t' '$1=="clades.step_id"{print $2}' "$OLD/SUMMARY.tsv")
[ -n "$STEP" ] || { echo "no clades.step_id in $OLD/SUMMARY.tsv" >&2; exit 1; }
export PC_ROOT PC_PATHWAY=$pathway PC_STEP=$STEP PC_OUT=$NEW PC_MIN_SUPPORT=0.5 \
  PC_CAND_DIR=${PC_CAND_DIR:?set PC_CAND_DIR} PC_EPA_ARGS="--filter-acc-lwr 0.999 --filter-max 10000"
DEPARG=; [ -n "$dep" ] && DEPARG="--dependency=afterany:$dep"
sbatch --parsable --export=ALL -p intel -c "$cpus" -t 4:00:00 -D "$NEW" $DEPARG "$PC_ROOT/slurm/step5_place_parallel.sbatch"
