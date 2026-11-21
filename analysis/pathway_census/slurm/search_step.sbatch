#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=pc_search
#SBATCH --output=pc_search_%j.log
# Search one protein FASTA with one step HMM. All paths come from the environment (no BASH_SOURCE).
# Required env: PC_ROOT (analysis/pathway_census), PC_PATHWAY, PC_STEP, PC_HMM_DIR, PC_PROTEINS,
#               PC_OUT_DIR, PC_MIN_SCORE. Optional: PC_PROVISIONAL=1.
# Usage: sbatch --export=ALL,PC_ROOT=...,PC_STEP=... slurm/search_step.sbatch
set -euo pipefail
: "${PC_ROOT:?}" "${PC_PATHWAY:?}" "${PC_STEP:?}" "${PC_HMM_DIR:?}" "${PC_PROTEINS:?}" "${PC_OUT_DIR:?}" "${PC_MIN_SCORE:?}"
SCR="${SCRATCH:?SCRATCH not set}"
module load hmmer/3.4
EXTRA=()
[ "${PC_PROVISIONAL:-0}" = "1" ] && EXTRA+=(--provisional)
# Stage the protein file on node-local scratch, write results to PC_OUT_DIR on shared storage.
cp "$PC_PROTEINS" "$SCR/"
/usr/bin/python3.12 "$PC_ROOT/bin/search_step.py" --pathway "$PC_PATHWAY" --step "$PC_STEP" \
  --hmm-dir "$PC_HMM_DIR" --proteins "$SCR/$(basename "$PC_PROTEINS")" --out-dir "$SCR/out" \
  --min-score "$PC_MIN_SCORE" --cpus "${SLURM_CPUS_PER_TASK:-8}" "${EXTRA[@]}"
mkdir -p "$PC_OUT_DIR"
zstd -T0 -q --rm "$SCR"/out/*.domtbl 2>/dev/null || true
cp "$SCR"/out/* "$PC_OUT_DIR"/
