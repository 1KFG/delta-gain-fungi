#!/bin/bash
#SBATCH -p short
#SBATCH --constraint="ryzen|cascade|broadwell"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=gl_null
#SBATCH --output=gl_null_%A_%a.log
# One task of the tree-aware simulation null (bin/run_null_sim.R): NREP replicates per array task.
# Required env: PC_ROOT (analysis/pathway_census), GL_PRIMARY (dir with primary.rds and calls.tsv), NULL_OUT.
# Optional env: NULL_NREP (default 250), NULL_ROOT (stationary|asr|flat; default stationary), NULL_TOL (prevalence window; 0 = off).
# No BASH_SOURCE is used. Task id = SLURM_ARRAY_TASK_ID.
set -euo pipefail
module unload R 2>/dev/null || true
module load R/4.2.0
Rscript "${PC_ROOT:?}/bin/run_null_sim.R" --primary "${GL_PRIMARY:?}" --outdir "${NULL_OUT:?}" --nrep "${NULL_NREP:-250}" \
  --task "${SLURM_ARRAY_TASK_ID:?}" --root "${NULL_ROOT:-stationary}" --match_tol "${NULL_TOL:-0}"
