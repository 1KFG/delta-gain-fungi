#!/bin/bash
#SBATCH -p short
#SBATCH --constraint="ryzen|cascade|broadwell"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=gl_pagel
#SBATCH --output=gl_pagel_%A_%a.log
# One Pagel test per array task (bin/run_pagel.R). Task id i picks line i of $GL_PAIRS (two tab-separated step names).
# Required env: PC_ROOT, GL_OUT (dir with primary.rds and calls.tsv), GL_PAIRS (file of step pairs).
set -euo pipefail
module unload R 2>/dev/null || true
module load R/4.2.0
read -r A B < <(sed -n "${SLURM_ARRAY_TASK_ID:?}p" "${GL_PAIRS:?}")
mkdir -p "${GL_OUT:?}/pagel"
Rscript "${PC_ROOT:?}/bin/run_pagel.R" "$GL_OUT/primary.rds" "$A" "$B" "$GL_OUT/pagel/pagel.${A}__${B}.tsv"
