#!/bin/bash
#SBATCH -p short
#SBATCH --constraint="ryzen|cascade|broadwell"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=gl_runr
#SBATCH --output=gl_runr_%j.log
# Ancestral state reconstruction, event counts and co-event permutation tests (bin/run_gainloss_roots.R).
# Required env: PC_ROOT (analysis/pathway_census), GL_TREE (FastTree treefile), GL_OUT (output dir containing calls.tsv).
# Optional env: GL_ROOTS (comma list, first is primary; default Opisthosporidia), GL_NPERM (default 10000). No BASH_SOURCE is used.
set -euo pipefail
module unload R 2>/dev/null || true
module load R/4.2.0
Rscript "${PC_ROOT:?}/bin/run_gainloss_roots.R" --calls "${GL_OUT:?}/calls.tsv" --tree "${GL_TREE:?}" --outdir "$GL_OUT" --nperm "${GL_NPERM:-10000}" --roots "${GL_ROOTS:-Opisthosporidia}"
