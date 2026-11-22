#!/bin/bash
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=pc_score
#SBATCH --output=pc_score_%j.log
# Step 2: score references and decoys against the step HMMs. No BASH_SOURCE; paths from the environment.
# Required env: PC_ROOT, PC_PATHWAY_ID, PC_OUT_DIR. Optional: PC_TAG (for example PROVISIONAL).
set -euo pipefail
: "${PC_ROOT:?}" "${PC_PATHWAY_ID:?}" "${PC_OUT_DIR:?}"
module load hmmer/3.4 mafft/7.505
/usr/bin/python3.12 "$PC_ROOT/bin/score_refs_decoys.py" \
  --refs-dir "$PC_ROOT/refs/$PC_PATHWAY_ID" \
  --hmm-dir "$PC_ROOT/results/hmms/$PC_PATHWAY_ID" \
  --out-dir "$PC_OUT_DIR" --tag "${PC_TAG:-}" --cpus "${SLURM_CPUS_PER_TASK:-4}"
