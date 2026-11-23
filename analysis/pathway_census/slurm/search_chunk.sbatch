#!/bin/bash
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=pc_chunk
#SBATCH --output=pc_chunk_%j.log
# Search one chunk of genomes with all step HMMs. Resources are UNVERIFIED until the timed job (step 3).
# No BASH_SOURCE. Required env: PC_ROOT, PC_PATHWAY, PC_HMM_DIR, PC_CHUNKS, PC_CHUNK_ID, PC_MIN_SCORES, PC_OUT_DIR.
# Optional: PC_PROVISIONAL=1.
set -euo pipefail
: "${PC_ROOT:?}" "${PC_PATHWAY:?}" "${PC_HMM_DIR:?}" "${PC_CHUNKS:?}" "${PC_CHUNK_ID:?}" "${PC_MIN_SCORES:?}" "${PC_OUT_DIR:?}"
: "${SCRATCH:?SCRATCH not set}"
module load hmmer/3.4
EXTRA=()
[ "${PC_PROVISIONAL:-0}" = "1" ] && EXTRA+=(--provisional)
/usr/bin/python3.12 "$PC_ROOT/bin/search_chunk.py" --pathway "$PC_PATHWAY" --hmm-dir "$PC_HMM_DIR" \
  --chunk-manifest "$PC_CHUNKS" --chunk-id "$PC_CHUNK_ID" --min-scores "$PC_MIN_SCORES" \
  --out-dir "$PC_OUT_DIR" --cpus "${SLURM_CPUS_PER_TASK:-16}" "${EXTRA[@]}"
