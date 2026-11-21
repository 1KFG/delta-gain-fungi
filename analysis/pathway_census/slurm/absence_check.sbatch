#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=pc_absence
#SBATCH --output=pc_absence_%j.log
# tblastn absence check for a chunk of genomes. All paths come from the environment (no BASH_SOURCE).
# Required env: PC_ROOT (analysis/pathway_census), PC_PATHWAY, PC_MATRIX, PC_GENOMES, PC_SEQDB_GENOMES,
#               PC_REFS_DIR, PC_OUT_DIR, PC_MIN_BITSCORE.
# Optional env: PC_MIN_QUERY_COVER (default 0.7, PROPOSED), PC_LOCUSTAGS (file of genomes for this chunk),
#               PC_PROVISIONAL=1.
# Sizing: UNVERIFIED. Aim for 1-1.5 h of real runtime per job (chunk genomes with PC_LOCUSTAGS). No job has
# run yet, so the genomes-per-chunk number is unknown. Measure one real job, then set the chunk size.
# Output names carry .chunk<ID> (SLURM_ARRAY_TASK_ID, else the PC_LOCUSTAGS basename stem), so chunks do not overwrite.
# Genome FASTA .gz files are streamed into $SCRATCH by the library and removed after each search.
set -euo pipefail
: "${PC_ROOT:?}" "${PC_PATHWAY:?}" "${PC_MATRIX:?}" "${PC_GENOMES:?}" "${PC_SEQDB_GENOMES:?}" \
  "${PC_REFS_DIR:?}" "${PC_OUT_DIR:?}" "${PC_MIN_BITSCORE:?}"
SCR="${SCRATCH:?SCRATCH not set}"
module load ncbi-blast/2.17.0+
EXTRA=()
[ "${PC_PROVISIONAL:-0}" = "1" ] && EXTRA+=(--provisional)
[ -n "${PC_LOCUSTAGS:-}" ] && EXTRA+=(--locustags "$PC_LOCUSTAGS")
/usr/bin/python3.12 "$PC_ROOT/bin/absence_check.py" --pathway "$PC_PATHWAY" --matrix "$PC_MATRIX" \
  --genomes "$PC_GENOMES" --seqdb-genomes "$PC_SEQDB_GENOMES" --refs-dir "$PC_REFS_DIR" \
  --out-dir "$SCR/out" --min-bitscore "$PC_MIN_BITSCORE" \
  --chunk-id "${SLURM_ARRAY_TASK_ID:-}" --min-query-cover "${PC_MIN_QUERY_COVER:-0.7}" --cpus "${SLURM_CPUS_PER_TASK:-8}" "${EXTRA[@]}"
mkdir -p "$PC_OUT_DIR"
cp "$SCR"/out/* "$PC_OUT_DIR"/
