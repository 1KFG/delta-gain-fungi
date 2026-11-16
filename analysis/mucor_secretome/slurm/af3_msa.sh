#!/bin/bash
# Stage 1 of the AF3 screen: MSAs, once per chain. CPU only, so it starts immediately.
#   sbatch analysis/mucor_secretome/slurm/af3_msa.sh            # chitin panel + controls
#   sbatch --export=ALL,SET=unknowns analysis/mucor_secretome/slurm/af3_msa.sh
# AF3's data pipeline is jackhmmer against uniref90/MGnify/BFD and is the expensive half; an MSA
# is a property of a chain, not of a pair, so computing it once here and embedding it later
# turns ~90 pair predictions into ~20 jackhmmer runs. Everything runs inside the AF3 container
# against the local databases - no sequence leaves the cluster.
# Reads shared folders read-only; writes mucor_secretome/out/af3/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -J af3_msa
#SBATCH -a 0-15
#SBATCH -o logs/af3_msa_%A_%a.out
#SBATCH -e logs/af3_msa_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
SET=${SET:-chitin}
OUT=$A/out/af3/$SET; mkdir -p $OUT/jobs $OUT/msas
AF3=/opt/linux/rocky/8.x/x86_64/pkgs/alphafold/3.0.2/alphafold
module load alphafold/3.0.2

TASK=${SLURM_ARRAY_TASK_ID:-0}
NTASK=${SLURM_ARRAY_TASK_COUNT:-1}

# AF3 3.0.2 writes its output under the job name AS GIVEN, e.g. pred/CTRL_4EBZ_03a7bf__CO4/
# CTRL_4EBZ_03a7bf__CO4_summary_confidences.json. These scripts assumed it lowercased the name.
# It does not, so every "is this already done?" test silently failed: nothing was ever skipped,
# every rerun refolded the whole panel from scratch, and the completeness gate never fired, so
# the report was never written even when all 196 and all 426 predictions existed. Accept either
# spelling, the way bin/af3_input.py's find_msa already does.
af3_done() {   # af3_done <dir> <name> <suffix>
  local d=$1 n=$2 sfx=$3
  [ -s "$d/$n/${n}${sfx}" ] || [ -s "$d/${n,,}/${n,,}${sfx}" ] || [ -s "$d/${n,,}/${n}${sfx}" ]
}
# Task 0 builds the job JSONs; the others wait for them, so one submission does both.
if [ "$TASK" = "0" ]; then
  $PY -W ignore -m unittest discover -s $A/test
  $PY $A/bin/af3_input.py single --fasta $OUT/panel.faa --out $OUT/jobs
  touch $OUT/jobs/.ready
else
  for _ in $(seq 1 60); do [ -f $OUT/jobs/.ready ] && break; sleep 10; done
fi

mapfile -t JOBS < <(ls $OUT/jobs/*.json | sort)
echo "task $TASK/$NTASK of ${#JOBS[@]} jobs"
for i in "${!JOBS[@]}"; do
  [ $(( i % NTASK )) -eq "$TASK" ] || continue
  j="${JOBS[$i]}"
  n=$(basename "$j" .json)
  # resumable: a chain whose data JSON exists is already done
  if af3_done "$OUT/msas" "$n" _data.json; then echo "skip $n"; continue; fi
  echo "=== $n $(date +%T)"
  $AF3 --json_path="$j" --output_dir="$OUT/msas" --db_dir="$ALPHAFOLD_DB" \
       --model_dir="$ALPHAFOLD_DB/model" --norun_inference --jackhmmer_n_cpu=8 \
       --force_output_dir || echo "FAILED $n" >&2
done
echo "done $(date +%T)"
