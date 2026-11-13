#!/bin/bash
# Stage 2 of the AF3 screen: inference only, reusing the stage-1 MSAs. GPU, resumable, shardable.
#   sbatch --array=0-7 analysis/mucor_secretome/slurm/af3_fold.sh          # chitin, 8 shards
#   sbatch --array=0-15 --export=ALL,SET=unknowns \
#       analysis/mucor_secretome/slurm/af3_fold.sh                         # unknowns, 16 shards
#   sbatch -p exfab --gres=gpu:ada6000:1 -t 1-00:00:00 --array=0-1 \
#       --export=ALL,SET=unknowns analysis/mucor_secretome/slurm/af3_fold.sh
#
# Card: no longer pinned to the a100. That pin was inherited, not measured - see
# scripts/af3_gpu_probe.sh, which checks jax device visibility and a matmul against a known
# answer on a given card. The AF3 container carries its own jax (0.9.1), so nothing proven about
# the torch env in scripts/setup_gpu_env.sh applies to it, and the probe is the evidence. The
# k80 (sm_37) and p100 (sm_60) nodes are too old and are in any case down/drained; everything
# Ampere or later is fair game. The job aborts below if the container cannot see a GPU, rather
# than silently running on CPU for hours.
#
# Sharding: work is split by `index % NTASK == TASK`, the same scheme as af3_msa.sh. Task 0
# builds the job JSONs and the rest wait for them, so one submission does everything. Because
# each prediction whose summary JSON already exists is skipped, shards never collide, a timeout
# costs only a resubmit, and the shard count can change between runs.
#
# Queue note: short_gpu caps at 2 h and is usually tens of jobs deep, while exfab is a single
# idle node (2x ada6000, 256 CPUs) with a 30-day limit. The unknowns set is ~426 predictions and
# does not fit in 2 h unsharded, so either shard it on short_gpu or send it to exfab with a
# longer -t.
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1 -n 1 -c 8
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -J af3_fold
#SBATCH -o logs/af3_fold_%j.out
#SBATCH -e logs/af3_fold_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
SET=${SET:-chitin}
OUT=$A/out/af3/$SET
AF3=/opt/linux/rocky/8.x/x86_64/pkgs/alphafold/3.0.2/alphafold
module load alphafold/3.0.2
mkdir -p $OUT/pairs $OUT/pred

TASK=${SLURM_ARRAY_TASK_ID:-0}
NTASK=${SLURM_ARRAY_TASK_COUNT:-1}

nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader || true
# Fail fast rather than spend the whole allocation on CPU. AF3 falls back silently otherwise,
# and a 30-day exfab job doing that is the expensive version of this mistake.
singularity exec --nv "$ALPHAFOLD_SING" python3 -c "
import sys, jax
d = jax.devices()
print('jax', jax.__version__, d)
sys.exit(0 if any(x.platform == 'gpu' for x in d) else 'ABORT: no GPU visible to jax')
" || exit 1

# Task 0 builds the job JSONs; the others wait for them, so one submission does both.
if [ "$TASK" = "0" ]; then
  if [ "$SET" = "chitin" ]; then
    $PY $A/bin/af3_input.py pairs --fasta $OUT/panel.faa --a $OUT/fungal.txt --b $OUT/host.txt \
        --msa-dir $OUT/msas --out $OUT/pairs
    $PY $A/bin/af3_input.py ligand --fasta $OUT/panel.faa --ids $OUT/ligand_ids.txt \
        --msa-dir $OUT/msas --nag 4,7 --out $OUT/pairs
  else
    # the unknowns are monomers: fold each with the MSA its stage-1 run produced. --msa-dir is
    # what makes that true; without it these are bare sequences and AF3 dies in featurising
    # under --norun_data_pipeline (job 32069699: 0 of 426 succeeded).
    $PY $A/bin/af3_input.py single --fasta $OUT/panel.faa --msa-dir $OUT/msas --out $OUT/pairs
  fi
  touch $OUT/pairs/.ready
else
  for _ in $(seq 1 60); do [ -f $OUT/pairs/.ready ] && break; sleep 10; done
fi

mapfile -t JOBS < <(ls $OUT/pairs/*.json | sort)
echo "task $TASK/$NTASK of ${#JOBS[@]} jobs"

for i in "${!JOBS[@]}"; do
  [ $(( i % NTASK )) -eq "$TASK" ] || continue
  j="${JOBS[$i]}"
  n=$(basename "$j" .json)
  if [ -s "$OUT/pred/${n,,}/${n,,}_summary_confidences.json" ]; then continue; fi
  echo "=== $n $(date +%T)"
  $AF3 --json_path="$j" --output_dir="$OUT/pred" --db_dir="$ALPHAFOLD_DB" \
       --model_dir="$ALPHAFOLD_DB/model" --norun_data_pipeline --force_output_dir \
       || echo "FAILED $n" >&2
done

# The report needs every prediction, so only the shard that finds the set complete writes it.
# Whichever task finishes last does that; if several finish together they write the same file
# from the same inputs, which is harmless. A set left incomplete by a timeout writes nothing and
# says so, instead of producing a report that silently covers part of the panel.
missing=0
for j in "${JOBS[@]}"; do
  n=$(basename "$j" .json)
  [ -s "$OUT/pred/${n,,}/${n,,}_summary_confidences.json" ] || missing=$(( missing + 1 ))
done
if [ "$missing" -gt 0 ]; then
  echo "task $TASK done; $missing of ${#JOBS[@]} predictions still missing - no report written"
  exit 0
fi
$PY $A/bin/af3_report.py --dir $OUT/pred --controls $OUT/controls.tsv --out $OUT/af3_scores.tsv
echo "=== best 20 ==="
head -21 $OUT/af3_scores.tsv | cut -f1,2,3,6,7,8 | column -t -s $'\t'
echo "done $(date +%T)"
