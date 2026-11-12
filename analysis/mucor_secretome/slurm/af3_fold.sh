#!/bin/bash
# Stage 2 of the AF3 screen: inference only, reusing the stage-1 MSAs. GPU, resumable.
#   sbatch analysis/mucor_secretome/slurm/af3_fold.sh              # chitin pairs + ligands
#   sbatch --export=ALL,SET=unknowns analysis/mucor_secretome/slurm/af3_fold.sh
# a100 for the same reason ESMFold needs one: AF3's triton flash attention wants Ampere or
# later, and the a100 is the card this stack is known to work on here. Every prediction whose
# output already exists is skipped, so a 2 h limit on a busy queue costs only a resubmit.
#SBATCH -p short_gpu
#SBATCH --gres=gpu:a100:1
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

if [ "$SET" = "chitin" ]; then
  $PY $A/bin/af3_input.py pairs --fasta $OUT/panel.faa --a $OUT/fungal.txt --b $OUT/host.txt \
      --msa-dir $OUT/msas --out $OUT/pairs
  $PY $A/bin/af3_input.py ligand --fasta $OUT/panel.faa --ids $OUT/ligand_ids.txt \
      --msa-dir $OUT/msas --nag 4,7 --out $OUT/pairs
else
  # the unknowns are monomers: fold each with the MSA its stage-1 run produced
  $PY $A/bin/af3_input.py single --fasta $OUT/panel.faa --out $OUT/pairs
fi
echo "jobs: $(ls $OUT/pairs/*.json | wc -l)"

for j in $OUT/pairs/*.json; do
  n=$(basename "$j" .json)
  if [ -s "$OUT/pred/${n,,}/${n,,}_summary_confidences.json" ]; then continue; fi
  echo "=== $n $(date +%T)"
  $AF3 --json_path="$j" --output_dir="$OUT/pred" --db_dir="$ALPHAFOLD_DB" \
       --model_dir="$ALPHAFOLD_DB/model" --norun_data_pipeline --force_output_dir \
       || echo "FAILED $n" >&2
done
$PY $A/bin/af3_report.py --dir $OUT/pred --controls $OUT/controls.tsv --out $OUT/af3_scores.tsv
echo "=== best 20 ==="
head -21 $OUT/af3_scores.tsv | cut -f1,2,3,6,7,8 | column -t -s $'\t'
echo "done $(date +%T)"
