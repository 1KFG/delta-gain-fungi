#!/bin/bash
# How long does ESMFold take per protein on CPU? A timing probe, not an analysis.
#   sbatch analysis/mucor_secretome/slurm/esmfold_cpu_timing.sh
# Every a100 is allocated and the GPU job may wait days, while CPU jobs here start in seconds.
# ESMFold runs on CPU without any architecture constraint - the question is only whether 640
# proteins are hours or weeks of it. Folds 5 of the real candidates and reports seconds each.
# Writes to $SCRATCH only; the models are thrown away.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 1:00:00
#SBATCH -J esm_cpu_t
#SBATCH -o logs/esm_cpu_timing_%j.out
#SBATCH -e logs/esm_cpu_timing_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
W="${SCRATCH:?}/esm_cpu_t"; mkdir -p $W
export OMP_NUM_THREADS=${SLURM_CPUS_ON_NODE:-32}
export CONDA_NO_PLUGINS=true
source /etc/profile.d/modules.sh 2>/dev/null || true
module load miniconda3
eval "$(conda shell.bash hook)"
module load esm/2.0.0
F=$A/out/effectors/to_fold.faa
[ -s "$F" ] || { echo "no $F" >&2; exit 1; }
python3 $A/bin/esmfold_batch.py --fasta $F --out $W/folds --device cpu --limit 5 --chunk 64
echo "--- seconds per protein ---"
cat $W/folds/plddt.tsv
echo "done $(date +%T)"
