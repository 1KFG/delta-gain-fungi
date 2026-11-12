#!/bin/bash
# SignalP 6 on the three non-BFD proteomes, on a GPU, the way BFD runs it.
#   sbatch analysis/mucor_secretome/slurm/signalp_gpu.sh
# Mirrors Fungi_BFD/nextflow/modules/BFD/SIGNALP/main.nf and the `signalp` label in
# conf/profile_BFD.config: short_gpu, one GPU, the signalp/6-gpu module, and the GFF3 output,
# which is what merge_signalp.py builds signalp.signal_peptide.parquet from.
#
# The module itself only pulls in miniconda3 and then emits a `conda activate` for the shell to
# evaluate, so the hook below is required - without it the module load fails with
# "Run 'conda init' before 'conda activate'", which is what killed job 32013761.
# Reads shared folders read-only; writes mucor_secretome/out/panels/.
# The BFD profile's `signalp` label uses short_gpu; exfab is this group's GPU queue (one node,
# 2 x ada6000, 256 cores), so jobs here go to exfab instead.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -N 1 -n 1 -c 2
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -J signalp_panels
#SBATCH -o logs/signalp_panels_%j.out
#SBATCH -e logs/signalp_panels_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
OUT=$A/out/panels; mkdir -p $OUT
FAA=$OUT/faa
W="${SCRATCH:?}/signalp_panels"; mkdir -p $W

source /etc/profile.d/modules.sh 2>/dev/null || true
export CONDA_NO_PLUGINS=true
module load miniconda3
eval "$(conda shell.bash hook)"
module load signalp/6-gpu

bash $A/bin/panel_proteomes.sh $FAA
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

for LT in JEM1MUCOR JEM2MUCOR MPTAK; do
  echo "=== $LT $(date +%T)"
  rm -rf $W/sp_$LT
  signalp6 -od $W/sp_$LT -org euk --mode fast -format none \
      -fasta $FAA/$LT.faa --write_procs 2 -bs 100
  $PY $A/bin/annotate_to_tables.py --locustag $LT --out $OUT \
      --signalp-gff3 $W/sp_$LT/output.gff3
  cp $W/sp_$LT/prediction_results.txt $OUT/${LT}.signalp_results.txt
done
wc -l $OUT/signalp.tsv
echo "done $(date +%T)"
