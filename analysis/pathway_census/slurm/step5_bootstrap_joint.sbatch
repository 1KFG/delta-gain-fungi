#!/bin/bash
# Bootstrap replicates for the joint tree (run after step5_prepare_joint.sbatch). One array task makes PC_BS_PER_TASK replicates
# of the trimmed alignment ref_set.trim.afa with raxml-ng --bootstrap (LG+G4, as the ML tree) and writes
# $BD/bs_<task>.raxml.bootstraps. Submit with the array range and CPU count you want, for example:
#   sbatch --export=ALL -p intel -c 1 --array=1-12 slurm/step5_bootstrap_joint.sbatch
# Required: PC_OUT (the run folder), PC_BS_PER_TASK. Optional: PC_BS_SEED_BASE (1000), PC_BS_DIR (bootstrap; use another folder per search mode),
#   PC_BS_EXTRA (extra raxml-ng options, for example '--opt-topology simplified --stop-rule kh-mult' = the --fast preset; changes the search, label the result). raxml-ng threads = SLURM_CPUS_PER_TASK.
# No BASH_SOURCE: all paths come from PC_OUT.
#SBATCH -p intel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 48:00:00
#SBATCH -J step5_bs_joint
#SBATCH -o step5_bs_joint_%A_%a.out
set -euo pipefail
: "${PC_OUT:?set PC_OUT}" "${PC_BS_PER_TASK:?set PC_BS_PER_TASK}"
TASK=${SLURM_ARRAY_TASK_ID:?run as an array job}
export PATH="$HOME/.pixi/bin:$PATH"
PX="pixi exec -c conda-forge -c bioconda -s raxml-ng==2.0.4 --"
BD=$PC_OUT/${PC_BS_DIR:-bootstrap}
mkdir -p "$BD"
SEED=$(( ${PC_BS_SEED_BASE:-1000} + TASK ))
T0=$(date +%s)
$PX raxml-ng --bootstrap --msa $PC_OUT/ref_set.trim.afa --model LG+G4 --prefix $BD/bs_$TASK \
  --bs-trees $PC_BS_PER_TASK ${PC_BS_EXTRA:-} --seed $SEED --threads ${SLURM_CPUS_PER_TASK:-1} --redo > $BD/bs_$TASK.stdout 2>&1
T1=$(date +%s)
N=$(wc -l < $BD/bs_$TASK.raxml.bootstraps)
echo "task $TASK seed $SEED replicates $N seconds $((T1 - T0)) threads ${SLURM_CPUS_PER_TASK:-1}" | tee $BD/bs_$TASK.time
rm -f $BD/bs_$TASK.raxml.rba
