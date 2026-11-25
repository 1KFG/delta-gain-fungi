#!/bin/bash
# Support values for the joint tree: all bootstrap trees onto the best ML tree, FBP and TBE (raxml-ng --support).
# Required: PC_OUT. Reads $PC_OUT/bootstrap/bs_*.raxml.bootstraps, $PC_OUT/ref.raxml.bestTree.
# Writes bootstrap/all.bootstraps, support.raxml.support (FBP), support.raxml.supportTBE (TBE).
#SBATCH -p intel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -J step5_sup_joint
#SBATCH -o step5_sup_joint_%j.out
set -euo pipefail
: "${PC_OUT:?set PC_OUT}"
export PATH="$HOME/.pixi/bin:$PATH"
PX="pixi exec -c conda-forge -c bioconda -s raxml-ng==2.0.4 --"
cat $PC_OUT/bootstrap/bs_*.raxml.bootstraps > $PC_OUT/bootstrap/all.bootstraps
N=$(wc -l < $PC_OUT/bootstrap/all.bootstraps)
echo "bootstrap trees: $N"
$PX raxml-ng --support --tree $PC_OUT/ref.raxml.bestTree --bs-trees $PC_OUT/bootstrap/all.bootstraps --bs-metric fbp,tbe \
  --prefix $PC_OUT/support --threads ${SLURM_CPUS_PER_TASK:-4} --redo > $PC_OUT/support.stdout 2>&1
ls $PC_OUT/support.raxml.support $PC_OUT/support.raxml.supportTBE
