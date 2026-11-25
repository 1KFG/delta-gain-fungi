#!/bin/bash
# Support values for the joint tree: all bootstrap trees onto the best ML tree, FBP and TBE (raxml-ng --support).
# Required: PC_OUT. Optional: PC_BS_DIR (bootstrap folder, default bootstrap; the output prefix gets _<dir>). Reads $BD/bs_*.raxml.bootstraps, $PC_OUT/ref.raxml.bestTree.
# Writes bootstrap/all.bootstraps, support.raxml.supportFBP (FBP), support.raxml.supportTBE (TBE).
#SBATCH -p intel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -J step5_sup_joint
#SBATCH -o step5_sup_joint_%j.out
set -euo pipefail
: "${PC_OUT:?set PC_OUT}"
BD=$PC_OUT/${PC_BS_DIR:-bootstrap}
export PATH="$HOME/.pixi/bin:$PATH"
PX="pixi exec -c conda-forge -c bioconda -s raxml-ng==2.0.4 --"
cat $BD/bs_*.raxml.bootstraps > $BD/all.bootstraps
N=$(wc -l < $BD/all.bootstraps)
echo "bootstrap trees: $N"
$PX raxml-ng --support --tree $PC_OUT/ref.raxml.bestTree --bs-trees $BD/all.bootstraps --bs-metric fbp,tbe \
  --prefix $PC_OUT/support${PC_BS_DIR:+_$PC_BS_DIR} --threads ${SLURM_CPUS_PER_TASK:-4} --redo > $PC_OUT/support${PC_BS_DIR:+_$PC_BS_DIR}.stdout 2>&1
ls $PC_OUT/support${PC_BS_DIR:+_$PC_BS_DIR}.raxml.supportFBP $PC_OUT/support${PC_BS_DIR:+_$PC_BS_DIR}.raxml.supportTBE
