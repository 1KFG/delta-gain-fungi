#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J upo_refbs
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/refbs_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/refbs_%j.err
set -euo pipefail

# Bootstrap support for the PF01328 reference tree.
#
# The question this answers: is MphUPO -- the sole Clade IV reference -- really nested INSIDE
# Clade V? The single ML search says yes. If bootstrap support agrees, then "Clade IV" is not a
# separate clade in this reference framework, and the G1 Finding 3 headline needs revising.
# Run on a compute node: the login node was too loaded and the first attempt stalled.

P=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/placement
module load raxml-ng/1.2.0
cd "$P"
rm -f refbs.raxml.*
raxml-ng --all --msa ref_aln.filtered.fa --model LG+G4 --prefix refbs \
         --threads "${SLURM_CPUS_PER_TASK:-8}" --seed 1337 --bs-trees 400 --force perf_threads
echo "bootstraps done $(date +%T)"
