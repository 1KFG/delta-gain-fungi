#!/bin/bash
# Step 5, part 1b, for a step with ONE reference: replace the clade definition made by step5_prepare.sbatch.
# Run after step5_prepare.sbatch, before step5_place.sbatch. See bin/step5_single_ref_clades.py.
# Required: PC_STEP. Optional: PC_PATHWAY, PC_ROOT, PC_OUT, PC_DECOYS_FA, PC_OUTGROUP_FA, PC_MEMBERS_FA, PC_OUTGROUP_ROLE.
#SBATCH -p short
#SBATCH --constraint="ryzen|cascade|broadwell"
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0:20:00
#SBATCH -J step5_sr_clades
#SBATCH -o step5_sr_clades_%j.out
set -euo pipefail
: "${PC_STEP:?set PC_STEP}"
PC_PATHWAY=${PC_PATHWAY:-carotenoid}
PC_ROOT=${PC_ROOT:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/pathway_census}
OUT=${PC_OUT:-$PC_ROOT/results/step5/$PC_PATHWAY/$PC_STEP}
REFS=$PC_ROOT/refs/$PC_PATHWAY/$PC_STEP.faa
DECOYS=${PC_DECOYS_FA:-$PC_ROOT/refs/$PC_PATHWAY/decoys/$PC_STEP.faa}
/usr/bin/python3.12 $PC_ROOT/bin/step5_single_ref_clades.py --tree $OUT/ref.raxml.bestTree --refs $REFS --decoys $DECOYS \
  ${PC_OUTGROUP_FA:+--outgroup $PC_OUTGROUP_FA} ${PC_MEMBERS_FA:+--members-fa $PC_MEMBERS_FA} \
  --outgroup-role ${PC_OUTGROUP_ROLE:-tree_outgroup} --reps $OUT/rep_ids.txt --step $PC_STEP --out-dir $OUT
