#!/bin/bash
# ML tree of the joint alignment plus clades. Same raxml-ng call as the tree stage of step5_prepare_joint.sbatch, as a separate job
# so that the CPUs of the tree job and of the bootstrap array can be set apart (run after the alignment stage of step5_prepare_joint).
# Required: PC_OUT, PC_ROOT, PC_STEP, PC_REFS_STEPS, PC_DECOYS_FA, PC_MEMBERS_FA. Optional: PC_RAXML_STARTS (pars{2}), PC_RAXML_WORKERS (2).
#SBATCH -p intel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -J step5_tree_joint
#SBATCH -o step5_tree_joint_%j.out
set -euo pipefail
: "${PC_OUT:?set PC_OUT}" "${PC_ROOT:?set PC_ROOT}" "${PC_STEP:?}" "${PC_REFS_STEPS:?}" "${PC_DECOYS_FA:?}"
export PATH="$HOME/.pixi/bin:$PATH"
PX="pixi exec -c conda-forge -c bioconda -s epa-ng==0.3.8 -s raxml-ng==2.0.4 --"
OUT=$PC_OUT
$PX raxml-ng --msa $OUT/ref_set.trim.afa --model LG+G4 --prefix $OUT/ref --seed 1 --tree ${PC_RAXML_STARTS:-pars{2}} --threads ${SLURM_CPUS_PER_TASK:-4} \
  --workers ${PC_RAXML_WORKERS:-2} --redo > $OUT/raxml.log 2>&1
rm -f $OUT/ref.raxml.startTree $OUT/ref.raxml.rba $OUT/ref.raxml.reduced.phy
# decoy scope: per step (default; PC_STEP_DECOYS="STEP=decoys.faa ...") or the first run's union (PC_DECOY_SCOPE=union)
if [ "${PC_DECOY_SCOPE:-per-step}" = union ]; then SCOPE="--decoy-scope union"; else
  : "${PC_STEP_DECOYS:?set PC_STEP_DECOYS (STEP=FASTA per step) or PC_DECOY_SCOPE=union}"; SCOPE="--decoy-scope per-step --step-decoys $PC_STEP_DECOYS"; fi
set +e
/usr/bin/python3.12 $PC_ROOT/bin/step5_joint.py clades $SCOPE --tree $OUT/ref.raxml.bestTree --refs $PC_REFS_STEPS --decoys $PC_DECOYS_FA \
  ${PC_MEMBERS_FA:+--members-fa $PC_MEMBERS_FA} --reps $OUT/rep_ids.txt --joint-id $PC_STEP --out-dir $OUT
rc=$?
echo "step5_tree_joint done rc=$rc $(date +%T)"
exit $rc
