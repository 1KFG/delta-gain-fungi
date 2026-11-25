#!/bin/bash
# Joint SDR tree (t4hn_reductase + t3hn_reductase): prepare on the joint candidate folder, derive two step folders, place each.
# usage: submit_sdr_joint.sh [min_seq_id] [partition]  (run from analysis/pathway_census)
set -euo pipefail
ID=${1:-0.5}; PART=${2:-short}
T=$PWD/results/step5/other_inputs/tree_inputs/sdr_joint
export PC_ROOT=$PWD PC_PATHWAY=melanin_dhn PC_STEP=t4hn_reductase PC_CAND_DIR=$PWD/results/search/sdr_joint
export PC_OUT=$PWD/results/step5/melanin_dhn/sdr_joint PC_DECOYS_FA=$T/decoys.faa PC_OUTGROUP_ROLE=tree_outgroup
export PC_RAXML_THREADS=4 PC_RAXML_STARTS='pars{2}' PC_MIN_SUPPORT=0.5 PC_MIN_SEQ_ID=$ID
export PC_MEMBERS_FA=$T/members.faa PC_OUTGROUP_FA=$T/outgroups.faa
mkdir -p $PC_OUT
# the prepare job may exit 3 (default rule: references not a clade free of decoys/outgroups); use afterany so the derive step still runs
DEPARG=""; [ -n "${DEP:-}" ] && DEPARG="--dependency=afterany:$DEP"
j1=$(sbatch --parsable --export=ALL -p $PART -c 15 -D $PC_OUT $DEPARG slurm/step5_prepare_ident.sbatch)
unset PC_MEMBERS_FA
export PC_CAND_DIR=$PWD/results/search/melanin_myco_full
prevp=""; jd=""; for S in t4hn_reductase t3hn_reductase; do
  d=$(sbatch --parsable --export=ALL -p $PART -c 1 --mem=4G -t 0:30:00 -D $PC_OUT --dependency=afterok:$j1 -J sdr_derive_$S --wrap="bash $PWD/results/step5/other_inputs/derive_sdr_step.sh $S")
  mkdir -p $PWD/results/step5/melanin_dhn/$S
  p=$(PC_STEP=$S PC_OUT=$PWD/results/step5/melanin_dhn/$S sbatch --parsable --export=ALL -p $PART -c 15 -D $PWD/results/step5/melanin_dhn/$S --dependency=afterok:$d${prevp:+,afterany:$prevp} slurm/step5_place.sbatch)
  echo "sdr $S derive=$d place=$p"
  jd="$jd $d $p"; prevp=$p
done
echo "sdr_joint prepare=$j1 $jd id=$ID part=$PART $(date +%F_%T)" | tee -a results/step5/other_steps_jobs.txt
echo $p > results/step5/.last_place_sdr_joint
