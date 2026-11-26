#!/bin/bash
# Submit step 5 for one step: prepare -> multi_ref_clades -> place (afterok chain). Run from analysis/pathway_census.
# usage: [DEP=<job id>] submit_step.sh <pathway> <step> [min_seq_id] [partition]   (DEP: start after that job ends, to cap the CPUs in use)
# Env set here: PC_ROOT=$PWD PC_CAND_DIR (both parts read it) PC_OUT, inputs from results/step5/other_inputs/tree_inputs/<step>.
set -euo pipefail
P=$1; S=$2; ID=${3:-0.5}; PART=${4:-short}
T=$PWD/results/step5/other_inputs/tree_inputs/$S
export PC_ROOT=$PWD PC_PATHWAY=$P PC_STEP=$S PC_CAND_DIR=${PC_CAND_DIR:-$PWD/results/search/melanin_myco_full}
export PC_OUT=$PWD/results/step5/$P/$S PC_DECOYS_FA=$T/decoys.faa PC_OUTGROUP_ROLE=tree_outgroup
export PC_RAXML_THREADS=4 PC_RAXML_STARTS='pars{2}' PC_MIN_SUPPORT=0.5 PC_MIN_SEQ_ID=$ID
[ -s $T/members.faa ] && export PC_MEMBERS_FA=$T/members.faa || unset PC_MEMBERS_FA
[ -s $T/outgroups.faa ] && export PC_OUTGROUP_FA=$T/outgroups.faa || unset PC_OUTGROUP_FA
mkdir -p $PC_OUT
SB="sbatch --parsable --export=ALL -p $PART -c 15 -D $PC_OUT"
DEPARG=""; [ -n "${DEP:-}" ] && DEPARG="--dependency=afterany:$DEP"
j1=$($SB $DEPARG slurm/step5_prepare_ident.sbatch)
j2=$(sbatch --parsable --export=ALL -p $PART -D $PC_OUT --dependency=afterok:$j1 slurm/step5_multi_ref_clades.sbatch)
j3=$($SB --dependency=afterok:$j2 slurm/step5_place.sbatch)
echo "$S prepare=$j1 clades=$j2 place=$j3 id=$ID part=$PART $(date +%F_%T)" | tee -a results/step5/other_steps_jobs.txt
echo $j3 > results/step5/.last_place_$S
