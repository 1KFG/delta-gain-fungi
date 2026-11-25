#!/bin/bash
# Submit a step 5 v2 run: prepare -> multi_ref_clades -> place (afterok chain). Run from analysis/pathway_census (2026-11-24).
# usage: [PLACE=slurm/step5_place_parallel.sbatch] submit_v2.sh <pathway> <step> <input dir> <out folder name> <min_seq_id> <max_reps> <cpus> <time> [min_cov]
#   input dir   a folder of results/step5/other_inputs/tree_inputs/ (decoys.faa, members.faa, outgroups.faa, optional refs.faa)
#   refs.faa    when present it replaces refs/<pathway>/<step>.faa (PC_REFS_FA; refs/ is not changed)
# Env: PC_ROOT=$PWD, PC_CAND_DIR (read by prepare AND place; the same value), PC_OUT=results/step5/<pathway>/<out name>.
set -euo pipefail
P=$1; S=$2; T=$PWD/results/step5/other_inputs/tree_inputs/$3; NAME=$4; ID=$5; CAP=$6; CPUS=$7; TL=$8; COV=${9:-0.8}
PLACE=${PLACE:-slurm/step5_place.sbatch}
export PC_ROOT=$PWD PC_PATHWAY=$P PC_STEP=$S PC_CAND_DIR=$PWD/results/search/melanin_myco_full
export PC_OUT=$PWD/results/step5/$P/$NAME PC_DECOYS_FA=$T/decoys.faa PC_OUTGROUP_ROLE=tree_outgroup
export PC_RAXML_THREADS=$CPUS PC_RAXML_STARTS='pars{2}' PC_MIN_SUPPORT=0.5 PC_MIN_SEQ_ID=$ID PC_MIN_COV=$COV PC_MAX_REPS=$CAP
[ -s $T/members.faa ] && export PC_MEMBERS_FA=$T/members.faa || unset PC_MEMBERS_FA
[ -s $T/outgroups.faa ] && export PC_OUTGROUP_FA=$T/outgroups.faa || unset PC_OUTGROUP_FA
[ -s $T/refs.faa ] && export PC_REFS_FA=$T/refs.faa || unset PC_REFS_FA
mkdir -p $PC_OUT
SB="sbatch --parsable --export=ALL -p intel -D $PC_OUT"
j1=$($SB -c $CPUS --mem=32G -t $TL slurm/step5_prepare_ident.sbatch)
j2=$($SB -c 1 --mem=4G -t 0:20:00 --dependency=afterok:$j1 slurm/step5_multi_ref_clades.sbatch)
j3=$($SB -c $CPUS --mem=32G -t $TL --dependency=afterok:$j2 $PLACE)
echo "$NAME prepare=$j1 clades=$j2 place=$j3 id=$ID cap=$CAP cpus=$CPUS $(date +%F_%T)" | tee -a results/step5/other_steps_jobs.v2.txt
