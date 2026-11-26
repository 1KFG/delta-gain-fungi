#!/bin/bash
# Rarefaction check (PROPOSED, approximate) for the joint tree. Run after step5_place_joint.sbatch (needs query_match_cols.fa.zst).
# For each fraction F in PC_RARE_FRACS and seed S in PC_RARE_SEEDS: keep all labelled leaves and the fraction F of the
# representatives, prune the ML tree and the alignment (branch lengths of removed leaves are summed into the neighbour; model
# parameters and the other branch lengths are NOT re-optimised), redefine the clades by the same rule, place the same query
# alignment with epa-ng, assign per step, and compare the calls with the full tree (bin/step5_bootstrap.py compare).
# Required: PC_OUT, PC_ROOT, PC_STEP_CAND_DIR_DUMMY not used. Needs PC_REFS_STEPS, PC_DECOYS_FA, PC_MEMBERS_FA, PC_STEPS (as env.sh).
# Optional: PC_RARE_FRACS ("0.5 0.75"), PC_RARE_SEEDS ("1 2 3"), PC_MIN_LWR (0.5), PC_MAX_PENDANT (1.0).
#SBATCH -p intel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -J step5_rare_joint
#SBATCH -o step5_rare_joint_%j.out
set -euo pipefail
: "${PC_OUT:?}" "${PC_ROOT:?}" "${PC_REFS_STEPS:?}" "${PC_DECOYS_FA:?}" "${PC_STEPS:?}"
export PATH="$HOME/.pixi/bin:$PATH"
module load zstd 2>/dev/null || true
PX="pixi exec -c conda-forge -c bioconda -s epa-ng==0.3.8 -s raxml-ng==2.0.4 --"
PY=/usr/bin/python3.12
UPO=$PC_ROOT/../enzymes_UPO/upo_place/bin/upo_place.py
B=$PC_ROOT/bin/step5_bootstrap.py
J=$PC_ROOT/bin/step5_joint.py
OUT=$PC_OUT
R=$OUT/rarefaction; mkdir -p $R
W="${SCRATCH:?}/rare_joint"; mkdir -p $W
if [ "${PC_DECOY_SCOPE:-per-step}" = union ]; then SCOPE="--decoy-scope union"; else
  : "${PC_STEP_DECOYS:?set PC_STEP_DECOYS or PC_DECOY_SCOPE=union}"; SCOPE="--decoy-scope per-step --step-decoys $PC_STEP_DECOYS"; fi
zstd -dc $OUT/query_match_cols.fa.zst > $W/q.fa
for F in ${PC_RARE_FRACS:-0.5 0.75}; do
  for S in ${PC_RARE_SEEDS:-1 2 3}; do
    D=$R/f${F}_s${S}; mkdir -p $D
    $PY $B subset --tree $OUT/ref.raxml.bestTree --msa $OUT/ref_set.trim.afa --rep-ids $OUT/rep_ids.txt --fraction $F --seed $S \
      --out-tree $D/sub.nwk --out-msa $D/sub.afa --out-reps $D/rep_ids.txt > $D/subset.json
    $PY $J clades $SCOPE --tree $D/sub.nwk --refs $PC_REFS_STEPS --decoys $PC_DECOYS_FA ${PC_MEMBERS_FA:+--members-fa $PC_MEMBERS_FA} \
      --reps $D/rep_ids.txt --joint-id rare --out-dir $D || echo "clades rc=$? (a step may have no target)"
    mkdir -p $W/e_${F}_$S
    $PX epa-ng --ref-msa $D/sub.afa --tree $D/sub.nwk --query $W/q.fa --model $OUT/ref.raxml.bestModel -w $W/e_${F}_$S --redo -T ${SLURM_CPUS_PER_TASK:-4} > $D/epa-ng.log 2>&1
    for ST in $PC_STEPS; do
      [ -s $D/clade_members.$ST.tsv ] || { echo "$ST: no target on this subset tree" >> $D/no_target.txt; continue; }
      $PY $UPO assign $W/e_${F}_$S/epa_result.jplace --clades $D/clade_members.$ST.tsv --min-lwr ${PC_MIN_LWR:-0.5} --max-pendant ${PC_MAX_PENDANT:-1.0} --out $D/assign.$ST.tsv
      [ -s $OUT/assign.$ST.tsv ] && $PY $B compare --old-assign $OUT/assign.$ST.tsv --new-assign $D/assign.$ST.tsv \
        --step-proteins $OUT/candidates_step.$ST.tsv --out-tsv $D/compare_vs_full.$ST.tsv --out-json $D/compare_vs_full.$ST.json
    done
    rm -rf $W/e_${F}_$S
  done
done
echo "rarefaction done $(date +%T)"
