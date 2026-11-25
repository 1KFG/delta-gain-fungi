#!/bin/bash
# Re-apply the clade rule with another decoy scope to a SAVED joint placement (no new tree, no new epa-ng run).
# Usage (all paths from environment, no BASH_SOURCE): PC_RUN=<run folder> PC_NEWDIR=<subfolder name> PC_ROOT=... PC_REFS_STEPS=... PC_DECOYS_FA=...
#   PC_MEMBERS_FA=... PC_STEP_DECOYS="STEP=fasta ..." PC_STEPS="..." PC_STEP_CAND_DIR=... PC_SCORES=... PC_GENOMES=... PC_MIN_SUPPORT=0.5 [PC_DECOY_SCOPE=per-step]
# Needs in PC_RUN: ref.raxml.bestTree, rep_ids.txt, epa_result.jplace.gz, clusters.tsv.zst. Writes PC_RUN/PC_NEWDIR/.
set -euo pipefail
: "${PC_RUN:?}" "${PC_NEWDIR:?}" "${PC_ROOT:?}" "${PC_REFS_STEPS:?}" "${PC_DECOYS_FA:?}" "${PC_STEPS:?}" "${PC_STEP_CAND_DIR:?}" "${PC_SCORES:?}" "${PC_GENOMES:?}"
PY=/usr/bin/python3.12; UPO=$PC_ROOT/../enzymes_UPO/upo_place/bin/upo_place.py; J=$PC_ROOT/bin/step5_joint.py
D=$PC_RUN/$PC_NEWDIR; mkdir -p $D
ln -sf $PC_RUN/ref.raxml.bestTree $PC_RUN/rep_ids.txt $D/
rm -f $D/clusters.tsv.zst; cp $PC_RUN/clusters.tsv.zst $D/   # zstd refuses a symbolic link
SCOPE="--decoy-scope ${PC_DECOY_SCOPE:-per-step}"
[ "${PC_DECOY_SCOPE:-per-step}" = per-step ] && SCOPE="$SCOPE --step-decoys $PC_STEP_DECOYS"
set +e
$PY $J clades $SCOPE --tree $PC_RUN/ref.raxml.bestTree --refs $PC_REFS_STEPS --decoys $PC_DECOYS_FA ${PC_MEMBERS_FA:+--members-fa $PC_MEMBERS_FA} \
  --reps $PC_RUN/rep_ids.txt --joint-id pks_joint --out-dir $D
rc=$?; set -e; [ $rc = 0 ] || [ $rc = 3 ] || exit $rc
JP=$D/epa_result.jplace; zcat $PC_RUN/epa_result.jplace.gz > $JP
for S in $PC_STEPS; do
  [ -s $D/clade_members.$S.tsv ] || { echo "no clade file for $S" | tee -a $D/place_skipped.txt; continue; }
  $PY $UPO assign $JP --clades $D/clade_members.$S.tsv --min-lwr ${PC_MIN_LWR:-0.5} --max-pendant ${PC_MAX_PENDANT:-1.0} --out $D/assign.$S.tsv
  $PY $J filter-placements --assign $D/assign.$S.tsv --step $S --cand-dir $PC_STEP_CAND_DIR --out-plc $D/placements.$S.tsv \
    --out-detail $D/placements_detail.$S.tsv --out-cand $D/candidates_step.$S.tsv | tee -a $D/filter_placements.jsonl
  $PY $PC_ROOT/bin/assign_steps.py --pathway $PC_ROOT/pathways/${PC_PATHWAY:-melanin_dhn}.yaml --candidates $D/candidates_step.$S.tsv \
    --placements $D/placements.$S.tsv --min-support ${PC_MIN_SUPPORT:-0.5} --out $D/assigned.$S.tsv --provisional
done
rm -f $JP
$PY $J report --out-dir $D --scores $PC_SCORES --genomes $PC_GENOMES
echo "reapply done $(date +%T)"
