#!/bin/bash
# Matched-identity validation of the 187-leaf placement reference (same design as the skylark
# session's job 31868181): per query and ceiling, drop the query and every reference above the
# ceiling, rebuild the HMM, re-align the raw query, place, assign.
#   sbatch analysis/enzymes_UPO/upo_place/slurm/place_matched.sh      (after prepare_ref.sh)
# Output: data/refset_zamocky2022/placement_ref/{matched_cases.tsv,MATCHED_IDENTITY.md}
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 8:00:00
#SBATCH -J upo_place_matched
#SBATCH -o logs/place_matched_%j.out
#SBATCH -e logs/place_matched_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source analysis/enzymes_UPO/upo_place/slurm/env.sh
W="${SCRATCH:?}/place_matched"; mkdir -p $W
MI=$UPO/upo_place/bin/matched_identity.py
$PY -W ignore -m unittest discover -s $UPO/upo_place/test
$PY $MI cases --msa $TREEDIR/aln.trim.afa --clades $TREEDIR/mfp.clade_members.tsv \
    --out $PREF/matched_cases.tsv
$PX epa-ng --version > /dev/null

one() {
  set -euo pipefail
  case_id=$1; leaf=$2; removed=$3; possible=$4
  [ "$possible" = yes ] || return 0
  d=$W/${case_id//@/__at__}; mkdir -p $d
  $PY $TOOL prune --tree $PREF/ref.raxml.bestTree --msa $TREEDIR/aln.trim.afa --leaf "$leaf" \
      --remove "$removed" --out-tree $d/tree.nwk --out-msa $d/ref.afa --out-query $d/self.fa
  hmmbuild --symfrac 0 --informat afa $d/ref.hmm $d/ref.afa > /dev/null
  awk -v L=">$leaf" '/^>/{p=($0==L)} p' $TREEDIR/tree_input.faa > $d/raw.fa
  hmmalign --outformat afa $d/ref.hmm $d/raw.fa > $d/aln.afa
  $PY $TOOL match-cols $d/aln.afa --length $NCOL --out $d/q.fa 2>/dev/null
  $PX epa-ng --ref-msa $d/ref.afa --tree $d/tree.nwk --query $d/q.fa \
      --model $PREF/ref.raxml.bestModel -w $d --redo -T 1 > $d/epa.log 2>&1
  $PY $TOOL assign $d/epa_result.jplace --clades $TREEDIR/mfp.clade_members.tsv \
      --out $d/assign.tsv 2>/dev/null
}
export -f one; export W PY TOOL PREF TREEDIR NCOL PX
# columns: case(1) query(2) ... possible(6) ... removed(10). `removed` goes LAST: it is empty for
# plain leave-one-out, and read with a tab IFS collapses empty fields in the middle.
tail -n +2 $PREF/matched_cases.tsv | awk -F'\t' '{print $1"\t"$2"\t"$6"\t"$10}' \
  | xargs -P ${SLURM_CPUS_PER_TASK:-4} -d '\n' -I{} bash -c \
      'IFS=$'"'"'\t'"'"' read -r c q p r <<< "$1"; one "$c" "$q" "$r" "$p" || echo "FAILED $c" >&2' _ {}
$PY $MI summarize --cases $PREF/matched_cases.tsv --assign-dir $W --out $PREF/MATCHED_IDENTITY.md
echo "done $(date +%T)"
