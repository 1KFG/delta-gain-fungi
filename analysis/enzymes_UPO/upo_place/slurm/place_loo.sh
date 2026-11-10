#!/bin/bash
# Leave-one-out validation of placement clade calls over all 187 reference leaves. Each leaf is
# pruned from the tree and MSA, the HMM is rebuilt without it, its raw sequence is re-aligned with
# hmmalign (the same path a BFD protein takes), placed with epa-ng, and assigned.
#   sbatch analysis/enzymes_UPO/upo_place/slurm/place_loo.sh      (after prepare_ref.sh)
# Output: data/refset_zamocky2022/placement_ref/loo.tsv (committed) + summary in the log.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 6:00:00
#SBATCH -J upo_place_loo
#SBATCH -o logs/place_loo_%j.out
#SBATCH -e logs/place_loo_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source analysis/enzymes_UPO/upo_place/slurm/env.sh
W="${SCRATCH:?}/place_loo"; mkdir -p $W
$PY -m unittest discover -s $UPO/upo_place/test
# warm the pixi environment once, so parallel workers do not race to solve it
$PX epa-ng --version > /dev/null

one() {
  set -euo pipefail
  leaf=$1; d=$W/$leaf; mkdir -p $d
  $PY $TOOL prune --tree $PREF/ref.raxml.bestTree --msa $TREEDIR/aln.trim.afa --leaf "$leaf" \
      --out-tree $d/tree.nwk --out-msa $d/ref.afa --out-query $d/self_aligned.fa
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
grep '^>' $TREEDIR/tree_input.faa | sed 's/^>//' \
  | xargs -P ${SLURM_CPUS_PER_TASK:-4} -I{} bash -c 'one "$@" || echo "FAILED $1" >&2' _ {}
$PY $TOOL loo $W/*/assign.tsv --clades $TREEDIR/mfp.clade_members.tsv --out $PREF/loo.tsv
echo "done $(date +%T)"
