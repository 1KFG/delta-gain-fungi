#!/bin/bash
# Build the placement reference from the in-house tree: an HMM whose match states are exactly the
# 224 reference columns (hmmbuild --symfrac 0), and raxml-ng model parameters + branch lengths
# for the mfp topology. Output: data/refset_zamocky2022/placement_ref/ (committed).
#   sbatch analysis/enzymes_UPO/upo_place/slurm/prepare_ref.sh
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH -J upo_place_ref
#SBATCH -o logs/place_ref_%j.out
#SBATCH -e logs/place_ref_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source analysis/enzymes_UPO/upo_place/slurm/env.sh
mkdir -p $PREF
hmmbuild --symfrac 0 --informat afa -n upo_ref $PREF/ref.hmm $TREEDIR/aln.trim.afa > $PREF/hmmbuild.log
LEN=$(awk '$1=="LENG"{print $2}' $PREF/ref.hmm)
[ "$LEN" = "$NCOL" ] || { echo "ref.hmm has $LEN match states, expected $NCOL" >&2; exit 1; }
# raxml-ng rejects support labels on some inputs; strip them
$PY - $TREEDIR/mfp.treefile $PREF/mfp.nolabels.nwk <<'PYEOF'
import re, sys
s = open(sys.argv[1]).read()
open(sys.argv[2], "w").write(re.sub(r"\)[^:;,()]+", ")", s))
PYEOF
$PX raxml-ng --evaluate --msa $TREEDIR/aln.trim.afa --tree $PREF/mfp.nolabels.nwk \
    --model "$MODEL" --prefix $PREF/ref --threads ${SLURM_CPUS_PER_TASK:-4} --redo --seed 20261108
rm -f $PREF/ref.raxml.startTree $PREF/ref.raxml.rba $PREF/ref.raxml.reduced.phy
echo "done $(date +%T)"
