#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 8 --mem 8G -t 0:59:00
#SBATCH -J upo_loo -o upo_loo.%j.log
# Leave-one-out test: for every reference used to train a family HMM (except the anchor AaeUPO),
# rebuild all models without it, then classify it.
# Usage: sbatch tests/run_loo.sh <upo_id_dir>     (run from anywhere; pass the tool folder)
set -euo pipefail
UPO_DIR=${1:?give the upo_id folder}
module load hmmer/3.4 mafft/7.505
PY=/usr/bin/python3.12
WORK=${SCRATCH:-${TMPDIR:?}}/upo_loo
mkdir -p "$WORK"
cd "$UPO_DIR"
OUT=tests/loo
mkdir -p $OUT
names=$(awk -F'\t' 'NR>1 && ($2=="FamII"||$2=="FamI_UPO"||$2=="CPO") && $1!="AaeUPO"{print $1}' ref/train_sets.tsv)
run_one() {
  n=$1
  $PY upo_id/upo_id.py build --refs ref/refs_all.fa --train ref/train_sets.tsv --manifest ref/manifest.tsv \
      --pfam-hmm ref/PF01328.hmm --out "$WORK/m_$n" --exclude "$n" --cpu 2 2> "$WORK/$n.build.log"
  $PY - "$n" ref/refs_all.fa > "$WORK/$n.fa" <<'EOF'
import sys
n, fa = sys.argv[1], sys.argv[2]
keep = False
for line in open(fa):
    if line.startswith('>'):
        keep = line[1:].split()[0] == n
    if keep:
        sys.stdout.write(line)
EOF
  $PY upo_id/upo_id.py classify -i "$WORK/$n.fa" -m "$WORK/m_$n" -o "$WORK/$n.tsv" --all --cpu 1 --tmpdir "$WORK" 2> "$WORK/$n.cls.log"
  rm -rf "$WORK/m_$n"
}
export -f run_one
export PY WORK
printf '%s\n' $names | xargs -P 4 -I{} bash -c 'run_one {}'
first=1
for n in $names; do
  if [ $first = 1 ]; then cat "$WORK/$n.tsv"; first=0; else tail -n +2 "$WORK/$n.tsv"; fi
done > $OUT/loo_results.tsv
echo "done: $(($(wc -l < $OUT/loo_results.tsv)-1)) rows"
