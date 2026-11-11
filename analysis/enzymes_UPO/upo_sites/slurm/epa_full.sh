#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 10:00:00
#SBATCH -J upo_epa_full
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/epafull_%A_%a.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/epafull_%A_%a.err
set -euo pipefail

# Phylogenetic placement across ALL PF01328 proteins in BFD, not just the 7,055 dereplicated
# representatives. The reference tree is FIXED, and PF01328 match-state columns are fixed by the
# profile, so each chunk can be placed independently and the results concatenated.

PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
P=$PROJ/data/placement
CH=$PROJ/data/bfd_full/chunks
OUT=$P/full
PER_TASK=${PER_TASK:-4}

module load hmmer/3.4 raxml/8.2.12
mkdir -p "$OUT"
export TMPDIR="${SCRATCH:?SCRATCH not set}/epafull.${SLURM_ARRAY_JOB_ID:-x}.${SLURM_ARRAY_TASK_ID:-1}"
mkdir -p "$TMPDIR"

TASK=${SLURM_ARRAY_TASK_ID:-1}
mapfile -t ALL < <(ls "$CH"/chunk_*.fa | sort)
START=$(( (TASK - 1) * PER_TASK ))
if [ "$START" -ge "${#ALL[@]}" ]; then echo "task $TASK: nothing to do"; exit 0; fi
MINE=( "${ALL[@]:$START:$PER_TASK}" )
echo "task $TASK handling ${#MINE[@]} chunks"

for f in "${MINE[@]}"; do
  base=$(basename "$f" .fa)
  final="$OUT/${base}.clades.tsv"
  if [ -s "$final" ]; then echo "skip $base"; continue; fi
  echo "=== $base $(date +%T) ==="
  w="$TMPDIR/$base"; mkdir -p "$w"

  /usr/bin/python3.12 "$PROJ/upo_sites/bin/make_placement_aln.py" \
      --hmm "$PROJ/upo_id/models/PF01328.hmm" \
      --refs "$P/refs.fa" --queries "$f" \
      --out-ref "$w/ref.fa" --out-all "$w/all.fa"

  # keep only the 93 references that survive the coverage filter, so the alignment matches the tree
  /usr/bin/python3.12 - "$w/all.fa" "$P/ref_aln.filtered.fa" "$w/epa.fa" <<'PY'
import sys
def rd(p):
    d={};n=None
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=''
        else: d[n]+=l.strip()
    return d
allf, reff, out = sys.argv[1], sys.argv[2], sys.argv[3]
keep=set(rd(reff))
a=rd(allf)
with open(out,'w') as fh:
    for n,s in a.items():
        if n.startswith('NEG_'): continue
        fh.write(f">{n}\n{s}\n")
PY

  cd "$w"
  raxmlHPC-PTHREADS -T "${SLURM_CPUS_PER_TASK:-8}" -f v -s epa.fa \
      -t "$P/reftree.raxml.bestTree" -m PROTGAMMALG -n epa -p 1337 > raxml.log 2>&1

  /usr/bin/python3.12 "$PROJ/upo_sites/bin/jplace_to_clade.py" \
      --jplace "$w/RAxML_portableTree.epa.jplace" \
      --labels "$P/ref_labels.tsv" --out "$final" 2> "$OUT/${base}.log"
  echo "=== $base done $(date +%T) ==="
  rm -rf "$w"
done
rm -rf "$TMPDIR"
echo "task $TASK complete"
