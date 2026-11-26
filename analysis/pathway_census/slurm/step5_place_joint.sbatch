#!/bin/bash
# Census step 5, part 2 for a JOINT tree (copy of step5_place.sbatch). Run after step5_prepare_joint.sbatch.
# Places every cut domain once; then, for each step in PC_STEPS (space separated), assigns with that step's clade file
# (clade_members.<S>.tsv), keeps the proteins that have a row of S in PC_STEP_CAND_DIR, and runs assign_steps for S.
# Required: PC_STEP (joint id), PC_STEPS, PC_ROOT, PC_PATHWAY, PC_CAND_DIR (the SAME union folder as part 1),
#   PC_STEP_CAND_DIR (rows per step), PC_MIN_SUPPORT, PC_SCORES (pks_scores.tsv), PC_GENOMES, PC_OUT. Optional: PC_MIN_LWR (0.5), PC_MAX_PENDANT (1.0),
#   PC_EPA_ARGS (extra epa-ng arguments; unset = epa-ng defaults, at most 7 placements per query).
#SBATCH -p intel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -J step5_place_joint
#SBATCH -o step5_place_joint_%j.out
set -euo pipefail
: "${PC_STEP:?set PC_STEP}"
: "${PC_MIN_SUPPORT:?set PC_MIN_SUPPORT}"
PC_PATHWAY=${PC_PATHWAY:-carotenoid}
PC_ROOT=${PC_ROOT:?set PC_ROOT}
PC_MIN_LWR=${PC_MIN_LWR:-0.5}
PC_MAX_PENDANT=${PC_MAX_PENDANT:-1.0}
PC_CAND_DIR=${PC_CAND_DIR:?set PC_CAND_DIR}
: "${PC_STEPS:?set PC_STEPS}" "${PC_STEP_CAND_DIR:?set PC_STEP_CAND_DIR}" "${PC_SCORES:?set PC_SCORES}" "${PC_GENOMES:?set PC_GENOMES}"
OUT=${PC_OUT:-$PC_ROOT/results/step5/$PC_PATHWAY/$PC_STEP}
PY=/usr/bin/python3.12
TOOL=$PC_ROOT/bin/step5_place.py
TOOL_J=$PC_ROOT/bin/step5_joint.py
UPO=$PC_ROOT/../enzymes_UPO/upo_place/bin/upo_place.py
export PATH="$HOME/.pixi/bin:$PATH"
module load hmmer/3.4
module load zstd 2>/dev/null || true
PX="pixi exec -c conda-forge -c bioconda -s epa-ng==0.3.8 -s raxml-ng==2.0.4 --"
W="${SCRATCH:?SCRATCH not set}/step5_place_${PC_STEP}"; mkdir -p "$W"
CPUS=${SLURM_CPUS_PER_TASK:-8}
NCOL=$(awk '/^>/{n++; next} n==1{L+=length($0)} END{print L}' $OUT/ref_set.trim.afa)
hmmbuild --symfrac 0 --informat afa -n ${PC_STEP}_ref $OUT/ref.hmm $OUT/ref_set.trim.afa > $OUT/hmmbuild.log
LEN=$(awk '$1=="LENG"{print $2}' $OUT/ref.hmm)
[ "$LEN" = "$NCOL" ] || { echo "ref.hmm has $LEN match states, trimmed MSA has $NCOL" >&2; exit 1; }
hmmalign --outformat afa $OUT/ref.hmm $OUT/domains.faa > $W/aln.afa 2> $OUT/hmmalign.log
# match-cols from upo_place: drop insert columns (lowercase and '.')
$PY $UPO match-cols $W/aln.afa --length $NCOL --out $W/q.fa
zstd -qf -T0 -c $W/q.fa > $OUT/query_match_cols.fa.zst   # kept so that the rarefaction runs need no new hmmalign
$PX epa-ng --ref-msa $OUT/ref_set.trim.afa --tree $OUT/ref.raxml.bestTree --query $W/q.fa \
  --model $OUT/ref.raxml.bestModel -w $W --redo -T $CPUS ${PC_EPA_ARGS:-} > $OUT/epa-ng.log 2>&1
$PX epa-ng --version > $OUT/epa-ng.version 2>&1 || true
N_DOM=$(grep -c '^>' $OUT/domains.faa)
gzip -c $W/epa_result.jplace > $OUT/epa_result.jplace.gz
for S in $PC_STEPS; do
  if [ ! -s $OUT/clade_members.$S.tsv ]; then echo "no clade file for $S (references not isolated): skipped" | tee -a $OUT/place_skipped.txt; continue; fi
  $PY $UPO assign $W/epa_result.jplace --clades $OUT/clade_members.$S.tsv --min-lwr $PC_MIN_LWR \
    --max-pendant $PC_MAX_PENDANT --out $OUT/assign.$S.tsv
  N_Q=$(($(wc -l < $OUT/assign.$S.tsv) - 1))
  [ "$N_Q" = "$N_DOM" ] || { echo "assign.$S.tsv has $N_Q queries, domains.faa has $N_DOM" >&2; exit 1; }
  $PY $TOOL_J filter-placements --assign $OUT/assign.$S.tsv --step $S --cand-dir $PC_STEP_CAND_DIR \
    --out-plc $OUT/placements.$S.tsv --out-detail $OUT/placements_detail.$S.tsv --out-cand $OUT/candidates_step.$S.tsv | tee -a $OUT/filter_placements.jsonl
  # lesson: candidates_step and placements must have the same proteins
  NC=$(tail -n +2 $OUT/candidates_step.$S.tsv | cut -f1 | sort -u | wc -l)
  NP=$(($(wc -l < $OUT/placements.$S.tsv) - 1))
  [ "$NC" = "$NP" ] || { echo "$S: $NC candidate proteins, $NP placements" >&2; exit 1; }
  $PY $PC_ROOT/bin/assign_steps.py --pathway $PC_ROOT/pathways/$PC_PATHWAY.yaml \
    --candidates $OUT/candidates_step.$S.tsv --placements $OUT/placements.$S.tsv \
    --min-support $PC_MIN_SUPPORT --out $OUT/assigned.$S.tsv --provisional
done
# the joint candidate folder (PC_CAND_DIR) has one row per protein: its count must equal the domains placed
NJ=$($PY - "$PC_CAND_DIR" $PC_STEP <<'PYEOF'
import glob, subprocess, sys
d, step = sys.argv[1:3]
n = 0
for p in sorted(glob.glob(f"{d}/candidates.chunk*.PROVISIONAL.tsv.zst")):
    t = subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout.decode().splitlines()
    rows = [l for l in t if not l.startswith("#")][1:]
    n += sum(1 for l in rows if l.split("\t")[-1] == step)
print(n)
PYEOF
)
[ "$NJ" = "$N_DOM" ] || { echo "joint candidates $NJ rows, domains placed $N_DOM" >&2; exit 1; }
$PY $TOOL_J report --out-dir $OUT --scores $PC_SCORES --genomes $PC_GENOMES
echo "step5_place_joint done $(date +%T): domains placed $N_DOM, joint candidate rows $NJ"
