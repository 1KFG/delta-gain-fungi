#!/bin/bash
# Census step 5, part 2 (parallel hmmalign variant of step5_place.sbatch, for a very large candidate set): align every cut domain to the reference, epa-ng placement, clade call,
# assign_steps. Run after step5_prepare.sbatch and only if the references form a clade.
#   PC_STEP=phytoene_synthase PC_MIN_SUPPORT=0.5 sbatch analysis/pathway_census/slurm/step5_place.sbatch
# Required: PC_STEP, PC_MIN_SUPPORT. Optional: PC_PATHWAY, PC_ROOT, PC_MIN_LWR (0.5), PC_MAX_PENDANT (1.0).
#SBATCH -p short
#SBATCH --constraint="ryzen|cascade|broadwell"
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J step5_place_par
#SBATCH -o step5_place_%j.out
set -euo pipefail
: "${PC_STEP:?set PC_STEP}"
: "${PC_MIN_SUPPORT:?set PC_MIN_SUPPORT}"
PC_PATHWAY=${PC_PATHWAY:-carotenoid}
PC_ROOT=${PC_ROOT:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/pathway_census}
PC_MIN_LWR=${PC_MIN_LWR:-0.5}
PC_MAX_PENDANT=${PC_MAX_PENDANT:-1.0}
PC_CAND_DIR=${PC_CAND_DIR:-$PC_ROOT/results/search/${PC_PATHWAY}_full}
OUT=${PC_OUT:-$PC_ROOT/results/step5/$PC_PATHWAY/$PC_STEP}
PY=/usr/bin/python3.12
TOOL=$PC_ROOT/bin/step5_place.py
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
# hmmalign is single-threaded: split the domains into CPUS parts and align the parts in parallel (each part against the same HMM)
mkdir -p $W/parts
awk -v n=$CPUS -v d=$W/parts '/^>/{i++} {print > (d "/part" (i % n) ".faa")}' $OUT/domains.faa
ls $W/parts/part*.faa | xargs -P $CPUS -I{} sh -c 'hmmalign --outformat afa '$OUT'/ref.hmm {} > {}.aln 2>> '$OUT'/hmmalign.log'
cat $W/parts/part*.faa.aln > $W/aln.afa
[ "$(grep -c '^>' $W/aln.afa)" = "$(grep -c '^>' $OUT/domains.faa)" ] || { echo "hmmalign parts lost sequences" >&2; exit 1; }
# match-cols from upo_place: drop insert columns (lowercase and '.')
$PY $UPO match-cols $W/aln.afa --length $NCOL --out $W/q.fa
$PX epa-ng --ref-msa $OUT/ref_set.trim.afa --tree $OUT/ref.raxml.bestTree --query $W/q.fa \
  --model $OUT/ref.raxml.bestModel -w $W --redo -T $CPUS > $OUT/epa-ng.log 2>&1
$PX epa-ng --version > $OUT/epa-ng.version 2>&1 || true
$PY $UPO assign $W/epa_result.jplace --clades $OUT/clade_members.tsv --min-lwr $PC_MIN_LWR \
  --max-pendant $PC_MAX_PENDANT --out $OUT/assign.tsv
gzip -c $W/epa_result.jplace > $OUT/epa_result.jplace.gz
$PY $TOOL placements --assign $OUT/assign.tsv --step $PC_STEP --out $OUT/placements.tsv \
  --out-detail $OUT/placements_detail.tsv
# lineage subclades (v2 only): placements with the subclade as the clade call, counts, nearest labelled leaf
if [ -s $OUT/labelled_leaves.tsv ]; then
  $PY $TOOL placements --assign $OUT/assign.tsv --step $PC_STEP --system fig1 --out $OUT/placements_subclade.tsv \
    --out-detail $OUT/placements_subclade_detail.tsv
  $PY $TOOL subclade-report --tree $OUT/ref.raxml.bestTree --clade-members $OUT/clade_members.tsv \
    --labels $OUT/labelled_leaves.tsv --detail $OUT/placements_detail.tsv --clusters $OUT/clusters.tsv.zst --out-dir $OUT
fi
# candidates for this step only (assign_steps takes a table with a step_id column)
$PY - "$PC_CAND_DIR" $PC_STEP $OUT/candidates_step.tsv <<'PYEOF'
import glob, subprocess, sys, io
d, step, out = sys.argv[1:4]
first = True
with open(out, "w") as fo:
    for p in sorted(glob.glob(f"{d}/candidates.chunk*.PROVISIONAL.tsv.zst")):
        txt = subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout.decode().splitlines()
        hdr = [l for l in txt if not l.startswith("#")][0]
        if first:
            fo.write(hdr + "\n"); first = False
        fo.writelines(l + "\n" for l in txt if not l.startswith("#") and l != hdr and l.split("\t")[-1] == step)
PYEOF
$PY $PC_ROOT/bin/assign_steps.py --pathway $PC_ROOT/pathways/$PC_PATHWAY.yaml \
  --candidates $OUT/candidates_step.tsv --placements $OUT/placements.tsv \
  --min-support $PC_MIN_SUPPORT --out $OUT/assigned.tsv --provisional
$PY $TOOL summarize --out-dir $OUT --assigned $OUT/assigned.PROVISIONAL.tsv --out $OUT/SUMMARY.tsv
echo "step5_place done $(date +%T)"
