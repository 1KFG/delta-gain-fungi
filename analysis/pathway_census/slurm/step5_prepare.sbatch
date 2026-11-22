#!/bin/bash
# Census step 5, part 1: cut domains, cluster, reference alignment, reference tree, clades.
# All paths come from environment variables (no BASH_SOURCE). Submit from the worktree root:
#   PC_STEP=phytoene_synthase sbatch analysis/pathway_census/slurm/step5_prepare.sbatch
# Required: PC_STEP. Optional: PC_PATHWAY (carotenoid), PC_MAX_REPS (300), PC_OUTGROUP_FA (a FASTA path,
# else the tree is rooted on the decoys), PC_ROOT, PC_SEQDB.
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J step5_prepare
#SBATCH -o step5_prepare_%j.out
set -euo pipefail
: "${PC_STEP:?set PC_STEP}"
PC_PATHWAY=${PC_PATHWAY:-carotenoid}
PC_ROOT=${PC_ROOT:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/pathway_census}
PC_SEQDB=${PC_SEQDB:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/seqdb}
PC_MAX_REPS=${PC_MAX_REPS:-300}
PC_CAND_DIR=${PC_CAND_DIR:-$PC_ROOT/results/search/${PC_PATHWAY}_full}
OUT=${PC_OUT:-$PC_ROOT/results/step5/$PC_PATHWAY/$PC_STEP}
PY=/usr/bin/python3.12
TOOL=$PC_ROOT/bin/step5_place.py
export PATH="$HOME/.pixi/bin:$PATH"
# module scripts are not safe under set -e / set -u; check the tools afterwards
set +eu
module load samtools/1.19.2
module load mmseqs2/17-b804f
module load mafft/7.505
module load clipkit/1.3.0
module load hmmer/3.4
set -eu
# the clipkit module loads miniconda3, whose `conda activate` fails in a batch shell (MEASURED, job 32376289);
# the install has its own interpreter, so fall back to its bin directory
command -v clipkit > /dev/null || export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/clipkit/1.3.0/bin:$PATH
for t in samtools mmseqs mafft clipkit hmmbuild zstd pixi; do
  command -v $t > /dev/null || { echo "tool not found after module load: $t" >&2; exit 1; }
done
PX="pixi exec -c conda-forge -c bioconda -s epa-ng==0.3.8 -s raxml-ng==2.0.4 --"
W="${SCRATCH:?SCRATCH not set}/step5_${PC_STEP}"; mkdir -p "$W" "$OUT"
REFS=$PC_ROOT/refs/$PC_PATHWAY/$PC_STEP.faa
DECOYS=$PC_ROOT/refs/$PC_PATHWAY/decoys/$PC_STEP.faa
OG=${PC_OUTGROUP_FA:-}
CPUS=${SLURM_CPUS_PER_TASK:-8}

# 1. cut domains from the BFD protein files
$PY $TOOL extract --candidates $PC_CAND_DIR/candidates.chunk*.PROVISIONAL.tsv.zst --step $PC_STEP \
  --seqdb $PC_SEQDB --out-fasta $OUT/domains.faa --out-coords $OUT/domain_coords.tsv \
  --out-fail $OUT/extract_failures.tsv --out-json $OUT/extract_summary.json --workers $CPUS

# 2. cluster at 50% identity, 80% coverage (both directions)
mmseqs easy-cluster $OUT/domains.faa $W/clu $W/tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads $CPUS
cp $W/clu_cluster.tsv $OUT/clusters.tsv
mmseqs version > $OUT/mmseqs.version
zstd -qf -T0 $OUT/clusters.tsv
$PY $TOOL select-reps --clusters $OUT/clusters.tsv.zst --genomes $PC_ROOT/results/genomes.tsv.zst \
  --max-reps $PC_MAX_REPS --seed 1 --out-ids $OUT/rep_ids.txt --out-report $OUT/reps_by_phylum.tsv \
  --out-json $OUT/reps_summary.json
$PY - $OUT/domains.faa $OUT/rep_ids.txt $OUT/reps.faa <<'PYEOF'
import sys
ids = {l.strip() for l in open(sys.argv[2]) if l.strip()}
keep = False
with open(sys.argv[1]) as fi, open(sys.argv[3], "w") as fo:
    for line in fi:
        if line.startswith(">"):
            keep = line[1:].split()[0] in ids
        if keep:
            fo.write(line)
PYEOF

# 3. alignment set: refs + decoys + outgroup + representatives; ids made safe for raxml/epa
$PY - $REFS $DECOYS "$OG" $OUT/reps.faa $OUT/ref_set.faa <<'PYEOF'
import sys
refs, decoys, og, reps, out = sys.argv[1:6]
with open(out, "w") as fo:
    for p in (refs, decoys, og, reps):
        if not p:
            continue
        for line in open(p):
            fo.write((">" + line[1:].split()[0] + "\n") if line.startswith(">") else line)
PYEOF
mafft --auto --thread $CPUS $OUT/ref_set.faa > $OUT/ref_set.aln.faa 2> $OUT/mafft.log
mafft --version 2> $OUT/mafft.version || true
clipkit $OUT/ref_set.aln.faa -m smart-gap -o $OUT/ref_set.trim.afa > $OUT/clipkit.log 2>&1
clipkit --version > $OUT/clipkit.version 2>&1 || true

# 4. reference tree, LG+G4 (epa-ng rejects some matrices). raxml-ng uses few threads on a small MSA.
$PX raxml-ng --check --msa $OUT/ref_set.trim.afa --model LG+G4 --prefix $W/chk --redo > $OUT/raxml_check.log 2>&1 || true
$PX raxml-ng --msa $OUT/ref_set.trim.afa --model LG+G4 --prefix $OUT/ref --threads 4 --seed 1 \
  --tree pars{5},rand{5} --redo > $OUT/raxml.log 2>&1
rm -f $OUT/ref.raxml.startTree $OUT/ref.raxml.rba $OUT/ref.raxml.reduced.phy

# 5. clades; exit code 3 means the references do not form a clade free of decoys (a result)
OGARG=""; [ -n "$OG" ] && OGARG="--outgroup $OG"
set +e
$PY $TOOL clades --tree $OUT/ref.raxml.bestTree --refs $REFS --decoys $DECOYS $OGARG \
  --reps $OUT/rep_ids.txt --step $PC_STEP --out-dir $OUT
rc=$?
set -e
echo "step5_prepare done rc=$rc $(date +%T)"
exit 0
