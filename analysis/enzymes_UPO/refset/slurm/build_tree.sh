#!/bin/bash
# Steps 2-4 of the reference-set expansion: dedup + QC the Zamocky 2022 sequences against the
# 19-enzyme backbone, align, and infer the in-house reference tree. Run from a clone of the repo:
#   sbatch analysis/enzymes_UPO/refset/slurm/build_tree.sh
# Output (small, committed): analysis/enzymes_UPO/data/refset_zamocky2022/tree/
#   dedup.tsv qc.tsv tree_input.faa aln.afa aln.trim.afa
#   mfp.*      IQ-TREE, ModelFinder + 1000 UFBoot + 1000 SH-aLRT (the tree we use)
#   lgig.*     IQ-TREE, fixed LG+I+G4 (Zamocky's reported model; a check, not the result)
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -J upo_refset_tree
#SBATCH -o logs/refset_tree_%j.out
#SBATCH -e logs/refset_tree_%j.err
set -euo pipefail
# No BASH_SOURCE under SLURM: run from the repo root (sbatch's working directory).
REPO=${SLURM_SUBMIT_DIR:-$PWD}
cd "$REPO"
UPO=analysis/enzymes_UPO
REF=$UPO/data/refset_zamocky2022
OUT=$REF/tree
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_PER_TASK:-4}
SEED=20261108

module load hmmer/3.4 mafft/7.505 trimal/1.5.1 iqtree/3.0.1
mkdir -p "$OUT"
TMP="${SCRATCH:?SCRATCH not set}/refset_tree"
mkdir -p "$TMP"

$PY -m unittest discover -s $UPO/refset/test

# 2. combine + dedup by sequence (backbone copy kept)
$PY $UPO/refset/bin/build_refset.py combine \
    --zamocky  $REF/fetched/zamocky2022_refs.faa \
    --table    $REF/tableS1_resolved.tsv \
    --backbone $UPO/data/clade_backbone/raw.fa \
    --out-fasta "$TMP/combined.faa" --out-dedup $OUT/dedup.tsv

# 3. QC: PF01328 hit (inclusive E-value, not the gathering threshold, so divergent HTPs are
#    scored rather than silently dropped) and HspUPO-numbered catalytic residues
$PY $UPO/upo_sites/bin/upo_sites.py --proteome "$TMP/combined.faa" \
    --hmm $UPO/upo_id/models/PF01328.hmm --backbone $UPO/data/clade_backbone/raw.fa \
    --no-ga --evalue 1e-5 --cpu "$CPU" --out $OUT/upo_sites.tsv
$PY $UPO/refset/bin/build_refset.py select \
    --combined "$TMP/combined.faa" --dedup $OUT/dedup.tsv --sites $OUT/upo_sites.tsv \
    --table $REF/tableS1_resolved.tsv --out-qc $OUT/qc.tsv --out-fasta $OUT/tree_input.faa

# 4. align (L-INS-i), trim to 75 % site coverage (Zamocky's criterion), ML tree
mafft --localpair --maxiterate 1000 --thread "$CPU" $OUT/tree_input.faa > $OUT/aln.afa
trimal -in $OUT/aln.afa -out $OUT/aln.trim.afa -gt 0.75
echo "trimmed alignment: $(grep -c '>' $OUT/aln.trim.afa) seqs," \
     "$(awk '/^>/{if(n++)exit;next}{l+=length($0)}END{print l}' $OUT/aln.trim.afa) columns"

iqtree3 -s $OUT/aln.trim.afa -m MFP -B 1000 -alrt 1000 -T "$CPU" --seed $SEED \
    --prefix $OUT/mfp -redo
iqtree3 -s $OUT/aln.trim.afa -m LG+I+G4 -B 1000 -alrt 1000 -T "$CPU" --seed $SEED \
    --prefix $OUT/lgig -redo
rm -f $OUT/*.ckp.gz $OUT/*.uniqueseq.phy
echo "done $(date +%T)"
