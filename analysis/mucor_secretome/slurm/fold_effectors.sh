#!/bin/bash
# Fold-first: what ARE the candidate effectors of JEM1, JEM2 and the BFD Mucor reference?
#   sbatch analysis/mucor_secretome/slurm/fold_effectors.sh
# These proteins were selected because their sequence says nothing - small, secreted, no Pfam
# domain. ESMFold is single-sequence, so an orphan effector folds as well as a conserved one and
# no sequence leaves the cluster; Foldseek against PDB, CATH and the AlphaFold proteome set then
# says whether the fold is a known effector scaffold (LysM, MAX, cystatin, KP6, ...) or new.
#
# Runs on exfab (2 x RTX 6000 Ada). The Foldseek half is CPU and could be split out, but the
# whole thing is short enough to keep in one job, and ESMFold must finish first anyway.
# Reads shared folders read-only; writes mucor_secretome/out/effectors/.
# short_gpu, where the free GPUs are (see slurm/signalp_gpu.sh); a multi-partition request is
# refused by the account association. short_gpu caps at 2 h, which is ample for ESMFold on a few
# hundred short proteins, and the job is resumable - models already written are skipped - so a
# timeout costs only a resubmit.
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1 -n 1 -c 16
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -J fold_eff
#SBATCH -o logs/fold_eff_%j.out
#SBATCH -e logs/fold_eff_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-16}
PANELS=$A/out/panels
OUT=$A/out/effectors; mkdir -p $OUT
FSDB=/srv/projects/db/foldseek
W="${SCRATCH:?}/fold_eff"; mkdir -p $W

$PY -W ignore -m unittest discover -s $A/test

# --- 1. candidates, from the tables the predictors wrote
$PY $A/bin/call_effectors.py --dir $PANELS --locustags JEM1MUCOR,JEM2MUCOR \
    --fasta $PANELS/faa/JEM1MUCOR.faa,$PANELS/faa/JEM2MUCOR.faa \
    --codons $PANELS/jem_codon_check.tsv --out $OUT/jem

# The BFD Mucor reference goes through the same code on the BFD export, as a non-Marchantia
# outgroup. Its tables are the sql/export.sql output, which has the same columns.
if [ -s $A/out/signalp.tsv ] && [ -s $A/out/candidates.faa ]; then
  $PY $A/bin/call_effectors.py --dir $A/out --locustags F00CDB58 \
      --fasta $A/out/candidates.faa --out $OUT/bfd
else
  echo "no BFD export in $A/out - skipping the reference genome" >&2
fi

cat $OUT/jem/to_fold.faa $OUT/bfd/to_fold.faa 2>/dev/null > $W/to_fold.faa || \
  cp $OUT/jem/to_fold.faa $W/to_fold.faa
cp $W/to_fold.faa $OUT/to_fold.faa
echo "to fold: $(grep -c '>' $OUT/to_fold.faa)"

# --- 2. ESMFold, resumable: models already in out/effectors/folds are not refolded
module load esm/2.0.0
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
python3 $A/bin/esmfold_batch.py --fasta $OUT/to_fold.faa --out $OUT/folds --max-len 400
echo "models: $(ls $OUT/folds/*.pdb 2>/dev/null | wc -l)"

# --- 3. Foldseek against the three structure databases
module unload esm/2.0.0 2>/dev/null || true
module load foldseek
FMT="query,target,theader,fident,alnlen,qcov,tcov,evalue,bits,alntmscore,prob"
HITARG=""
for db in pdb cath afproteome; do
  [ -s $FSDB/$db.dbtype ] || { echo "no foldseek db $db at $FSDB" >&2; continue; }
  foldseek easy-search $OUT/folds $FSDB/$db $W/hits_$db.tsv $W/tmp_$db \
      --format-output "$FMT" --alignment-type 1 -e 10 --max-seqs 1000 \
      --threads $CPU -v 1
  cp $W/hits_$db.tsv $OUT/hits_$db.tsv
  HITARG="$HITARG,$db=$OUT/hits_$db.tsv"
done

# --- 4. one fold call per candidate
$PY $A/bin/fold_report.py --hits "${HITARG#,}" --plddt $OUT/folds/plddt.tsv \
    --candidates $OUT/jem/candidates.tsv --out $OUT/fold_calls.tsv
echo "=== fold calls, best 25 by TM-score ==="
head -26 $OUT/fold_calls.tsv | cut -f1,2,5,6,9 | column -t -s $'\t'

# --- 5. host panel, built identically for both liverworts (the isolation host is unrecorded)
$PY $A/bin/plant_panel.py --dir $PANELS --locustags MPTAK,MPALE --out $OUT/plant_panel || \
  echo "plant panel skipped (Marchantia tables not ready)" >&2
if [ -s $OUT/plant_panel/compare.tsv ]; then
  echo "=== host target classes: M. polymorpha vs M. paleacea ==="
  column -t -s $'\t' $OUT/plant_panel/compare.tsv | cut -c1-150
fi
echo "done $(date +%T)"
