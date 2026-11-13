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
# short_gpu on an a100. Every a100 is allocated, so this waits whichever queue it is in, but the
# preempt account's priority is ~2.8M against stajichlab's ~103M and its association caps
# gres/gpu at 1, which put an identical job three weeks later in preempt_gpu than in short_gpu.
# The gpu partition is three weeks out. exfab cannot run this at all: its node has only ada6000
# (sm_89), which this torch was not built for, so `-p exfab --gres=gpu:a100:1` is rejected.
# short_gpu caps at 2 h; the ESMFold step skips models already written, so if 640 proteins do not
# finish in one go, resubmitting continues where it stopped.
#SBATCH -p short_gpu
# Any GPU. This used to be `--gres=gpu:a100:1`, because the esm/2.0.0 module's torch is
# 1.13.1+cu117, built for sm_37..sm_86: on the ada6000 (sm_89), h100 (sm_90) or blackwell6000
# (sm_120) nodes it does not refuse - it warns and returns garbage (job 32033882: 487 errors and
# 72 equally worthless "successes"). The default backend is now the project env built by
# scripts/setup_gpu_env.sh (torch 2.8.0+cu128, transformers), which runs on all of them, so the
# queue is 84 GPUs wide instead of 24. bin/esmfold_batch.py still refuses a mismatched card.
# For a fair-esm run, override on the command line: sbatch --gres=gpu:a100:1 ... ESM_BACKEND=fair-esm
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

# --- which ESMFold implementation, and therefore which output folder
# The two backends are the same weights but not bit-identical, so their models are kept in
# separate folders and never mixed in one set. `folds/` is the original fair-esm run (559 models,
# a100 only); `folds_hf/` is the transformers run, which is the default from 2026-11-13 because
# it reaches all 84 GPUs. Downstream files carry the same suffix, so both runs can coexist and
# be compared. See scripts/setup_gpu_env.sh.
ESM_BACKEND=${ESM_BACKEND:-transformers}
ESMFOLD_ENV=${ESMFOLD_ENV:-/bigdata/stajichlab/jstajich/envs/esmfold-torch28}
ESMFOLD_CACHE=${ESMFOLD_CACHE:-/bigdata/stajichlab/jstajich/cache/huggingface}
case $ESM_BACKEND in
  transformers) FOLDS=$OUT/folds_hf; SUF=_hf ;;
  fair-esm)     FOLDS=$OUT/folds;    SUF= ;;
  *) echo "unknown ESM_BACKEND=$ESM_BACKEND (transformers|fair-esm)" >&2; exit 2 ;;
esac
echo "backend $ESM_BACKEND -> $FOLDS"

$PY -W ignore -m unittest discover -s $A/test

# --- 1. candidates, from the tables the predictors wrote
$PY $A/bin/call_effectors.py --dir $PANELS --locustags JEM1MUCOR,JEM2MUCOR \
    --fasta $PANELS/faa/JEM1MUCOR.faa,$PANELS/faa/JEM2MUCOR.faa \
    --codons $PANELS/jem_codon_check.tsv --out $OUT/jem

# The BFD Mucor reference goes through the same code on the BFD export, as a non-Marchantia
# outgroup. Its tables are the sql/export.sql output, which has the same columns.
# Prefer the restored export tables, so the reference goes through exactly the code path the
# JEM genomes do; fall back to the per-protein table from the earlier run if the cache is not
# rebuilt yet (slurm/secretome.sh).
if [ -s $A/cache/signalp.tsv ] && [ -s $A/cache/candidates.faa ]; then
  $PY $A/bin/call_effectors.py --dir $A/cache --locustags F00CDB58 \
      --fasta $A/cache/candidates.faa --codons $A/out/codon_check.tsv --out $OUT/bfd
elif [ -s $A/out/focal_secretome.tsv ] && [ -s $A/out/screen_set.faa ]; then
  echo "cache not built; using the per-protein table from the earlier run" >&2
  $PY $A/bin/call_effectors.py --secretome $A/out/focal_secretome.tsv \
      --fasta $A/out/screen_set.faa --out $OUT/bfd
else
  echo "no BFD tables and no focal secretome table - skipping the reference genome" >&2
fi

cat $OUT/jem/to_fold.faa $OUT/bfd/to_fold.faa 2>/dev/null > $W/to_fold.faa || \
  cp $OUT/jem/to_fold.faa $W/to_fold.faa
cp $W/to_fold.faa $OUT/to_fold.faa
echo "to fold: $(grep -c '>' $OUT/to_fold.faa)"

# --- 2. ESMFold, resumable: models already in $FOLDS are not refolded
# The esm module activates a conda env with `conda activate --stack`, which is the shell
# function, not the conda executable - without the hook it fails with
# "conda: error: unrecognized arguments: --stack" (job 32033764).
export CONDA_NO_PLUGINS=true
module load miniconda3
eval "$(conda shell.bash hook)"
if [ "$ESM_BACKEND" = transformers ]; then
  conda activate "$ESMFOLD_ENV"
  # Weights are prefetched to $HF_HOME on the login node by scripts/setup_gpu_env.sh, because
  # compute nodes may have no outbound network.
  export HF_HOME="$ESMFOLD_CACHE"
  export TOKENIZERS_PARALLELISM=false
  RUNPY=python
else
  module load esm/2.0.0
  RUNPY=python3
fi
nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader || true
$RUNPY $A/bin/esmfold_batch.py --fasta $OUT/to_fold.faa --out $FOLDS --max-len 400 \
    --backend $ESM_BACKEND
echo "models: $(ls $FOLDS/*.pdb 2>/dev/null | wc -l)"

# If the other backend's run exists, compare them on the proteins both folded. Same weights,
# different implementation: pLDDT should track closely, and a systematic gap is a bug, not a
# finding. Costs nothing and is the cheapest check that the swap was sound.
OTHER=$OUT/folds/plddt.tsv
[ "$ESM_BACKEND" = fair-esm ] && OTHER=$OUT/folds_hf/plddt.tsv
if [ -s "$OTHER" ] && [ -s "$FOLDS/plddt.tsv" ]; then
  echo "=== backend comparison (pLDDT, 0-100, shared proteins) ==="
  $PY - "$FOLDS/plddt.tsv" "$OTHER" <<'PYEOF'
import sys, statistics as st
def load(p):
    d = {}
    with open(p) as fh:
        head = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                d[f[0]] = float(f[2])
    return d
a, b = load(sys.argv[1]), load(sys.argv[2])
both = sorted(set(a) & set(b))
print(f"this backend {len(a)} models, other {len(b)}, shared {len(both)}")
if both:
    xa = [a[k] for k in both]; xb = [b[k] for k in both]
    d = [x - y for x, y in zip(xa, xb)]
    print(f"mean pLDDT this {st.mean(xa):.1f} vs other {st.mean(xb):.1f}; "
          f"mean diff {st.mean(d):+.1f}, median diff {st.median(d):+.1f}, "
          f"max |diff| {max(abs(v) for v in d):.1f}")
    # the number that actually matters downstream: does the >= 70 call flip?
    flip = sum((x >= 70) != (y >= 70) for x, y in zip(xa, xb))
    print(f"proteins crossing the pLDDT 70 line between backends: {flip} of {len(both)}")
PYEOF
fi

# --- 3. Foldseek against the three structure databases
module unload esm/2.0.0 2>/dev/null || true
module load foldseek
FMT="query,target,theader,fident,alnlen,qcov,tcov,evalue,bits,alntmscore,prob"
HITARG=""
for db in pdb cath afproteome; do
  [ -s $FSDB/$db.dbtype ] || { echo "no foldseek db $db at $FSDB" >&2; continue; }
  foldseek easy-search $FOLDS $FSDB/$db $W/hits_$db$SUF.tsv $W/tmp_$db \
      --format-output "$FMT" --alignment-type 1 -e 10 --max-seqs 1000 \
      --threads $CPU -v 1
  cp $W/hits_$db$SUF.tsv $OUT/hits_$db$SUF.tsv
  HITARG="$HITARG,$db=$OUT/hits_$db$SUF.tsv"
done

# --- 4. one fold call per candidate
$PY $A/bin/fold_report.py --hits "${HITARG#,}" --plddt $FOLDS/plddt.tsv \
    --candidates $OUT/jem/candidates.tsv --out $OUT/fold_calls$SUF.tsv
echo "=== fold calls, best 25 by TM-score ==="
head -26 $OUT/fold_calls$SUF.tsv | cut -f1,2,5,6,9 | column -t -s $'\t'

# --- 5. host panel, built identically for both liverworts (the isolation host is unrecorded)
$PY $A/bin/plant_panel.py --dir $PANELS --locustags MPTAK,MPALE --out $OUT/plant_panel || \
  echo "plant panel skipped (Marchantia tables not ready)" >&2
if [ -s $OUT/plant_panel/compare.tsv ]; then
  echo "=== host target classes: M. polymorpha vs M. paleacea ==="
  column -t -s $'\t' $OUT/plant_panel/compare.tsv | cut -c1-150
fi
echo "done $(date +%T)"
