#!/bin/bash
# DeepLoc 2.1 for the four panel proteomes - the GPU quarter of the four predictors that
# annotate_panels.sh never ran. The other three are in slurm/annotate_panels_part2.sh.
#   sbatch analysis/mucor_secretome/slurm/deeploc_panels.sh
#
# Invocation copied from Fungi_BFD/nextflow/modules/BFD/DEEPLOC/main.nf, with the production
# parameters from nextflow/conf/profile_BFD.config, not test.config:
#   deeploc_model = 'Fast'   (ESM1b)
#   deeploc_gpu   = true     -> device cuda, queue short_gpu, 4 cpus, 32 GB
#   deeploc_sif   = ${singularity_cache}/deeploc-2.1.sif
# so these proteomes get the same model as every BFD genome. Running 'Accurate' here would
# produce numbers that cannot be compared with the database, which is the whole point of
# running it at all.
#
# DeepLoc is licensed software in a local container and the MorenoLab strains are unpublished:
# this runs entirely on cluster GPUs and no sequence leaves the machine. There is no hosted
# DeepLoc endpoint in this pipeline and there must not be.
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J deeploc_panels
#SBATCH -o logs/deeploc_panels_%j.out
#SBATCH -e logs/deeploc_panels_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

A=analysis/mucor_secretome
PY=/usr/bin/python3.12
BFD_BIN=${BFD_BIN:-$HOME/projects/Fungi_BFD/nextflow/bin}
SIF=${DEEPLOC_SIF:-/bigdata/stajichlab/shared/singularity_cache/deeploc-2.1.sif}
MODEL=Fast
OUT=$A/out/panels; mkdir -p $OUT
FAA=$OUT/faa
W="${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}/deeploc_panels"
mkdir -p "$W"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer
[ -s "$SIF" ] || { echo "[ERROR] DeepLoc SIF not found: $SIF" >&2; exit 1; }
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader || true

LOCUSTAGS="${PANEL_LOCUSTAGS:-JEM1MUCOR JEM2MUCOR MPTAK MPALE}"
export OMP_NUM_THREADS=${SLURM_CPUS_ON_NODE:-4}

for LT in $LOCUSTAGS; do
  F=$FAA/$LT.faa
  [ -s "$F" ] || { echo "[ERROR] $F missing - run annotate_panels.sh first" >&2; exit 1; }
  NSEQ=$(grep -c '^>' "$F")
  echo "=== $LT  $NSEQ proteins  $(date +%T)"

  PROT_REAL_DIR=$(dirname "$(readlink -f "$F")")
  rm -rf "$W/out_$LT"
  apptainer exec --nv \
      --bind "$W:$W,$PROT_REAL_DIR:$PROT_REAL_DIR" \
      "$SIF" \
      deeploc2 -f "$F" -o "$W/out_$LT" -m "$MODEL" -d cuda

  # DeepLoc keys results by FASTA id, so duplicate ids silently collapse rows. BFD's module
  # requires one row per input sequence; the same check is made here rather than trusted.
  shopt -s nullglob
  RES=("$W/out_$LT"/results_*.csv)
  [ ${#RES[@]} -eq 1 ] || { echo "[ERROR] expected 1 results CSV, found ${#RES[@]}" >&2; exit 1; }
  NROW=$(( $(wc -l < "${RES[0]}") - 1 ))
  [ "$NROW" -eq "$NSEQ" ] || {
      echo "[ERROR] $LT: $NSEQ sequences in, $NROW rows out" >&2; exit 1; }
  echo "    $NROW rows, one per sequence"

  cp "${RES[0]}" "$W/$LT.deeploc_${MODEL,,}.csv"
  gzip -f "$W/$LT.deeploc_${MODEL,,}.csv"
done

echo "=== merging $(date +%T)"
module load duckdb 2>/dev/null || true
$PY "$BFD_BIN/merge_deeploc.py" $W/*.deeploc_*.csv.gz -o "$OUT/deeploc.parquet"
[ -s "$OUT/deeploc.parquet" ] || { echo "[ERROR] deeploc.parquet is empty" >&2; exit 1; }
ls -l "$OUT/deeploc.parquet"

echo "done $(date +%T)"
