#!/bin/bash
# Phase B: fetch the A. alternata stress-series RNA-seq.
# Login-node download per CLAUDE.md: NCBI prefetch, not ENA FTP.
# Do NOT use BASH_SOURCE here; the project root is passed in as UPO_ROOT.
set -u
: "${UPO_ROOT:?set UPO_ROOT to the enzymes_UPO directory}"
MAN="$UPO_ROOT/data/g3/rnaseq/runs.tsv"
OUT="$UPO_ROOT/data/g3/rnaseq/sra"
module load sratoolkit/3.2.0
mkdir -p "$OUT"
fail=0
while IFS=$'\t' read -r run cond nreads proj; do
  [ -z "$run" ] && continue
  if [ -s "$OUT/$run/$run.sra" ]; then echo "have $run"; continue; fi
  echo "=== prefetch $run ($cond) ==="
  if ! prefetch "$run" -O "$OUT" --max-size 30G; then
     echo "FAILED $run"; fail=$((fail+1))
  fi
done < "$MAN"
echo "=== done, $fail failures ==="
ls "$OUT" | wc -l
