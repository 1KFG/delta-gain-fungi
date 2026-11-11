#!/usr/bin/env bash
# Regression test for upo_scan.py.
#
# WHAT IT PROVES
#   1. upo_scan runs end to end on a proteome FASTA.
#   2. Its t1_pocket_matches and calib_domain calls match, protein for protein, what the full
#      BFD pipeline produced for the same proteins by a different route (upo_sites.py ->
#      merge_final.py -> apply_t1_score.py -> applicability.py).
#   3. A raw backbone and an aligned (.afa) backbone give byte-identical output.
#
# Point 2 is the one that matters. The 87 fixture proteins are real BFD proteins, so their
# expected values come from an independent run of the pipeline rather than from this script.
#
# UPO_ROOT must be set. Do not use BASH_SOURCE: it breaks under SLURM.
set -euo pipefail
: "${UPO_ROOT:?set UPO_ROOT to the enzymes_UPO directory}"
cd "$UPO_ROOT"

PY=/usr/bin/python3.12
module load hmmer/3.4 2>/dev/null || true

WORK="${SCRATCH:-${TMPDIR:-/tmp}}/upo_scan_test.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

FIXTURE=upo_sites/test/bfd87.fa
EXPECTED=data/bfd_full/upo_sites_bfd_domain.tsv.zst

for bb in data/clade_backbone/raw.fa data/clade_backbone/backbone.afa.zst; do
  $PY -P upo_sites/bin/upo_scan.py \
      --proteome "$FIXTURE" \
      --hmm upo_id/models/PF01328.hmm \
      --backbone "$bb" \
      --s2 data/calib_dl/active_site_S2.tsv \
      --out "$WORK/$(basename "$bb").tsv" \
      --cpu 2 >"$WORK/$(basename "$bb").log" 2>&1
done

# 3. both backbone forms agree
if ! diff <(grep -v '^#' "$WORK/raw.fa.tsv") \
          <(grep -v '^#' "$WORK/backbone.afa.zst.tsv") >/dev/null; then
  echo "FAIL: raw and aligned backbones give different output"; exit 1
fi
echo "PASS: raw and aligned backbone forms give identical output"

# 2. agreement with the full BFD pipeline
col () { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
                     {print $h["protein_id"]"\t"$h["t1_pocket_matches"]"\t"$h["calib_domain"]}'; }
grep -v '^#' "$WORK/raw.fa.tsv" | col | awk -F'\t' '{split($1,a,"|"); print a[1]"\t"$2"\t"$3}' \
  | sort > "$WORK/got.txt"
zstd -dcq "$EXPECTED" | grep -v '^#' | col | sort > "$WORK/want.txt"

RES=$(join -t$'\t' "$WORK/got.txt" "$WORK/want.txt" \
      | awk -F'\t' '{n++; if($2==$4 && $3==$5) ok++; else print "  MISMATCH",$1,$2,$3,"vs",$4,$5}
                    END{print "AGREE "ok" "n}')
echo "$RES" | grep MISMATCH || true
STAT=$(echo "$RES" | awk '/^AGREE/{print $2" "$3}')
OK=${STAT% *}; N=${STAT#* }
if [ -z "$N" ] || [ "$N" -eq 0 ]; then
  echo "FAIL: no proteins joined against the expected table"; exit 1
fi
if [ "$OK" != "$N" ]; then
  echo "FAIL: only $OK/$N proteins match the BFD pipeline"; exit 1
fi
echo "PASS: $OK/$N proteins match the full BFD pipeline exactly"
echo "ALL TESTS PASSED"
