#!/usr/bin/env bash
#SBATCH -J compress_hits
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/logs/compress_hits.%j.out
# Preserve raw DIAMOND hit files (zstd) before the work/ dirs are deleted.
#  - n2000 run (nextflow run focused_jepsen): all 177 nr_cluster + 177 uniprot hit files
#  - full run  (nextflow run dreamy_maxwell, STOPPED): finished DIAMOND_NOVELTY_SCREEN (nr) hit files only
# Each file is zstd-compressed, then decompressed md5 is compared to the original md5.
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow >/dev/null 2>&1 || true
F=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi

do_one() {   # $1 src  $2 destdir
  src=$1; dst=$2/$(basename "$src").zst
  if [[ -s $dst ]]; then :; else zstd -q -3 -f "$src" -o "$dst.part" && mv "$dst.part" "$dst"; fi
  a=$(md5sum < "$src" | cut -d' ' -f1); b=$(zstd -dc "$dst" | md5sum | cut -d' ' -f1)
  [[ $a == $b ]] && ok=1 || ok=0
  printf '%s\t%s\t%s\t%s\t%s\n' "$(basename "$src")" "$(stat -c %s "$src")" "$(stat -c %s "$dst")" "$a" "$ok"
}
export -f do_one

run_set() {  # $1 launchdir  $2 runname  $3 process regex  $4 destdir  $5 manifest
  mkdir -p "$4"
  ( cd "$1"; nextflow log "$2" -f process,status,exit,workdir 2>/dev/null ) \
    | awk -F'\t' -v re="$3" '$1 ~ re && ($2=="CACHED"||$2=="COMPLETED") && $3==0 {print $4}' \
    | while read -r d; do ls "$d"/*.tsv 2>/dev/null; done \
    | xargs -P 8 -I{} bash -c 'do_one "$1" "$2"' _ {} "$4" > "$5"
  echo "$4: $(wc -l < "$5") files, md5 mismatches: $(awk -F'\t' '$5!=1' "$5" | wc -l)"
}

run_set $F focused_jepsen 'DIAMOND_NOVELTY_SCREEN' $F/results/protein_novelty_n2000/raw_hits/nr_cluster  $F/results/protein_novelty_n2000/raw_hits/nr_cluster.manifest.tsv
run_set $F focused_jepsen 'DIAMOND_SEARCH_UNIPROT' $F/results/protein_novelty_n2000/raw_hits/uniprot     $F/results/protein_novelty_n2000/raw_hits/uniprot.manifest.tsv
run_set $F/runs/full dreamy_maxwell 'DIAMOND_NOVELTY_SCREEN' $F/results/protein_novelty_full/raw_hits_partial/nr_cluster $F/results/protein_novelty_full/raw_hits_partial/nr_cluster.manifest.tsv
echo "all done $(date)"
