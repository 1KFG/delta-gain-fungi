#!/bin/bash
# Merge the place_bfd.sh chunk outputs and compare with the old upo_sites BFD calls.
#   sbatch analysis/enzymes_UPO/upo_place/slurm/merge_bfd.sh     (after place_bfd.sh)
# Output (this clone): data/bfd_place/{bfd_place.tsv.zst,COMPARE.md,G1_BFD_SUMMARY.md}
# G1_BFD_SUMMARY.md is copied to analysis/enzymes_UPO/ by hand after review.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH -J upo_place_merge
#SBATCH -o logs/place_merge_%j.out
#SBATCH -e logs/place_merge_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source analysis/enzymes_UPO/upo_place/slurm/env.sh
module load zstd 2>/dev/null || true
OLD=${OLD:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full/upo_sites_bfd_all.tsv.zst}
OUT=$UPO/data/bfd_place
$PY $UPO/upo_place/bin/compare_bfd.py --placed $OUT/chunk_*.tsv --old "$OLD" \
    --out-merged $OUT/bfd_place.tsv.zst --out-report $OUT/COMPARE.md
$PY $UPO/upo_place/bin/summarize_place.py --table $OUT/bfd_place.tsv.zst \
    --source-label "data/bfd_full/upo_place_bfd_all.tsv.zst" --out $OUT/G1_BFD_SUMMARY.md
echo "done $(date +%T)"
