#!/bin/bash
# UPO clade summary: Figures 1, 3, 4 and Table T1.
#   sbatch analysis/enzymes_UPO/figures/slurm/figures.sh
# Output: analysis/enzymes_UPO/figures/data/*.tsv (+ table_T1.md), figures/out/*.{pdf,png}
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH -J upo_figures
#SBATCH -o logs/figures_%j.out
#SBATCH -e logs/figures_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
BFD=${BFD:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/data/bfd_full/upo_place_bfd_all.tsv.zst}
export PATH="$HOME/.pixi/bin:$PATH"
module load zstd 2>/dev/null || true
$PY -W ignore -m unittest discover -s $UPO/figures/test
$PY $UPO/figures/bin/prep_figures.py --root $UPO --bfd "$BFD" --out $UPO/figures/data
# Pinned R plotting environment (bioconda/conda-forge through pixi).
pixi exec -c conda-forge -c bioconda -s r-base=4.4 -s bioconductor-ggtree=4.0.4 \
    -s r-ggplot2 -s r-ggalluvial=0.12.6 -s r-patchwork -s r-ape -s r-readr -s r-dplyr \
    -s r-tidyr -s r-scales -- \
    Rscript $UPO/figures/bin/plot_figures.R $UPO $UPO/figures/data $UPO/figures/out
echo "done $(date +%T)"
