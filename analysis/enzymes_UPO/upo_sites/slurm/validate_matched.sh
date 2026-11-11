#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 6:00:00
#SBATCH -J upo_valmatch
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/valmatch_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/valmatch_%j.err
set -euo pipefail
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
module load raxml-ng/1.2.0 raxml/8.2.12
export TMPDIR="${SCRATCH:?SCRATCH not set}/valmatch.$SLURM_JOB_ID"; mkdir -p "$TMPDIR"
/usr/bin/python3.12 "$PROJ/upo_sites/bin/validate_placement_matched.py" \
  --ref-aln "$PROJ/data/placement/ref_aln.filtered.fa" \
  --labels  "$PROJ/data/placement/ref_labels.tsv" \
  --thresholds 95,70,55,45,40,35 \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --out "$PROJ/data/placement/matched_identity_validation.tsv"
rm -rf "$TMPDIR"
echo "matched validation done $(date +%T)"
