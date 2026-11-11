#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 8:00:00
#SBATCH -J upo_valexp
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/valexp_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/valexp_%j.err
set -euo pipefail
# Does expanding the reference set with high-confidence BFD members actually improve clade
# calling at low identity? Queries are the ORIGINAL 19 characterized enzymes only. Extended
# references above the identity ceiling are removed along with everything else, so an extended
# reference cannot leak the answer for the query it was derived from.
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
module load raxml-ng/1.2.0 raxml/8.2.12
export TMPDIR="${SCRATCH:?SCRATCH not set}/valexp.$SLURM_JOB_ID"; mkdir -p "$TMPDIR"
/usr/bin/python3.12 "$PROJ/upo_sites/bin/validate_placement_matched.py" \
  --ref-aln "$PROJ/data/placement/expanded/ref_aln.filtered.fa" \
  --labels  "$PROJ/data/placement/expanded/ref_labels.tsv" \
  --query-list "$PROJ/data/placement/expanded/orig19.txt" \
  --thresholds 95,70,55,45,40,35 \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --out "$PROJ/data/placement/expanded/matched_identity_validation_expanded.tsv"
rm -rf "$TMPDIR"
echo "expanded validation done $(date +%T)"
