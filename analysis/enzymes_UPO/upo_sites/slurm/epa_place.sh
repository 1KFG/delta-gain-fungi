#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -J upo_epa
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/epa_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/epa_%j.err
set -euo pipefail
# Evolutionary Placement Algorithm: place each dereplicated BFD UPO onto the 93-reference
# tree. This replaces nearest-neighbour clade assignment, which PLAN.md section 3.1 flagged as
# validated only at 52.7-96.2% identity while BFD sits at 25-35%.
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
P=$PROJ/data/placement
module load raxml/8.2.12
export TMPDIR="${SCRATCH:?SCRATCH not set}/upo_epa.$SLURM_JOB_ID"; mkdir -p "$TMPDIR"
cd "$TMPDIR"
cp "$P/epa_input.fa" aln.fa
cp "$P/reftree.raxml.bestTree" ref.tre
raxmlHPC-PTHREADS -T "${SLURM_CPUS_PER_TASK:-16}" -f v -s aln.fa -t ref.tre \
    -m PROTGAMMALG -n epa -p 1337
cp RAxML_classification.epa          "$P/" 2>/dev/null || true
cp RAxML_portableTree.epa.jplace     "$P/" 2>/dev/null || true
cp RAxML_labelledTree.epa            "$P/" 2>/dev/null || true
cp RAxML_info.epa                    "$P/" 2>/dev/null || true
rm -rf "$TMPDIR"
echo "EPA done $(date +%T)"
