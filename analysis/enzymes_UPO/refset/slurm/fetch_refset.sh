#!/bin/bash
# Fetch the 174 Zamocky 2022 reference sequences (step 1 of the reference-set expansion).
# Network-bound and light: one core; ~15-30 min, paced by the Wayback Machine (3 s between requests). Run from a clone of the repo:
#   sbatch analysis/enzymes_UPO/refset/slurm/fetch_refset.sh
# Output: analysis/enzymes_UPO/data/refset_zamocky2022/fetched/{zamocky2022_refs.faa,fetch_manifest.tsv}
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 2:00:00
#SBATCH -J upo_refset_fetch
#SBATCH -o logs/refset_fetch_%j.out
#SBATCH -e logs/refset_fetch_%j.err
set -euo pipefail
# No BASH_SOURCE under SLURM: run from the repo root (sbatch's working directory).
REPO=${SLURM_SUBMIT_DIR:-$PWD}
cd "$REPO"
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12

$PY -m unittest discover -s $UPO/refset/test
CACHE="${SCRATCH:?SCRATCH not set}/refset_cache"
$PY $UPO/refset/bin/fetch_refset.py \
    --table $UPO/data/refset_zamocky2022/tableS1_resolved.tsv \
    --out   $UPO/data/refset_zamocky2022/fetched \
    --cache "$CACHE"
echo "done $(date +%T)"
