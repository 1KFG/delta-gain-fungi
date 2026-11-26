#!/bin/bash
# Placement-only rerun of a joint PKS run with epa-ng keeping all placements, per-step decoy rule, into <run>_allplace.
# Usage: PC_ROOT=<analysis/pathway_census> submit_allplace_joint.sh <old run folder> <cpus>
# Copies the tree, alignment, representatives and clusters (not the old placements), runs `step5_joint.py clades --decoy-scope per-step`,
# then submits slurm/step5_place_joint.sbatch with PC_EPA_ARGS. Prints the job id. Paths from the environment (no BASH_SOURCE).
set -euo pipefail
: "${PC_ROOT:?}"
OLD=${1%/}; cpus=$2; NEW=${OLD}_allplace
[ ! -e "$NEW" ] || { echo "$NEW exists; not overwritten" >&2; exit 1; }
mkdir -p "$NEW/inputs"
cp "$OLD"/inputs/refs_*.faa "$OLD"/inputs/decoys_all.faa "$OLD"/inputs/named_members.faa "$OLD"/inputs/named_leaves.tsv "$NEW/inputs/"
for f in ref_set.trim.afa ref_set.aln.faa ref_set.faa ref.raxml.bestTree ref.raxml.bestModel ref.raxml.log rep_ids.txt reps.faa reps_by_phylum.tsv \
         reps_summary.json clusters.tsv.zst domain_coords.tsv extract_summary.json extract_failures.tsv; do cp "$OLD/$f" "$NEW/$f"; done
if [ -s "$OLD/domains.faa.zst" ]; then zstd -dc "$OLD/domains.faa.zst" > "$NEW/domains.faa"; else cp "$OLD/domains.faa" "$NEW/domains.faa"; fi
I=$NEW/inputs; D=$PC_ROOT/refs/melanin_dhn/decoys
export PC_PATHWAY=melanin_dhn PC_STEP=pks_joint PC_OUT=$NEW PC_ROOT
export PC_CAND_DIR=$PC_ROOT/results/search/pks_joint_union PC_STEP_CAND_DIR=$PC_ROOT/results/search/pks_joint
export PC_STEPS="thn_synthase hexaketide_pks heptaketide_pks" PC_MIN_SUPPORT=0.5 PC_SCORES=$PC_ROOT/results/search/pks_joint/pks_scores.tsv
export PC_GENOMES=$PC_ROOT/results/manifests/genomes.tsv.zst PC_EPA_ARGS="--filter-acc-lwr 0.999 --filter-max 10000"
/usr/bin/python3.12 "$PC_ROOT/bin/step5_joint.py" clades --decoy-scope per-step \
  --step-decoys thn_synthase=$D/thn_synthase.faa hexaketide_pks=$D/hexaketide_pks.faa heptaketide_pks=$D/heptaketide_pks.faa \
  --tree "$NEW/ref.raxml.bestTree" --refs thn_synthase=$I/refs_thn_synthase.faa hexaketide_pks=$I/refs_hexaketide_pks.faa heptaketide_pks=$I/refs_heptaketide_pks.faa \
  --decoys "$I/decoys_all.faa" --members-fa "$I/named_members.faa" --reps "$NEW/rep_ids.txt" --joint-id pks_joint --out-dir "$NEW" > "$NEW/clades.log" 2>&1 || {
  rc=$?; [ $rc = 3 ] || { echo "clades failed rc $rc, see $NEW/clades.log" >&2; exit $rc; }; }
sbatch --parsable --export=ALL -p intel -c "$cpus" -t 6:00:00 -D "$NEW" "$PC_ROOT/slurm/step5_place_joint.sbatch"
