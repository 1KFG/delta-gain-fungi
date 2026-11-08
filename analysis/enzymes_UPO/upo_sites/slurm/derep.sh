#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -J upo_derep
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/derep_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/derep_%j.err
set -euo pipefail
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
OUT=$PROJ/data/bfd_full/novelty
mkdir -p "$OUT"
export TMPDIR="${SCRATCH:?SCRATCH not set}/upo_derep.$SLURM_JOB_ID"; mkdir -p "$TMPDIR"
module load mmseqs2/13-45111

# Dereplicate the CATALYTIC CORE set, not all 97,520. Proteins lacking the UPO catalytic
# residues are not what the novelty question is about. See G1_BFD_SUMMARY.md section 6.
IN=$PROJ/data/bfd_full/core_pf01328.faa
for ID in 0.90 0.70; do
  tag=${ID/./}
  mmseqs easy-cluster "$IN" "$TMPDIR/clu$tag" "$TMPDIR/tmp$tag" \
      --min-seq-id "$ID" -c 0.8 --cov-mode 0 --threads "${SLURM_CPUS_PER_TASK:-16}" >/dev/null
  cp "$TMPDIR/clu${tag}_rep_seq.fasta" "$OUT/reps_${tag}.faa"
  cp "$TMPDIR/clu${tag}_cluster.tsv"   "$OUT/clusters_${tag}.tsv"
  echo "id=$ID representatives: $(grep -c '^>' "$OUT/reps_${tag}.faa")"
done
zstd -qf --rm "$OUT"/clusters_*.tsv 2>/dev/null || true
rm -rf "$TMPDIR"
