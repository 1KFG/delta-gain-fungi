#!/usr/bin/env bash
# run_genome_classify.sh — sbatch launcher for DeltaGain Stage 1
# (genome_classify.nf): label-free BFD genome classification via mash + skani.
#
# Lives in DeltaGain_Fungi (the project instance), not DeltaGain (the
# pipeline code) — see both repos' READMEs for the split. Submit from
# WITHIN a DeltaGain_Fungi checkout, since results/work/logs land relative
# to $PWD (= Nextflow's launchDir) via conf/profile_genome_classify.config:
#
#   cd DeltaGain_Fungi && sbatch run_genome_classify.sh   [extra nextflow args...]
#
# The Nextflow head process needs to outlive your shell and must not run on
# the login node — it uses a light allocation and submits the real work as
# separate SLURM jobs via the slurm executor. Head partition time limit must
# cover the WHOLE pipeline (it stays alive to submit/monitor every child job).

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name nf-DeltaGain-genome_classify
#SBATCH -o logs/slurm/nf_%j.out
#SBATCH -e logs/slurm/nf_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow
# apptainer FIRST, not singularity: profile_genome_classify.config sets
# apptainer.enabled=true / singularity.enabled=false, and the two modules
# provide different binaries (`apptainer` vs `singularity`) -- module load
# singularity succeeding does NOT put `apptainer` on PATH, so the original
# `singularity || apptainer` fallback order silently short-circuited and
# left every task failing with "env: 'apptainer': No such file or
# directory" (hit for real 2026-09-08, job 28199425).
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

# Container cache — same convention as Fungi_BFD/nextflow. Set one of these
# before submitting if not already exported in your environment.
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}"

# PIPELINE resolves the DeltaGain CODE checkout (nextflow.config, *.nf,
# modules/, bin/, lib/ all live there) -- default assumes the standard
# sibling-folder layout under BFD/. Override if your checkout is elsewhere.
PIPELINE="${PIPELINE:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}"

# -c (not -profile): this project's conf/profile_genome_classify.config is a
# bare config body (not wrapped in a profiles{} block) -- see DeltaGain's
# nextflow.config header comment for why it's loaded this way rather than
# via a named profile baked into the pipeline repo.
nextflow run "${PIPELINE}/genome_classify.nf" \
    -c conf/profile_genome_classify.config \
    -resume \
    "$@"

# Cheap validation of the PIPELINE graph itself (no real data/resources) --
# run from within the DeltaGain checkout, not here:
#   nextflow run . -profile test -stub-run --n_test 2
