#!/usr/bin/env bash
# run_protein_novelty.sh — sbatch launcher for DeltaGain Stage 3
# (protein_novelty.nf): stage the comprehensive UniProtKB-fungi reference set
# (FTP bulk download + dat->fasta/xrefs conversion + MMseqs2 linclust
# dereplication) and run the small REST smoke test alongside it.
#
# Lives in DeltaGain_Fungi (the project instance), not DeltaGain (the
# pipeline code) — see both repos' READMEs for the split. Submit from
# WITHIN a DeltaGain_Fungi checkout, since results/work/logs land relative
# to $PWD (= Nextflow's launchDir) via conf/profile_protein_novelty.config:
#
#   cd DeltaGain_Fungi && sbatch run_protein_novelty.sh   [extra nextflow args...]
#
# IMPORTANT if resuming a run started BEFORE this script moved here
# (2026-09-07): Nextflow's -resume cache is keyed to workDir, which is
# `${launchDir}/work/protein_novelty` -- i.e. tied to WHERE `nextflow run`
# was invoked from. A run launched from the old DeltaGain-rooted location
# has its cache under DeltaGain/work/protein_novelty; resuming it from HERE
# would look under DeltaGain_Fungi/work/protein_novelty instead, find
# nothing, and silently restart from scratch (re-downloading the 13GB
# TrEMBL file, re-running the ~45min conversion, etc.) rather than actually
# resuming. Any such in-flight run must keep being resumed from its
# ORIGINAL launch location until it completes; only NEW runs should use
# this script.
#
# This is a genuinely long-running pipeline: uniprot_trembl_fungi.dat.gz is
# ~13 GB (resumable download), and MMSEQS_LINCLUST on the combined ~17.5M-
# sequence set is configured for up to 48h on highmem. The head process must
# survive the whole thing (it stays alive to submit/monitor every child job),
# hence the long time limit below despite the head itself being light on
# cpu/mem.

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name nf-DeltaGain-protein_novelty
#SBATCH -o logs/slurm/nf_%j.out
#SBATCH -e logs/slurm/nf_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow lib/uniprot_fungi

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true

export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}"

# PIPELINE resolves the DeltaGain CODE checkout (nextflow.config, *.nf,
# modules/, bin/, lib/ all live there) -- default assumes the standard
# sibling-folder layout under BFD/. Override if your checkout is elsewhere.
PIPELINE="${PIPELINE:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}"

# -c (not -profile): this project's conf/profile_protein_novelty.config is a
# bare config body (not wrapped in a profiles{} block) -- see DeltaGain's
# nextflow.config header comment for why it's loaded this way rather than
# via a named profile baked into the pipeline repo.
nextflow run "${PIPELINE}/protein_novelty.nf" \
    -c conf/profile_protein_novelty.config \
    -resume \
    "$@"

# Cheap validation of the PIPELINE graph itself (no real data/resources) --
# run from within the DeltaGain checkout, not here:
#   nextflow run . -main-script protein_novelty.nf -profile test -stub-run
