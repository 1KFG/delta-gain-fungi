#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 8 --mem 24gb --time 2:00:00
#SBATCH --out /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/bgc_subset/logs/antismash_%x.%A.log
# One genome through antiSMASH 8.0.4 in the BFD apptainer image, then reduce to the subset.
# Mirrors Fungi_BFD/nextflow/modules/BFD/ANTISMASH/main.nf (gff3 + genome fasta, --genefinding-gff3),
# plus --cb-knownclusters (MIBiG) which that module does not turn on.
# --pfam2go, --fullhmmer, --clusterhmmer dropped (not used by the subset). --cb-general kept as the explicit flag set
# was asked for; the subset's blast_top table needs it. All work in $SCRATCH; only the subset is copied to /bigdata.
# Set KEEP_RAW=1 to also keep the gzipped raw JSON (about 20 MB) for a timing/comparison run.
# Usage: sbatch -J <name> run_antismash_one.sbatch <genome_folder_name> <short_id>
set -euo pipefail
NAME=${1:?genome folder name}; ID=${2:?short id}
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation
SIF=/bigdata/stajichlab/shared/singularity/antismash-standalone-8.0.4.sif
DB=/bigdata/stajichlab/shared/lib/antismash_db
HERE=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/bgc_subset
OUT=${BGC_OUT:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi_wt_realdata/analysis/bgc_subset/results}
mkdir -p $OUT
module load apptainer
W=${SCRATCH:?SCRATCH not set}/as_$ID; mkdir -p $W; cd $W
export TMPDIR=$W/tmp; mkdir -p $TMPDIR
zcat $RUNS/$NAME/predict_results/$NAME.gff3.gz > $ID.gff3
zcat $RUNS/$NAME/predict_results/$NAME.scaffolds.fa.gz > $ID.fa
OUTD=$W/out
t0=$(date +%s)
apptainer exec --bind $DB:$DB,$W:$W,$TMPDIR:$TMPDIR $SIF \
  antismash --taxon fungi --databases $DB --output-dir $OUTD \
    --genefinding-gff3 $ID.gff3 \
    --cb-general --cb-knownclusters \
    -c ${SLURM_CPUS_ON_NODE:-8} $ID.fa
t1=$(date +%s)
J=$(ls $OUTD/*.json | head -1)
du -sk $OUTD | awk '{print $1" KB raw antismash dir"}'
/usr/bin/python3.12 $HERE/extract_antismash_subset.py $J $OUT $ID
[ "${KEEP_RAW:-0}" = 1 ] && gzip -9 -c $J > $OUT/$ID.antismash.json.gz
echo -e "$ID\t$NAME\t$((t1-t0))\t${SLURM_CPUS_ON_NODE:-8}" >> $OUT/timings.tsv
rm -rf $W
