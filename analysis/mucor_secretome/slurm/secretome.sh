#!/bin/bash
# Mucoromycota secretome and the Mucor hiemalis screen set.
#   sbatch analysis/mucor_secretome/slurm/secretome.sh
# Reads the shared BFD tables, seqdb and input assemblies read-only; writes mucor_secretome/out/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 6:00:00
#SBATCH -J mucor_secr
#SBATCH -o logs/mucor_secr_%j.out
#SBATCH -e logs/mucor_secr_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
FOCUS=F00CDB58                       # Mucor hiemalis, GCA_977110905.1_gzMucHiem1
OUT=$A/out; mkdir -p $OUT
W="${SCRATCH:?}/mucor_secr"; mkdir -p $W
module load duckdb/1.1.3
module load MMseqs2/17-b804f

$PY -W ignore -m unittest discover -s $A/test
$PY -W ignore -m unittest discover -s $UPO/bfd_findings/test

# 1. annotation tables for every Mucoromycota genome
sed -e "s#\${OUT}#$W#" $A/sql/export.sql > $W/q.sql
sed -i "1i SET memory_limit='52GB'; SET threads=16; SET temp_directory='$W/duck_tmp';" $W/q.sql
( cd $RUNS/tables && duckdb < $W/q.sql )
wc -l $W/*.tsv

# 2. sequences for the candidate pool (every protein with a signal peptide at p >= 0.5)
tail -n +2 $W/genomes.tsv | cut -f1 | sort -u > $W/locustags.txt
tail -n +2 $W/signalp.tsv | cut -f1 | sort -u > $W/ids.txt
$PY $A/bin/extract_proteins.py --seqdb-index $RUNS/seqdb/genomes.tsv \
    --locustags $W/locustags.txt --ids $W/ids.txt --threads 16 --out $W/candidates.faa

# 3. protein families over the candidate pool (50 % identity, 80 % coverage of the shorter)
mmseqs easy-cluster $W/candidates.faa $W/clu $W/mmtmp \
    --min-seq-id 0.5 -c 0.8 --cov-mode 5 --threads 16 -v 1
wc -l $W/clu_cluster.tsv

# 4. real 5' ends: BFD's start/stop flags are wrong in these genomes before Fungi_BFD fb1d7cd
#    (enzymes_UPO/bfd_findings/REPORT.md), and a model without a true start loses its signal peptide
$PY $UPO/bfd_findings/bin/check_codons.py --cds $W/codon_cds.tsv --models $W/codon_models.tsv \
    --genomes-dir $RUNS/input_clean_genomes --proteins $W/candidates.faa --threads 16 \
    --out $OUT/codon_check.tsv

# 5. classify, tier and write the screen set
$PY $A/bin/build_secretome.py --dir $W --codons $OUT/codon_check.tsv \
    --clusters $W/clu_cluster.tsv --fasta $W/candidates.faa --focus $FOCUS --out $OUT
cp $W/genomes.tsv $OUT/mucoromycota_genomes.tsv
echo "done $(date +%T)"
