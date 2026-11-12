#!/bin/bash
# Is any JEM protein actually a Fusarium protein? The check that matters for effector calling.
#   sbatch analysis/mucor_secretome/slurm/contam_proteins.sh
# The ITS check (slurm/contam_check.sh) found no Ascomycota rDNA in the cleaned assemblies, but
# it can only see a contaminant contig that carries an rRNA array. A Fusarium contig without one
# would pass it, and its proteins would be scored as Mucor effectors. This compares the two JEM
# proteomes directly against Fusarium proteomes from BFD: a surviving contaminant protein is
# near-identical to its Fusarium counterpart, where a genuine Mucor orthologue is not.
#
# JEM1 is the negative control - its assembly was never Fusarium-contaminated - and the fungal
# candidates of the BFD Mucor reference are a second negative. A real contaminant shows as a
# cluster of >= 90 % identity hits; orthologous conservation across a phylum sits far lower.
# Reads shared folders read-only; writes mucor_secretome/out/contam/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -J contam_prot
#SBATCH -o logs/contam_prot_%j.out
#SBATCH -e logs/contam_prot_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-32}
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
PANELS=$A/out/panels
OUT=$A/out/contam; mkdir -p $OUT
W="${SCRATCH:?}/contam_prot"; mkdir -p $W
module load duckdb/1.1.3 diamond
$PY -W ignore -m unittest discover -s $A/test

# --- Fusarium reference proteomes from BFD
cat > $W/q.sql <<'SQL'
.mode csv
.separator "\t"
copy (select LOCUSTAG, ASMID, SPECIES from read_parquet('species.parquet')
      where GENUS = 'Fusarium')
  to '${W}/fusarium.tsv' (header, delimiter '\t');
SQL
sed -i "s#\${W}#$(readlink -f $W)#" $W/q.sql
( cd $RUNS/tables && duckdb < $W/q.sql )
echo "Fusarium genomes in BFD: $(( $(wc -l < $W/fusarium.tsv) - 1 ))"
# A handful is plenty - this asks "is this protein Fusarium", not "which Fusarium". Written with
# awk rather than `tail | head`, because head exits early, tail takes SIGPIPE, and under
# `set -o pipefail` that non-zero status ends the job (it did: 32050463).
awk -F'\t' 'NR>1 && n<6 {print $1; n++}' $W/fusarium.tsv | sort -u > $W/fus_locustags.txt
cat $W/fus_locustags.txt
$PY $A/bin/extract_proteins.py --seqdb-index $RUNS/seqdb/genomes.tsv \
    --locustags $W/fus_locustags.txt --threads $CPU --out $W/fusarium.faa
grep -c '>' $W/fusarium.faa

diamond makedb --in $W/fusarium.faa -d $W/fus --threads $CPU --quiet
for S in JEM1MUCOR JEM2MUCOR; do
  diamond blastp -q $PANELS/faa/$S.faa -d $W/fus -o $W/$S.tsv --threads $CPU --quiet \
      --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
      --max-target-seqs 1 --evalue 1e-5 --query-cover 80
  n90=$(awk -F'\t' '$3>=90 && $4/$5>=0.8' $W/$S.tsv | wc -l)
  n95=$(awk -F'\t' '$3>=95 && $4/$5>=0.8' $W/$S.tsv | wc -l)
  tot=$(grep -c '>' $PANELS/faa/$S.faa)
  echo "$S: $tot proteins, $(wc -l < $W/$S.tsv) with a Fusarium hit, $n90 at >=90% id, $n95 at >=95%"
  cp $W/$S.tsv $OUT/${S}_vs_fusarium.tsv
done

echo "=== identity distribution of best Fusarium hits ==="
for S in JEM1MUCOR JEM2MUCOR; do
  printf "%s " $S
  awk -F'\t' '{b=int($3/10)*10; c[b]++} END{for(i=30;i<=90;i+=10) printf "%d-%d:%d ", i, i+9, c[i]+0; print ""}' $W/$S.tsv
done
echo "done $(date +%T)"
