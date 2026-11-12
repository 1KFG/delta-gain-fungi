#!/bin/bash
# Which species are MorenoLab JEM1 and JEM2? fastANI against BFD Mucor references.
#   sbatch analysis/mucor_secretome/slurm/identify_jem.sh
# Reads shared folders read-only; writes mucor_secretome/out/jem_ani.tsv.
# Species boundary for fungi is about 95 % ANI.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 1:00:00
#SBATCH -J jem_ani
#SBATCH -o logs/jem_ani_%j.out
#SBATCH -e logs/jem_ani_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes
OUT=$A/out; mkdir -p $OUT
W="${SCRATCH:?}/jem_ani"; mkdir -p $W/ref
module load duckdb/1.1.3 fastani/1.34

# reference assemblies: every BFD Mucor genome of the two candidate species, plus near relatives
cat > $W/q.sql <<'SQL'
.mode csv
.separator "\t"
copy (select ASMID, SPECIES, STRAIN from read_parquet('species.parquet')
      where SPECIES in ('Mucor hiemalis','Mucor circinelloides','Mucor lusitanicus',
                        'Mucor racemosus','Mucor flavus','Mucor piriformis'))
  to '/dev/stdout' (header, delimiter '\t');
SQL
( cd $RUNS/tables && duckdb < $W/q.sql ) > $W/refs_meta.tsv
wc -l $W/refs_meta.tsv

: > $W/refs.txt
while IFS=$'\t' read -r asm sp strain; do
  [ "$asm" = "ASMID" ] && continue
  src=$RUNS/input_clean_genomes/$asm.fa.gz
  [ -f "$src" ] || continue
  zcat "$src" > $W/ref/$asm.fa
  echo "$W/ref/$asm.fa" >> $W/refs.txt
done < $W/refs_meta.tsv
wc -l $W/refs.txt

printf '%s\n%s\n' $MOR/ANI/JEM1.sorted.fasta $MOR/ANI/JEM2.sorted.fasta > $W/queries.txt
fastANI --ql $W/queries.txt --rl $W/refs.txt -t 8 -o $W/ani.tsv
# label the references with their species
awk -F'\t' 'NR==FNR{if(FNR>1){n=$1; s[n]=$2" "$3} next}
  {ref=$2; sub(/.*\//,"",ref); sub(/\.fa$/,"",ref);
   q=$1; sub(/.*\//,"",q); sub(/\.sorted\.fasta$/,"",q);
   print q"\t"ref"\t"(s[ref]?s[ref]:"?")"\t"$3"\t"$4"\t"$5}' \
  $W/refs_meta.tsv $W/ani.tsv | sort -k1,1 -k4,4gr > $W/labelled.tsv
{ echo -e "query\treference_asmid\treference_species\tANI\tfragments_mapped\tfragments_total";
  cat $W/labelled.tsv; } > $OUT/jem_ani.tsv
echo "--- best hits ---"
awk -F'\t' 'NR>1 && !seen[$1]++' $OUT/jem_ani.tsv
echo "--- all >= 90 % ---"
awk -F'\t' 'NR>1 && $4>=90' $OUT/jem_ani.tsv
echo "done $(date +%T)"
