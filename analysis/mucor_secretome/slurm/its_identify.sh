#!/bin/bash
# What species are JEM1 and JEM2? ITS against UNITE, with the BFD reference genomes as controls.
#   sbatch analysis/mucor_secretome/slurm/its_identify.sh
# fastANI placed JEM2 in M. circinelloides (99.8 %) but left JEM1 with no relative above 80 %,
# and the only M. hiemalis genome in BFD is a single assembly. ITS reaches the many described
# Mucor species that have no sequenced genome, and tests that single reference as well.
# Reads shared folders read-only; writes mucor_secretome/out/its/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -J its_id
#SBATCH -o logs/its_id_%j.out
#SBATCH -e logs/its_id_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
UNITE=/srv/projects/db/UNITE/current/sh_refs_qiime_ver10_dynamic_all_19.02.2025.fasta
UNITETAX=/srv/projects/db/UNITE/current/sh_taxonomy_qiime_ver10_dynamic_all_19.02.2025.txt
OUT=$A/out/its; mkdir -p $OUT
W="${SCRATCH:?}/its_id"; mkdir -p $W
module load ncbi-blast/2.16.0+ ITSx/1.1.3

# --- control genomes: the BFD assemblies this identification leans on
#     (the single M. hiemalis, plus M. circinelloides, which JEM2 matched at 99.8 % ANI)
declare -A CTRL=( [BFD_M_hiemalis]=GCA_977110905.1_gzMucHiem1
                  [BFD_M_circinelloides]=GCA_973907015.1_gzMucCirc1 )
for name in "${!CTRL[@]}"; do
  asm=${CTRL[$name]}
  zcat $RUNS/input_clean_genomes/$asm.fa.gz > $W/$name.fa
done
cp $MOR/ANI/JEM1.sorted.fasta $W/JEM1.fa
cp $MOR/ANI/JEM2.sorted.fasta $W/JEM2.fa

# --- ITS from each assembly, by the same method for all four
for s in JEM1 JEM2 BFD_M_hiemalis BFD_M_circinelloides; do
  ITSx -i $W/$s.fa -o $W/$s.itsx --cpu 8 --save_regions all -t F --preserve T > $W/$s.itsx.log 2>&1 || true
  if [ -s $W/$s.itsx.full.fasta ]; then
    $PY -c "
import sys
sys.path.insert(0, '$A/bin')
from extract_proteins import read_fasta
recs = list(read_fasta('$W/$s.itsx.full.fasta'))
best = max(recs, key=lambda r: len(r[1])) if recs else None
if best:
    print('>' + '$s' + '\n' + best[1])
" > $W/$s.its.fa
  else
    : > $W/$s.its.fa
  fi
  echo "$s ITS: $(grep -v '>' $W/$s.its.fa | tr -d '\n' | wc -c) bp (ITSx found $(grep -c '>' $W/$s.itsx.full.fasta 2>/dev/null || echo 0))"
done
# the MorenoLab ITSx run already produced a 910 bp JEM1 ITS; keep it if ours is shorter
if [ "$(grep -v '>' $W/JEM1.its.fa | tr -d '\n' | wc -c)" -lt 400 ]; then
  echo "using the MorenoLab ITSx result for JEM1"
  sed 's/^>.*/>JEM1/' $MOR/ITS/JEM1.sorted.full.fasta > $W/JEM1.its.fa
fi
cat $W/JEM1.its.fa $W/JEM2.its.fa $W/BFD_M_hiemalis.its.fa $W/BFD_M_circinelloides.its.fa \
  > $OUT/query_its.fa
grep '>' $OUT/query_its.fa

# --- UNITE, restricted to Mucoromycota for a readable ranking, and whole-database for a check
makeblastdb -in $UNITE -dbtype nucl -out $W/unite -title unite > /dev/null
$PY - "$UNITE" "$W/mucoro.fa" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
keep, n = False, 0
with open(src) as fh, open(dst, 'w') as fo:
    for line in fh:
        if line.startswith('>'):
            keep = 'Mucoromycota' in line
            n += keep
        if keep:
            fo.write(line)
print(f'{n} Mucoromycota reference sequences', file=sys.stderr)
PYEOF
makeblastdb -in $W/mucoro.fa -dbtype nucl -out $W/mucoro -title mucoro > /dev/null

for db in unite mucoro; do
  blastn -query $OUT/query_its.fa -db $W/$db -num_threads 8 -max_target_seqs 50 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
      > $W/$db.blast.tsv
  { printf "query\tsubject\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\n";
    cat $W/$db.blast.tsv; } > $OUT/its_vs_${db}.tsv
  echo "=== $db: best 5 per query ==="
  awk -F'\t' 'NR>1{if(c[$1]++<5) print $1"\t"$3"% cov"$13"\t"$2}' $OUT/its_vs_${db}.tsv
done
echo "done $(date +%T)"
