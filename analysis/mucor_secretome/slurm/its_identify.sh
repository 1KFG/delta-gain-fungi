#!/bin/bash
# What species are JEM1 and JEM2? rDNA/ITS against UNITE, with the BFD genomes as controls.
#   sbatch analysis/mucor_secretome/slurm/its_identify.sh
# fastANI placed JEM2 in M. circinelloides (99.8 %) but left JEM1 with no relative above 80 %,
# and the only M. hiemalis genome in BFD is a single assembly. ITS reaches the many described
# Mucor species that have no sequenced genome, and tests that single reference as well.
#
# The rDNA array is collapsed or broken in these assemblies, so ITSx recovers little or nothing
# (job 32019432 found 0 regions; the MorenoLab ITSx run's JEM1 "full" sequence is SSU, reported
# with a negative ITS length). The primary method here therefore does not depend on ITSx: every
# Mucoromycota UNITE reference is blasted against each assembly, and the best hit per species is
# ranked. ITSx is still run, with the --nhmmer settings that work on genome-length contigs, and
# whatever it extracts is blasted against the whole of UNITE as a cross-check.
# Reads shared folders read-only; writes mucor_secretome/out/its/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -J its_id
#SBATCH -o logs/its_id_%j.out
#SBATCH -e logs/its_id_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-8}
MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes
RUNS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
U=/srv/projects/db/UNITE/current
UNITE=$U/sh_refs_qiime_ver10_dynamic_all_19.02.2025.fasta
UNITETAX=$U/sh_taxonomy_qiime_ver10_dynamic_all_19.02.2025.txt
OUT=$A/out/its; mkdir -p $OUT
W="${SCRATCH:?}/its_id"; mkdir -p $W
module load ncbi-blast/2.16.0+ ITSx/1.1.3

# --- assemblies: the two strains by both assemblers, and the BFD references as controls
#     (the single M. hiemalis, plus M. circinelloides, which JEM2 matched at 99.8 % ANI)
# JEM2's AAFTF run landed in a stray asm/AAFTF/x/ subdirectory, which is why ITS/ holds four
# broken JEM2 symlinks; $MOR/genomes/ is the stable path for the sorted assemblies (it is what
# the lab's own ANI run reads). The spades assemblies are taken wherever they are, and a
# missing one is skipped rather than failing the job.
declare -A ASM=( [BFD_M_hiemalis]=$RUNS/input_clean_genomes/GCA_977110905.1_gzMucHiem1.fa.gz
                 [BFD_M_circinelloides]=$RUNS/input_clean_genomes/GCA_973907015.1_gzMucCirc1.fa.gz
                 [JEM1_sorted]=$MOR/genomes/JEM1.sorted.fasta
                 [JEM2_sorted]=$MOR/genomes/JEM2.sorted.fasta
                 [JEM1_spades]=$MOR/asm/AAFTF/JEM1.spades.fasta
                 [JEM2_spades]=$MOR/asm/AAFTF/x/JEM2.spades.fasta )
NAMES=""
for s in BFD_M_hiemalis BFD_M_circinelloides JEM1_sorted JEM2_sorted JEM1_spades JEM2_spades; do
  src=${ASM[$s]}
  if [ ! -s "$src" ]; then echo "skipping $s: no assembly at $src" >&2; continue; fi
  case $src in *.gz) zcat "$src" > $W/$s.fa ;; *) cp "$src" $W/$s.fa ;; esac
  NAMES="$NAMES $s"
done
echo "assemblies:$NAMES"

# --- UNITE restricted to Mucoromycota. The lineage is in the taxonomy file, not the FASTA
#     header, so this goes through the accession map (see bin/unite_tax.py).
$PY $A/bin/unite_tax.py subset --fasta $UNITE --tax $UNITETAX --phylum Mucoromycota \
    --out $W/mucoro.fa
grep -c '>' $W/mucoro.fa

# --- primary: every Mucoromycota reference against each assembly, best hit per species
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
for s in $NAMES; do
  makeblastdb -in $W/$s.fa -dbtype nucl -out $W/db_$s -title $s > /dev/null
  blastn -query $W/mucoro.fa -db $W/db_$s -num_threads $CPU -evalue 1e-20 \
      -max_target_seqs 5 -outfmt "$FMT" > $W/$s.hits.tsv
  $PY $A/bin/unite_tax.py rank --blast $W/$s.hits.tsv --tax $UNITETAX --top 25 \
      --out $OUT/rank_$s.tsv
  echo "=== $s: top 8 species by rDNA ==="
  head -9 $OUT/rank_$s.tsv | cut -f1-6 | column -t -s $'\t'
done

# --- cross-check: ITSx with the settings that work on genome-length contigs, then the
#     extracted regions against the whole of UNITE (not only Mucoromycota)
: > $W/query_its.fa
for s in $NAMES; do
  ITSx -i $W/$s.fa -o $W/$s.itsx --cpu $CPU --nhmmer T --allow_reorder T --save_regions all \
      -t fungi --preserve T > $W/$s.itsx.log 2>&1 || true
  cp -f $W/$s.itsx.summary.txt $OUT/itsx_$s.summary.txt 2>/dev/null || true
  for reg in full ITS1 ITS2; do
    f=$W/$s.itsx.$reg.fasta
    [ -s "$f" ] || continue
    $PY -c "
import sys; sys.path.insert(0, '$A/bin')
from extract_proteins import read_fasta
for i, (h, seq) in enumerate(read_fasta('$f')):
    if len(seq) >= 100:
        print(f'>{\"$s\"}|{\"$reg\"}|{i}\n{seq}')
" >> $W/query_its.fa
  done
  echo "$s: $(grep -c '>' $W/query_its.fa) ITSx regions kept so far"
done
cp -f $W/query_its.fa $OUT/query_its.fa

if [ -s $W/query_its.fa ]; then
  makeblastdb -in $UNITE -dbtype nucl -out $W/unite -title unite > /dev/null
  blastn -query $W/query_its.fa -db $W/unite -num_threads $CPU -max_target_seqs 20 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
      > $W/unite.blast.tsv
  $PY - "$W/unite.blast.tsv" "$UNITETAX" "$OUT/itsx_vs_unite.tsv" <<'PYEOF'
import sys
sys.path.insert(0, 'analysis/mucor_secretome/bin')
from unite_tax import parse_taxonomy, rank_of, read_blast
blast, taxf, out = sys.argv[1:4]
with open(taxf) as fh:
    tax = parse_taxonomy(fh)
with open(blast) as fh, open(out, 'w') as fo:
    fo.write("query\tspecies\tpident\taln_len\tbitscore\tsubject\n")
    seen = {}
    for r in read_blast(fh):
        q = r['qseqid']
        if seen.get(q, 0) >= 5:
            continue
        seen[q] = seen.get(q, 0) + 1
        lin = tax.get(r['sseqid'], '')
        sp = rank_of(lin, 's') or rank_of(lin, 'g') or r['sseqid']
        fo.write(f"{q}\t{sp}\t{r['pident']}\t{r['length']}\t{r['bitscore']}\t{r['sseqid']}\n")
PYEOF
  echo "=== ITSx regions vs all of UNITE (top 5 each) ==="
  column -t -s $'\t' $OUT/itsx_vs_unite.tsv
else
  echo "ITSx extracted nothing above 100 bp from any assembly"
fi
echo "done $(date +%T)"
