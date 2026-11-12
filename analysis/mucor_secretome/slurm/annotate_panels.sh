#!/bin/bash
# Secretome annotation for the genomes that are not in BFD: the two MorenoLab Mucor strains
# (JEM1, JEM2 - isolated from Marchantia) and the Marchantia polymorpha reference proteome.
#   sbatch analysis/mucor_secretome/slurm/annotate_panels.sh
# Runs the same predictors BFD uses, so the three proteomes are classified by one set of rules.
# Reads shared folders read-only; writes mucor_secretome/out/panels/ and data_external/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 24
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -J annot_panels
#SBATCH -o logs/annot_panels_%j.out
#SBATCH -e logs/annot_panels_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes/annotation
DB=/bigdata/stajichlab/shared/lib/funannotate_db
EXT=data_external
OUT=$A/out/panels; mkdir -p $OUT $EXT
W="${SCRATCH:?}/annot_panels"; mkdir -p $W
module load signalp/6.0i tmhmm/2.0c hmmer/3.4
$PY -W ignore -m unittest discover -s $A/test

# --- Marchantia polymorpha Tak-1 reference (GCA_039105155.1, NIG annotation, 18,162 proteins).
# Downloaded data goes in its own empty directory and is only ever read as FASTA.
MARCH=$EXT/marchantia_GCA_039105155.1
if [ ! -s $MARCH/protein.faa ]; then
  mkdir -p $MARCH
  curl -sS --max-time 1800 -o $W/march.zip \
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_039105155.1/download?include_annotation_type=PROT_FASTA"
  unzip -o -q $W/march.zip -d $W/march
  find $W/march -name "protein.faa" -exec cp {} $MARCH/protein.faa \;
fi
grep -c '>' $MARCH/protein.faa

# --- proteomes: locustag, fasta, (gff3, scaffolds) for the fungal ones
cp $MOR/JEM1/predict_results/Mucor_sp._JEM1.proteins.fa $W/JEM1MUCOR.faa
cp $MOR/JEM2/predict_results/Mucor_sp._JEM2.proteins.fa $W/JEM2MUCOR.faa
# Marchantia protein IDs are NCBI accessions with no shared prefix; give them one so the
# LOCUSTAG convention (prefix before the first '_') holds for every table in this project.
$PY -c "
import sys
inp, out = sys.argv[1], sys.argv[2]
with open(inp) as fh, open(out, 'w') as fo:
    for line in fh:
        if line.startswith('>'):
            fo.write('>MPTAK_' + line[1:].split()[0].replace('_', '-') + '\n')
        else:
            fo.write(line)
" $MARCH/protein.faa $W/MPTAK.faa
grep -c '>' $W/MPTAK.faa

for LT in JEM1MUCOR JEM2MUCOR MPTAK; do
  F=$W/$LT.faa
  echo "=== $LT $(grep -c '>' $F) proteins $(date +%T)"
  $PY -c "
import sys
sys.path.insert(0, '$A/bin')
from extract_proteins import read_fasta
with open('$W/$LT.lengths.tsv', 'w') as fo:
    for p, s in read_fasta('$F'):
        fo.write(f'{p}\t{len(s.rstrip(\"*\"))}\n')
"
  signalp6 --fastafile $F --output_dir $W/sp_$LT --format none --organism eukarya \
      --mode fast --torch_num_threads 24 --write_procs 8
  tmhmm --short $F > $W/$LT.tmhmm.txt 2> $W/$LT.tmhmm.err || true
  hmmsearch --cpu 24 --cut_ga --domtblout $W/$LT.pfam.domtbl $DB/Pfam-A.hmm $F > /dev/null
  DBC=""
  if [ "$LT" != "MPTAK" ]; then
    hmmsearch --cpu 24 -E 1e-5 --domtblout $W/$LT.dbcan.domtbl $DB/dbCAN.hmm $F > /dev/null
    DBC="--dbcan $W/$LT.dbcan.domtbl"
  fi
  $PY $A/bin/annotate_to_tables.py --locustag $LT \
      --signalp $W/sp_$LT/prediction_results.txt --tmhmm $W/$LT.tmhmm.txt \
      --pfam $W/$LT.pfam.domtbl $DBC --lengths $W/$LT.lengths.tsv --out $OUT
done

# --- real 5' ends for the two Mucor strains, from their own assemblies
mkdir -p $W/asm
for S in JEM1 JEM2; do
  LT=${S}MUCOR
  gzip -c $MOR/$S/predict_results/Mucor_sp._$S.scaffolds.fa > $W/asm/$S.fa.gz
  $PY $A/bin/gff3_to_cds.py --gff3 $MOR/$S/predict_results/Mucor_sp._$S.gff3 \
      --asmid $S --set $LT --out-cds $W/$S.cds.tsv --out-models $W/$S.models.tsv
done
head -1 $W/JEM1.models.tsv > $W/jem.models.tsv
tail -q -n +2 $W/JEM1.models.tsv $W/JEM2.models.tsv >> $W/jem.models.tsv
head -1 $W/JEM1.cds.tsv > $W/jem.cds.tsv
tail -q -n +2 $W/JEM1.cds.tsv $W/JEM2.cds.tsv >> $W/jem.cds.tsv
cat $W/JEM1MUCOR.faa $W/JEM2MUCOR.faa > $W/jem.faa
$PY $UPO/bfd_findings/bin/check_codons.py --cds $W/jem.cds.tsv --models $W/jem.models.tsv \
    --genomes-dir $W/asm --proteins $W/jem.faa --threads 24 --out $OUT/jem_codon_check.tsv

cp $W/JEM1MUCOR.faa $W/JEM2MUCOR.faa $W/MPTAK.faa $OUT/
wc -l $OUT/*.tsv
echo "done $(date +%T)"
