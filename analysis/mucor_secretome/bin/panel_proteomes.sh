#!/bin/bash
# Stage the three non-BFD proteomes, once, for whichever annotation job runs first.
#   source analysis/mucor_secretome/bin/panel_proteomes.sh <outdir>
# Writes <outdir>/{JEM1MUCOR,JEM2MUCOR,MPTAK}.faa and .lengths.tsv and does nothing if they are
# already there, so the GPU SignalP job and the CPU predictor job can each call it and the two
# can run in either order or at the same time. Idempotent, and writes only inside <outdir>.
set -euo pipefail
PANEL_OUT="${1:?usage: panel_proteomes.sh <outdir>}"
PANEL_PY=${PY:-/usr/bin/python3.12}
PANEL_A=${A:-analysis/mucor_secretome}
PANEL_MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes/annotation
PANEL_EXT=data_external/marchantia_GCA_039105155.1
mkdir -p "$PANEL_OUT"

# Marchantia polymorpha Tak-1 (GCA_039105155.1, 18,162 proteins). Downloaded data stays in its
# own directory and is only ever read as FASTA.
if [ ! -s $PANEL_EXT/protein.faa ]; then
  mkdir -p $PANEL_EXT
  TMPZ=$(mktemp -d)
  curl -sS --max-time 1800 -o $TMPZ/march.zip \
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_039105155.1/download?include_annotation_type=PROT_FASTA"
  unzip -o -q $TMPZ/march.zip -d $TMPZ/march
  find $TMPZ/march -name "protein.faa" -exec cp {} $PANEL_EXT/protein.faa \;
  rm -rf $TMPZ
fi

if [ ! -s "$PANEL_OUT/JEM1MUCOR.faa" ]; then
  cp $PANEL_MOR/JEM1/predict_results/Mucor_sp._JEM1.proteins.fa "$PANEL_OUT/JEM1MUCOR.faa"
fi
if [ ! -s "$PANEL_OUT/JEM2MUCOR.faa" ]; then
  cp $PANEL_MOR/JEM2/predict_results/Mucor_sp._JEM2.proteins.fa "$PANEL_OUT/JEM2MUCOR.faa"
fi
# Marchantia protein IDs are NCBI accessions with no shared prefix; give them one so the
# LOCUSTAG convention (the prefix before the first '_') holds for every table in this project.
if [ ! -s "$PANEL_OUT/MPTAK.faa" ]; then
  $PANEL_PY -c "
import sys
with open(sys.argv[1]) as fh, open(sys.argv[2], 'w') as fo:
    for line in fh:
        if line.startswith('>'):
            fo.write('>MPTAK_' + line[1:].split()[0].replace('_', '-') + '\n')
        else:
            fo.write(line)
" $PANEL_EXT/protein.faa "$PANEL_OUT/MPTAK.faa"
fi

for LT in JEM1MUCOR JEM2MUCOR MPTAK; do
  [ -s "$PANEL_OUT/$LT.lengths.tsv" ] && continue
  $PANEL_PY -c "
import sys
sys.path.insert(0, '$PANEL_A/bin')
from extract_proteins import read_fasta
with open(sys.argv[2], 'w') as fo:
    for p, s in read_fasta(sys.argv[1]):
        fo.write(f'{p}\t{len(s.rstrip(\"*\"))}\n')
" "$PANEL_OUT/$LT.faa" "$PANEL_OUT/$LT.lengths.tsv"
done
for LT in JEM1MUCOR JEM2MUCOR MPTAK; do
  echo "$LT: $(grep -c '>' "$PANEL_OUT/$LT.faa") proteins"
done
