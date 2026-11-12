#!/bin/bash
# Stage the non-BFD proteomes, once, for whichever annotation job runs first.
#   bash analysis/mucor_secretome/bin/panel_proteomes.sh <outdir> [LOCUSTAG ...]
# Writes <outdir>/<LOCUSTAG>.faa and .lengths.tsv and does nothing if they are already there, so
# the GPU SignalP job and the CPU predictor job can each call it and the two can run in either
# order or at the same time. Idempotent, and writes only inside <outdir> and data_external/.
#
# Panels:
#   JEM1MUCOR, JEM2MUCOR  the MorenoLab Mucor strains isolated from Marchantia
#   MPTAK                 Marchantia polymorpha Tak-1, GCA_039105155.1
#   MPALE                 Marchantia paleacea, GCA_014180765.2
# Both liverworts are staged because the host species of the isolates is not recorded, and the
# two differ in the way that matters: M. paleacea forms arbuscular mycorrhizal symbiosis, and
# M. polymorpha has lost it. Annotating both turns an unknown into a comparison.
set -euo pipefail
PANEL_OUT="${1:?usage: panel_proteomes.sh <outdir> [LOCUSTAG ...]}"; shift || true
PANEL_WANT="${*:-JEM1MUCOR JEM2MUCOR MPTAK MPALE}"
PANEL_PY=${PY:-/usr/bin/python3.12}
PANEL_A=${A:-analysis/mucor_secretome}
PANEL_MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes/annotation
mkdir -p "$PANEL_OUT"

wants() { case " $PANEL_WANT " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

# Downloaded data stays in its own directory and is only ever read as FASTA.
fetch_ncbi() {   # fetch_ncbi <accession> <destdir>
  local acc="$1" dest="$2" tmpz
  [ -s "$dest/protein.faa" ] && return 0
  mkdir -p "$dest"
  tmpz=$(mktemp -d)
  curl -sS --max-time 1800 -o "$tmpz/g.zip" \
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${acc}/download?include_annotation_type=PROT_FASTA"
  unzip -o -q "$tmpz/g.zip" -d "$tmpz/g"
  find "$tmpz/g" -name "protein.faa" -exec cp {} "$dest/protein.faa" \;
  rm -rf "$tmpz"
  [ -s "$dest/protein.faa" ] || { echo "no protein.faa for $acc" >&2; return 1; }
}

# Plant protein IDs are NCBI accessions with no shared prefix; give them one so the LOCUSTAG
# convention (the prefix before the first '_') holds for every table in this project.
prefix_ids() {   # prefix_ids <in.faa> <out.faa> <LOCUSTAG>
  $PANEL_PY -c "
import sys
lt = sys.argv[3]
with open(sys.argv[1]) as fh, open(sys.argv[2], 'w') as fo:
    for line in fh:
        if line.startswith('>'):
            fo.write('>' + lt + '_' + line[1:].split()[0].replace('_', '-') + '\n')
        else:
            fo.write(line)
" "$1" "$2" "$3"
}

if wants JEM1MUCOR && [ ! -s "$PANEL_OUT/JEM1MUCOR.faa" ]; then
  cp $PANEL_MOR/JEM1/predict_results/Mucor_sp._JEM1.proteins.fa "$PANEL_OUT/JEM1MUCOR.faa"
fi
if wants JEM2MUCOR && [ ! -s "$PANEL_OUT/JEM2MUCOR.faa" ]; then
  cp $PANEL_MOR/JEM2/predict_results/Mucor_sp._JEM2.proteins.fa "$PANEL_OUT/JEM2MUCOR.faa"
fi
if wants MPTAK && [ ! -s "$PANEL_OUT/MPTAK.faa" ]; then
  fetch_ncbi GCA_039105155.1 data_external/marchantia_GCA_039105155.1
  prefix_ids data_external/marchantia_GCA_039105155.1/protein.faa "$PANEL_OUT/MPTAK.faa" MPTAK
fi
if wants MPALE && [ ! -s "$PANEL_OUT/MPALE.faa" ]; then
  fetch_ncbi GCA_014180765.2 data_external/marchantia_paleacea_GCA_014180765.2
  prefix_ids data_external/marchantia_paleacea_GCA_014180765.2/protein.faa \
      "$PANEL_OUT/MPALE.faa" MPALE
fi

for LT in $PANEL_WANT; do
  [ -s "$PANEL_OUT/$LT.faa" ] || { echo "missing $PANEL_OUT/$LT.faa" >&2; exit 1; }
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
for LT in $PANEL_WANT; do
  echo "$LT: $(grep -c '>' "$PANEL_OUT/$LT.faa") proteins"
done
