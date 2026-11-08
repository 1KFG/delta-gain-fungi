#!/bin/bash
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -J upo_streme
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/streme_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO/logs/streme_%j.err
# G3 Phase A step 2: within-lineage clade contrasts, de novo motif discovery + known-TF naming.
#
# PRIMARY DESIGN: a clade's motifs are discovered separately inside each lineage where that clade
# is well represented. A motif that REPLICATES ACROSS INDEPENDENT LINEAGES is evidence of a clade
# motif. A motif found in only one lineage is a lineage motif. This is the control that would have
# caught the pocket-versus-clade confound found earlier in this project.
set -euo pipefail
PROJ=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/enzymes_UPO
G3=$PROJ/data/g3
SEQ=$G3/seqs
OUT=$G3/streme
MAP=$G3/rep_clade_lineage.tsv
JASPAR=$G3/motifdb/jaspar_fungi.meme
THREADS=${SLURM_CPUS_PER_TASK:-8}

# Set MEME paths directly. `module load meme/5.5.5` also loads miniconda3, whose
# `conda activate` fails in a non-interactive batch shell (job 31865837 died on this).
# The MEME binaries themselves do not need conda.
MEMEROOT=/opt/linux/rocky/8.x/x86_64/pkgs/meme/5.5.5
export PATH="$MEMEROOT/bin:$MEMEROOT/libexec/meme-5.5.5:$PATH"
export PERL5LIB="$MEMEROOT/lib/meme-5.5.5/perl${PERL5LIB:+:$PERL5LIB}"
export PYTHONPATH="$MEMEROOT/lib/meme-5.5.5/python${PYTHONPATH:+:$PYTHONPATH}"
streme --version >/dev/null || { echo "streme not runnable"; exit 1; }
mkdir -p "$OUT"

# Contrasts with enough representatives on both sides (counts measured 2026-11-08).
# format: LINEAGE:CLADE
CONTRASTS="
Agaricomycetes:I
Agaricomycetes:III
Agaricomycetes:V
Agaricomycetes:II
Sordariomycetes:V
Sordariomycetes:IV
Sordariomycetes:III
Dothideomycetes:V
Dothideomycetes:IV
Dothideomycetes:II
Eurotiomycetes:V
Eurotiomycetes:IV
Leotiomycetes:V
Leotiomycetes:IV
Lecanoromycetes:V
Lecanoromycetes:III
Lecanoromycetes:I
"

for c in $CONTRASTS; do
  LIN=${c%%:*}; CLADE=${c##*:}
  tag="${LIN}_clade${CLADE}"
  fg="$OUT/$tag.fg.fa"; bg="$OUT/$tag.bg.fa"

  # foreground: this clade, this lineage
  /usr/bin/python3.12 "$PROJ/upo_reg/bin/subset_fasta.py" \
      --fasta "$SEQ/target.fa" --map "$MAP" --lineage "$LIN" --clade "$CLADE" --out "$fg"
  # background: genome-matched non-UPO promoters from the SAME lineage
  /usr/bin/python3.12 "$PROJ/upo_reg/bin/subset_fasta.py" \
      --fasta "$SEQ/background.fa" --map "$MAP" --lineage "$LIN" --any-clade --out "$bg"

  nfg=$(grep -c '^>' "$fg" || true); nbg=$(grep -c '^>' "$bg" || true)
  if [ "$nfg" -lt 50 ] || [ "$nbg" -lt 50 ]; then
    echo "SKIP $tag (fg=$nfg bg=$nbg, under 50)"; rm -f "$fg" "$bg"; continue
  fi
  echo "=== $tag fg=$nfg bg=$nbg $(date +%T) ==="

  streme --p "$fg" --n "$bg" --dna --oc "$OUT/$tag" --nmotifs 8 --verbosity 1 || \
    { echo "streme failed for $tag"; continue; }

  # known-TF naming. JASPAR CORE fungi is 178 motifs and heavily Saccharomyces-biased,
  # so an unnamed motif means "no close match in a thin database", NOT "novel".
  if [ -s "$OUT/$tag/streme.xml" ]; then
    tomtom -no-ssc -oc "$OUT/$tag/tomtom" -verbosity 1 -min-overlap 5 \
           -dist pearson -thresh 0.05 "$OUT/$tag/streme.txt" "$JASPAR" || true
  fi
done

# Second control: dinucleotide shuffle. This build has no fasta-shuffle-letters, but STREME
# generates its own order-2 shuffled control when --n is omitted, which is the same control.
echo "=== shuffle control $(date +%T) ==="
streme --p "$SEQ/target.fa" --dna --order 2 \
       --oc "$OUT/pooled_vs_shuffle" --nmotifs 10 --verbosity 1 || true

# Pooled real-background run, for comparison against the shuffle control.
echo "=== pooled vs genome-matched background $(date +%T) ==="
streme --p "$SEQ/target.fa" --n "$SEQ/background.fa" --dna \
       --oc "$OUT/pooled_vs_background" --nmotifs 10 --verbosity 1 || true
echo "done $(date +%T)"
