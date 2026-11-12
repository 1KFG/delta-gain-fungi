#!/bin/bash
# Did the Fusarium contamination survive into the assemblies the JEM annotations were built on?
#   sbatch analysis/mucor_secretome/slurm/contam_check.sh
# FCS-GX flagged 2,233 of JEM2's 13,555 spades contigs as Fusarium venenatum and marked them
# EXCLUDE; the cleaned assembly has 2,584 contigs, and the ITS run found a 100 %-identity
# Fusarium graminearum ITS in JEM2.spades. The proteomes are predicted from the cleaned
# assemblies, so contaminant proteins would otherwise be scored as Mucor effectors - this
# checks the cleaned assemblies directly rather than trusting that the purge worked.
#
# Every Ascomycota UNITE reference is blasted against each cleaned assembly, with the spades
# assemblies as positive controls (they should light up) and the BFD genomes as negatives.
# Reads shared folders read-only; writes mucor_secretome/out/contam/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J contam
#SBATCH -o logs/contam_%j.out
#SBATCH -e logs/contam_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-8}
MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes
U=/srv/projects/db/UNITE/current
UNITE=$U/sh_refs_qiime_ver10_dynamic_all_19.02.2025.fasta
UNITETAX=$U/sh_taxonomy_qiime_ver10_dynamic_all_19.02.2025.txt
OUT=$A/out/contam; mkdir -p $OUT
W="${SCRATCH:?}/contam"; mkdir -p $W
module load ncbi-blast/2.16.0+

$PY $A/bin/unite_tax.py subset --fasta $UNITE --tax $UNITETAX --rank p --phylum Ascomycota \
    --out $W/asco.fa

declare -A ASM=( [JEM1_clean]=$MOR/genomes/JEM1.sorted.fasta
                 [JEM2_clean]=$MOR/genomes/JEM2.sorted.fasta
                 [JEM1_masked]=$MOR/genomes/JEM1.masked.fasta
                 [JEM2_masked]=$MOR/genomes/JEM2.masked.fasta
                 [JEM2_spades_POSITIVE]=$MOR/asm/AAFTF/x/JEM2.spades.fasta
                 [JEM1_spades]=$MOR/asm/AAFTF/JEM1.spades.fasta )
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
printf "assembly\tn_hits_ge_97pct\tn_contigs_hit\ttop_species\ttop_pident\n" > $OUT/summary.tsv
for s in JEM1_clean JEM2_clean JEM1_masked JEM2_masked JEM1_spades JEM2_spades_POSITIVE; do
  src=${ASM[$s]}
  [ -s "$src" ] || { echo "skipping $s: no assembly at $src" >&2; continue; }
  makeblastdb -in "$src" -dbtype nucl -out $W/db_$s -title $s > /dev/null
  blastn -query $W/asco.fa -db $W/db_$s -num_threads $CPU -evalue 1e-20 -perc_identity 90 \
      -max_target_seqs 5 -outfmt "$FMT" > $W/$s.hits.tsv
  $PY $A/bin/unite_tax.py rank --blast $W/$s.hits.tsv --tax $UNITETAX --top 10 \
      --out $OUT/rank_$s.tsv
  # a real contaminant shows as near-identical ITS over most of the reference, not as the
  # 85-90 % similarity any fungal rDNA shares with any other
  n=$(awk -F'\t' '$3>=97 && $4>=300' $W/$s.hits.tsv | wc -l)
  nc=$(awk -F'\t' '$3>=97 && $4>=300 {print $2}' $W/$s.hits.tsv | sort -u | wc -l)
  top=$(sed -n '2p' $OUT/rank_$s.tsv | cut -f1,2)
  printf "%s\t%s\t%s\t%s\n" "$s" "$n" "$nc" "$top" >> $OUT/summary.tsv
done
column -t -s $'\t' $OUT/summary.tsv
echo "done $(date +%T)"
