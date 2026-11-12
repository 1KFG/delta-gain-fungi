#!/bin/bash
# The CPU half of the secretome annotation for the genomes that are not in BFD: the two
# MorenoLab Mucor strains (JEM1, JEM2 - isolated from Marchantia) and the Marchantia polymorpha
# reference proteome.
#   sbatch analysis/mucor_secretome/slurm/annotate_panels.sh
# SignalP runs on a GPU in slurm/signalp_gpu.sh; the two jobs write different tables into the
# same directory and may run in either order.
#
# Each predictor here is invoked as the matching Fungi_BFD/nextflow/modules/BFD/* process
# invokes it, so these three proteomes are annotated by the same tools, versions and flags as
# every BFD genome and build_secretome.py can classify them by one set of rules:
#   TMHMM    tmhmm --short --noplot          (RUN_TMHMM)
#   PredGPI  predgpi.py -f X -m gff3         (RUN_PREDGPI)
#   TargetP  targetp -format short -org non-pl (RUN_TARGETP)
#   Pfam     hmmsearch --cut_ga              (RUN_PFAM; the module runs it from the HMMER SIF)
#   dbCAN    dbcanlight search + conclude    (RUN_CAZY - this is what carries the substrate
#            names that tier T2 is defined on; an hmmsearch against dbCAN.hmm has none)
# predgpi, targetp and dbcanlight are conda environments behind their modules, so the same
# `conda shell.bash hook` the BFD profile uses is needed before loading them.
# Reads shared folders read-only; writes mucor_secretome/out/panels/ and data_external/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 24
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J annot_panels
#SBATCH -o logs/annot_panels_%j.out
#SBATCH -e logs/annot_panels_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
UPO=analysis/enzymes_UPO
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-24}
MOR=/bigdata/stajichlab/shared/projects/MorenoLab/Mucor_genomes/annotation
OUT=$A/out/panels; mkdir -p $OUT
FAA=$OUT/faa
W="${SCRATCH:?}/annot_panels"; mkdir -p $W

source /etc/profile.d/modules.sh 2>/dev/null || true
export CONDA_NO_PLUGINS=true
module load miniconda3
eval "$(conda shell.bash hook)"
module load tmhmm predgpi targetp dbcanlight hmmer/3.4 db-pfam
$PY -W ignore -m unittest discover -s $A/test

# PANEL_LOCUSTAGS lets a later genome be added without re-running the ones already in the
# tables; annotate_to_tables.py appends, so a second run of a genome would duplicate its rows.
LOCUSTAGS="${PANEL_LOCUSTAGS:-JEM1MUCOR JEM2MUCOR MPTAK MPALE}"
bash $A/bin/panel_proteomes.sh $FAA $LOCUSTAGS

for LT in $LOCUSTAGS; do
  F=$FAA/$LT.faa
  echo "=== $LT $(grep -c '>' $F) proteins $(date +%T)"
  tmhmm --short --noplot < $F > $W/$LT.tmhmm.txt 2> $W/$LT.tmhmm.err
  predgpi.py -f $F -m gff3 -o $W/$LT.predgpi.gff3
  TMPD=$(mktemp -d)
  targetp -batch 50 -tmp $TMPD -format short -fasta $F -org non-pl -prefix $W/$LT
  rm -rf $TMPD
  hmmsearch --cut_ga --noali --cpu $CPU --domtblout $W/$LT.pfam.domtbl \
      $PFAM_DB/Pfam-A.hmm $F > /dev/null
  CZ=""
  case "$LT" in MPTAK|MPALE) PLANT=1 ;; *) PLANT=0 ;; esac
  if [ "$PLANT" = "0" ]; then   # dbCAN is a fungal/microbial CAZyme set
    rm -rf $W/cz_$LT; mkdir -p $W/cz_$LT
    dbcanlight search -i $F -m cazyme -o $W/cz_$LT -t $CPU
    dbcanlight search -i $F -m sub    -o $W/cz_$LT -t $CPU
    dbcanlight conclude $W/cz_$LT
    CZ="--cazy-overview $W/cz_$LT/overview.tsv"
  fi
  $PY $A/bin/annotate_to_tables.py --locustag $LT --out $OUT \
      --tmhmm $W/$LT.tmhmm.txt --predgpi $W/$LT.predgpi.gff3 \
      --targetp ${W}/${LT}_summary.targetp2 \
      --pfam $W/$LT.pfam.domtbl $CZ --lengths $FAA/$LT.lengths.tsv
done

# --- real 5' ends for the two Mucor strains, from their own assemblies
mkdir -p $W/asm
case " $LOCUSTAGS " in *JEM1MUCOR*|*JEM2MUCOR*) ;; *) echo "no JEM panel, skipping codon check"; exit 0 ;; esac
for S in JEM1 JEM2; do
  gzip -c $MOR/$S/predict_results/Mucor_sp._$S.scaffolds.fa > $W/asm/$S.fa.gz
  $PY $A/bin/gff3_to_cds.py --gff3 $MOR/$S/predict_results/Mucor_sp._$S.gff3 \
      --asmid $S --set ${S}MUCOR --out-cds $W/$S.cds.tsv --out-models $W/$S.models.tsv
done
head -1 $W/JEM1.models.tsv > $W/jem.models.tsv
tail -q -n +2 $W/JEM1.models.tsv $W/JEM2.models.tsv >> $W/jem.models.tsv
head -1 $W/JEM1.cds.tsv > $W/jem.cds.tsv
tail -q -n +2 $W/JEM1.cds.tsv $W/JEM2.cds.tsv >> $W/jem.cds.tsv
cat $FAA/JEM1MUCOR.faa $FAA/JEM2MUCOR.faa > $W/jem.faa
$PY $UPO/bfd_findings/bin/check_codons.py --cds $W/jem.cds.tsv --models $W/jem.models.tsv \
    --genomes-dir $W/asm --proteins $W/jem.faa --threads $CPU --out $OUT/jem_codon_check.tsv
wc -l $OUT/*.tsv
echo "done $(date +%T)"
