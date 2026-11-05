#!/usr/bin/env bash
# Fill the prior-annotation gap in NCBI_fungi/source/NCBI_ASM.
#  (a) 44 annotated assemblies in lib/ncbi_accessions.csv but never downloaded -> project sync script
#  (b) 158 folders that have a genome but no protein file -> fetch protein+gff3+cds only, add to the existing folder
set -uo pipefail
N=/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi
D=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi
export PATH=$N/bin:$PATH
OUT=$N/source/NCBI_ASM
LOG=$D/logs/fetch_gap.log; mkdir -p $D/logs
echo "start $(date)" | tee -a $LOG
echo "== (a) new assemblies"
tail -n +2 $D/gap_notdownloaded.csv | while IFS=, read -r ACC SP ST TX BP LEN N50 ASMNAME FOLDER; do
  $N/scripts/sync_ncbi_assembly.sh --method datasets "$ACC" "$ASMNAME" "$FOLDER" "$OUT" >>$LOG 2>&1 && echo "ok $ACC" || echo "FAIL $ACC" | tee -a $LOG
  sleep 1
done
echo "== (b) add proteins to existing folders"
while IFS=$'\t' read -r ACC FOLDER; do
  T=$(mktemp -d); T2=$OUT/$FOLDER
  if datasets download genome accession "$ACC" --filename $T/d.zip --include gff3,protein,cds >>$LOG 2>&1 && unzip -q $T/d.zip -d $T; then
     S=$T/ncbi_dataset/data/$ACC
     for pair in protein.faa:protein.faa.gz genomic.gff:genomic.gff.gz cds_from_genomic.fna:cds_from_genomic.fna.gz; do
        src=${pair%%:*}; dst=${pair##*:}
        [[ -s $S/$src && ! -e $T2/${FOLDER}_$dst ]] && gzip -c $S/$src > $T2/${FOLDER}_$dst
     done
     [[ -s $T2/${FOLDER}_protein.faa.gz ]] && echo "ok $ACC" || echo "NOPROT $ACC" | tee -a $LOG
  else echo "FAIL $ACC" | tee -a $LOG; fi
  rm -rf $T; sleep 1
done < $D/gap_noprotein_folders.tsv
echo "done $(date)" | tee -a $LOG
