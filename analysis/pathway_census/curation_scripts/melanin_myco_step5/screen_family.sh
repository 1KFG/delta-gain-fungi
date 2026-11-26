#!/bin/bash
# Alignability screens used before make_leaves.sh (record of the commands that were run; needs module hmmer/3.4 on PATH).
# 1. family FASTA (fungal + other reviewed entries, fetch_family.py) against the step HMM
# 2. decoys of refs/<pathway>/decoys against the step HMM (best and union coverage), filter_decoys_for_tree.py, min 10 bits and 0.5
# 3. references of sibling steps against the step HMM (laccase/MCO and tyrosinase groups)
I=results/step5/other_inputs; mkdir -p $I/screen $I/decoys $I/sibling
declare -A FAM=( [ayg1]="melanin_dhn PF06500" [scytalone_dehydratase]="melanin_dhn PF02982" [t4hn_reductase]="melanin_dhn PF13561" [t3hn_reductase]="melanin_dhn PF13561" [laccase_abr2]="melanin_dhn PF00394" [mco_abr1]="melanin_dhn PF00394" [laccase_diphenol]="melanin_ldopa PF00394" [tyrosinase]="melanin_ldopa PF00264" [asp_melanin_tyrosinase]="melanin_ldopa PF00264" [hmga]="melanin_ldopa PF04209" [hppd]="melanin_ldopa PF00903" [aspulvinone_e_synthase]="melanin_ldopa PF00501" [ddgs]="mycosporine PF01761" [omt]="mycosporine PF01596" [atp_grasp_ligase]="mycosporine PF02655" )
for s in "${!FAM[@]}"; do
  set -- ${FAM[$s]}; p=$1; pf=$2
  cat $I/family/$pf.fungal.faa $I/family/$pf.other.faa > $I/screen/$s.family.faa
  /usr/bin/python3.12 bin/screen_outgroups.py --fasta $I/screen/$s.family.faa --hmm results/hmms/$p/$s.PROVISIONAL.hmm --out $I/screen/$s.family_screen.tsv
  for m in best union; do
    /usr/bin/python3.12 bin/filter_decoys_for_tree.py --decoys refs/$p/decoys/$s.faa --hmm results/hmms/$p/$s.PROVISIONAL.hmm --min-score 10 \
      --coverage-mode $m --out-fasta $I/decoys/$s.$m.faa --out-table $I/decoys/$s.$m.tsv
  done
done
# atp_grasp_ligase: five more ATP-grasp Pfam families (PF02222 PF02786 PF07478 PF08443 PF13535; fetch_family.py --max-other 300), 1,100 proteins, 0 pass
cat $I/family/PF0{2222,2786,7478}.*.faa $I/family/PF{08443,13535}.*.faa > $I/screen/atp_grasp_ligase.atpfam.faa
/usr/bin/python3.12 bin/screen_outgroups.py --fasta $I/screen/atp_grasp_ligase.atpfam.faa --hmm results/hmms/mycosporine/atp_grasp_ligase.PROVISIONAL.hmm --out $I/screen/atp_grasp_ligase.atpfam_screen.tsv
# sibling references (refs of the other steps of the group, file $I/sibling/<step>.refs.faa = cat of those refs files)
for s in laccase_abr2 mco_abr1 laccase_diphenol tyrosinase asp_melanin_tyrosinase; do
  p=${FAM[$s]%% *}
  /usr/bin/python3.12 bin/screen_outgroups.py --fasta $I/sibling/$s.refs.faa --hmm results/hmms/$p/$s.PROVISIONAL.hmm --out $I/sibling/$s.screen.tsv
done
