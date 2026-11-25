#!/bin/bash
# Select family leaves per step (bin/select_family_leaves.py). Run from analysis/pathway_census. step|pathway|pfam|exclude ids|exclude EC
I=results/step5/other_inputs
ALLREFS=Q4WZB3,Q6BCB5,A2R3S8,E9RBR0,Q4WZB4,A2R3T0,Q00455,O14434,P56221,A7DWM5,Q12634,P87025,A7DWM6,Q9HFV6,E9QUT3,W3XAA0,F8WTI4,A7DWM4,Q0CRX0,Q0CRX1,Q4WHT9,Q00667,Q4WHU1,Q55P57,J9VY90,J9VQZ4,P00440,A0A1S9DK56,Q00234,C7FF05,C7FF04,A7BHQ9,Q2TL94,A0A9Y1GGR7,B2J6X7,B2J6X9,Q3M6C3,B2J6X8
while IFS='|' read -r s p pf ec; do
  [ -z "$s" ] && continue
  /usr/bin/python3.12 bin/select_family_leaves.py --screen $I/screen/$s.family_screen.tsv \
    --meta $I/family/$pf.fungal.tsv $I/family/$pf.other.tsv --fasta $I/family/$pf.fungal.faa $I/family/$pf.other.faa \
    --hmm $s --mode best --max-score $(awk -F"\t" -v s=$s '$1==s{print $2}' results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv) --exclude-gene "$(awk -F"\t" -v s=$s '$1==s{print $2}' $I/ref_gene_tokens.tsv)" --exclude-ids $ALLREFS --exclude-ec "$ec" --n-named 10 --n-outgroup 6 \
    --out-named $I/leaves/$s.named.faa --out-outgroups $I/leaves/$s.outgroups.faa \
    --out-provenance $I/leaves/$s.leaves_provenance.tsv --out-exclusions $I/leaves/$s.exclusions.tsv
done <<'LIST'
ayg1|melanin_dhn|PF06500|
scytalone_dehydratase|melanin_dhn|PF02982|4.2.1.94
t4hn_reductase|melanin_dhn|PF13561|1.1.1.252,1.1.1.251
t3hn_reductase|melanin_dhn|PF13561|1.1.1.252,1.1.1.251
laccase_abr2|melanin_dhn|PF00394|
mco_abr1|melanin_dhn|PF00394|
laccase_diphenol|melanin_ldopa|PF00394|
tyrosinase|melanin_ldopa|PF00264|
asp_melanin_tyrosinase|melanin_ldopa|PF00264|
hmga|melanin_ldopa|PF04209|1.13.11.5
hppd|melanin_ldopa|PF00903|1.13.11.27
aspulvinone_e_synthase|melanin_ldopa|PF00501|
ddgs|mycosporine|PF01761|4.2.3.154
omt|mycosporine|PF01596|
LIST
