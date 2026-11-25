#!/bin/bash
# Assemble tree inputs per step: decoys (alignable, cut), named members (family + alignable sibling-step references, cut),
# outgroups (cut). Run from analysis/pathway_census after make_leaves.sh. Needs module hmmer/3.4.
I=results/step5/other_inputs
ids() { grep '>' "$1" | sed 's/^>//' | awk '{print $1}'; }
drop_ids() { # drop_ids <in.faa> <ids file> <out.faa>
  /usr/bin/python3.12 - "$1" "$2" "$3" <<'PY'
import sys
drop = {l.strip() for l in open(sys.argv[2]) if l.strip()}
keep = False
with open(sys.argv[1]) as fi, open(sys.argv[3], "w") as fo:
    for line in fi:
        if line.startswith(">"):
            keep = line[1:].split()[0] not in drop
        if keep:
            fo.write(line)
PY
}
declare -A PW=( [ayg1]=melanin_dhn [scytalone_dehydratase]=melanin_dhn [t4hn_reductase]=melanin_dhn [t3hn_reductase]=melanin_dhn [laccase_abr2]=melanin_dhn [mco_abr1]=melanin_dhn [laccase_diphenol]=melanin_ldopa [tyrosinase]=melanin_ldopa [asp_melanin_tyrosinase]=melanin_ldopa [hmga]=melanin_ldopa [hppd]=melanin_ldopa [aspulvinone_e_synthase]=melanin_ldopa [ddgs]=mycosporine [omt]=mycosporine [atp_grasp_ligase]=mycosporine )
for s in ${ONLY:-ayg1 scytalone_dehydratase laccase_abr2 mco_abr1 laccase_diphenol tyrosinase asp_melanin_tyrosinase hppd aspulvinone_e_synthase ddgs omt}; do
  p=${PW[$s]}; H=results/hmms/$p/$s.PROVISIONAL.hmm; T=$I/tree_inputs/$s; mkdir -p $T
  # decoys that pass the alignability rule (best-domain coverage), cut to the domain
  /usr/bin/python3.12 bin/cut_leaves_to_domain.py --fasta $I/decoys/$s.best.faa --hmm $H --out-fasta $T/decoys.faa --out-table $T/decoys_cut.tsv | tail -1
  ids $T/decoys.faa > $T/decoy_ids.txt
  # sibling-step references that pass the screen in best mode
  : > $T/sibling_pass.faa
  if [ -s $I/sibling/$s.screen.tsv ]; then
    /usr/bin/python3.12 bin/select_by_screen.py --fasta $I/sibling/$s.refs.faa --screen $I/sibling/$s.screen.tsv --hmm $s --mode best --out $T/sibling_pass.faa | tail -1
  fi
  cat $I/leaves/$s.named.faa $T/sibling_pass.faa > $T/members_raw.faa
  drop_ids $T/members_raw.faa $T/decoy_ids.txt $T/members_nodecoy.faa
  drop_ids $I/leaves/$s.outgroups.faa $T/decoy_ids.txt $T/outgroups_nodecoy.faa
  /usr/bin/python3.12 bin/cut_leaves_to_domain.py --fasta $T/members_nodecoy.faa --hmm $H --out-fasta $T/members.faa --out-table $T/members_cut.tsv | tail -1
  /usr/bin/python3.12 bin/cut_leaves_to_domain.py --fasta $T/outgroups_nodecoy.faa --hmm $H --out-fasta $T/outgroups.faa --out-table $T/outgroups_cut.tsv | tail -1
  echo "$s: decoys $(grep -c '>' $T/decoys.faa) members $(grep -c '>' $T/members.faa) outgroups $(grep -c '>' $T/outgroups.faa)"
done
