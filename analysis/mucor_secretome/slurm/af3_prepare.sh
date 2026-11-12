#!/bin/bash
# Assemble the two AF3 panels, then hand them to af3_msa.sh / af3_fold.sh.
#   sbatch analysis/mucor_secretome/slurm/af3_prepare.sh
# Panel 1 "chitin": the proteins whose ESMFold structure was confidently assigned to a LysM,
# deacetylase, cerato-platanin or chitinase fold, against the LysM receptor ectodomains and
# chitinases of both Marchantia species, plus chito-oligosaccharides at CO4 and CO7.
# Panel 2 "unknowns": the candidates ESMFold could not fold (pLDDT < 70). Most are
# mucorales_wide, so homologues exist and an MSA-based predictor has something to work with -
# which is exactly what ESMFold, being single-sequence, did not.
# Control sequences come from RCSB (public structures, no project data sent).
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -J af3_prep
#SBATCH -o logs/af3_prep_%j.out
#SBATCH -e logs/af3_prep_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
E=$A/out/effectors
CH=$A/out/af3/chitin;   mkdir -p $CH
UN=$A/out/af3/unknowns; mkdir -p $UN
$PY -W ignore -m unittest discover -s $A/test

# --- panel 1: the chitin axis
$PY $A/bin/build_chitin_panel.py --fold-calls $E/fold_calls.tsv --folds-faa $E/to_fold.faa \
    --plant-panel $E/plant_panel/panel.tsv --plant-faa $E/plant_panel/panel.faa \
    --out $CH

# Controls, fetched from RCSB: the construct in the solved structure is the one worth
# reproducing. 4EBZ = AtCERK1 LysM ectodomain with chito-oligosaccharide (protein/ligand);
# 6Q76 = AVR-PikD effector with the Pikp-1 HMA domain (small fungal effector/plant target).
for pdb in 4EBZ 6Q76; do
  curl -sS --max-time 300 "https://www.rcsb.org/fasta/entry/$pdb/display" -o $CH/$pdb.fa.raw
  $PY - "$CH/$pdb.fa.raw" "$pdb" "$CH/$pdb.fa" <<'PYEOF'
import sys, re
src, pdb, dst = sys.argv[1:4]
# one entry per distinct polymer chain; name them CTRL_<pdb>_<n> and drop duplicates
seqs, name = [], None
for line in open(src):
    if line.startswith(">"):
        name = f"CTRL_{pdb}_{len(seqs) + 1}"
        seqs.append([name, ""])
    elif seqs:
        seqs[-1][1] += line.strip()
seen = set()
with open(dst, "w") as fo:
    for n, s in seqs:
        s = re.sub(r"[^A-Z]", "", s.upper())
        if len(s) < 40 or s in seen:
            continue
        seen.add(s)
        fo.write(f">{n}\n{s}\n")
print(f"{pdb}: {len(seen)} control chains", file=sys.stderr)
PYEOF
  cat $CH/$pdb.fa >> $CH/panel.faa
done

# ligand jobs: every fungal chain plus the 4EBZ control receptor, at CO4 and CO7
cut -f1 $CH/fungal.txt > $CH/ligand_ids.txt
grep '^>' $CH/4EBZ.fa | tr -d '>' >> $CH/ligand_ids.txt
# protein-protein controls: the 6Q76 pair is the positive, and crossing a control chain with an
# unrelated host chain is the negative that says what "no interaction" scores in this same run
printf 'job\trole\n' > $CH/controls.tsv
A6=$(grep '^>' $CH/6Q76.fa | tr -d '>' | sed -n 1p)
B6=$(grep '^>' $CH/6Q76.fa | tr -d '>' | sed -n 2p)
[ -n "${A6:-}" ] && [ -n "${B6:-}" ] && printf '%s__%s\tpositive\n' "$A6" "$B6" >> $CH/controls.tsv
C4=$(grep '^>' $CH/4EBZ.fa | tr -d '>' | sed -n 1p)
[ -n "${C4:-}" ] && printf '%s__CO7\tpositive_ligand\n%s__CO4\tpositive_ligand\n' "$C4" "$C4" >> $CH/controls.tsv
HOST1=$(sed -n 1p $CH/host.txt)
[ -n "${A6:-}" ] && [ -n "${HOST1:-}" ] && printf '%s__%s\tnegative\n' "$A6" "$HOST1" >> $CH/controls.tsv
# the control chains must also be in the pair lists for those jobs to be generated
{ [ -n "${A6:-}" ] && echo "$A6"; [ -n "${C4:-}" ] && echo "$C4"; } >> $CH/fungal.txt
[ -n "${B6:-}" ] && echo "$B6" >> $CH/host.txt
sort -u -o $CH/fungal.txt $CH/fungal.txt
sort -u -o $CH/host.txt $CH/host.txt
wc -l $CH/fungal.txt $CH/host.txt $CH/ligand_ids.txt $CH/controls.tsv

# --- panel 2: the unfoldable unknowns, where an MSA is the thing that was missing
$PY - "$E/fold_calls.tsv" "$E/to_fold.faa" "$UN/panel.faa" <<'PYEOF'
import csv, sys
calls, faa, out = sys.argv[1:4]
want = set()
with open(calls) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("fold_call") in ("low_confidence", "no_match"):
            want.add(r["protein_id"])
n, name, buf = 0, None, []
with open(faa) as fh, open(out, "w") as fo:
    def flush():
        global n
        if name and name in want and 40 <= len("".join(buf)) <= 400:
            fo.write(f">{name}\n{''.join(buf)}\n"); n += 1
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush(); name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    flush()
print(f"{n} unfoldable candidates -> {out}", file=sys.stderr)
PYEOF
grep -c '>' $UN/panel.faa
echo "done $(date +%T)"
