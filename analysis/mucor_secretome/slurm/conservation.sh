#!/bin/bash
# How lineage-restricted is each candidate effector? CPU only - no GPU queue.
#   sbatch analysis/mucor_secretome/slurm/conservation.sh
# Clusters the JEM candidate secretomes together with the 195,284 candidates of all 292
# Mucoromycota genomes in BFD, at the same 50 % identity / 80 % coverage the Mucoromycota
# secretome used, and asks how far each JEM candidate's family reaches. A secreted protein found
# in one strain and nowhere else is the sequence-independent effector signal; paired with the
# fold calls it gives the two axes that matter, neither of which stands alone.
# Reads shared folders read-only; writes mucor_secretome/out/effectors/.
#SBATCH -p epyc,intel
#SBATCH --constraint="ryzen|broadwell|cascade|sapphire"
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 8:00:00
#SBATCH -J eff_cons
#SBATCH -o logs/eff_cons_%j.out
#SBATCH -e logs/eff_cons_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
A=analysis/mucor_secretome
PY=/usr/bin/python3.12
CPU=${SLURM_CPUS_ON_NODE:-32}
CACHE=${CACHE:-$A/cache}; CACHE=$(readlink -f "$CACHE")
PANELS=$A/out/panels
OUT=$A/out/effectors; mkdir -p $OUT
W="${SCRATCH:?}/eff_cons"; mkdir -p $W
module load MMseqs2/17-b804f
$PY -W ignore -m unittest discover -s $A/test

[ -s $CACHE/candidates.faa ] || { echo "no $CACHE/candidates.faa - run slurm/secretome.sh" >&2; exit 1; }

# The JEM candidates go in whole, not just the fold set: a family's reach is only right if every
# member of it is present to be clustered.
$PY $A/bin/call_effectors.py --dir $PANELS --locustags JEM1MUCOR,JEM2MUCOR \
    --fasta $PANELS/faa/JEM1MUCOR.faa,$PANELS/faa/JEM2MUCOR.faa \
    --codons $PANELS/jem_codon_check.tsv --max-len 100000 --out $W/all
cat $CACHE/candidates.faa $W/all/to_fold.faa > $W/pool.faa
grep -c '>' $W/pool.faa

# The two JEM strains are not in BFD, so they need a genome row of their own for the species,
# genus and order a cluster's reach is counted over. ITS put JEM1 in Mucor hiemalis and JEM2 in
# Mucor circinelloides (out/its/), which is what these rows record.
printf 'LOCUSTAG\tASMID\tSPECIES\tSTRAIN\tORD\tFAMILY\tGENUS\n' > $W/jem_genomes.tsv
printf 'JEM1MUCOR\tJEM1\tMucor hiemalis\tJEM1\tMucorales\tMucoraceae\tMucor\n' >> $W/jem_genomes.tsv
printf 'JEM2MUCOR\tJEM2\tMucor circinelloides\tJEM2\tMucorales\tMucoraceae\tMucor\n' >> $W/jem_genomes.tsv

if [ ! -s $OUT/pool_cluster.tsv ]; then
  mmseqs easy-cluster $W/pool.faa $W/clu $W/mmtmp \
      --min-seq-id 0.5 -c 0.8 --cov-mode 5 --threads $CPU -v 1
  cp $W/clu_cluster.tsv $OUT/pool_cluster.tsv
else
  echo "reusing $OUT/pool_cluster.tsv"
fi
wc -l $OUT/pool_cluster.tsv

$PY $A/bin/conservation.py --clusters $OUT/pool_cluster.tsv --genomes $CACHE/genomes.tsv \
    --extra-genomes $W/jem_genomes.tsv --focus JEM1MUCOR,JEM2MUCOR \
    --out $OUT/jem_conservation.tsv

# Control. JEM2 (M. circinelloides) has 24 conspecific genomes in BFD; JEM1 (M. hiemalis) has
# none, because BFD's single "M. hiemalis" assembly is mislabelled (out/its/). So JEM1's
# genome_specific count is inflated relative to JEM2's by reference sampling alone. Removing
# JEM2's conspecifics from the genome table puts the two on the same footing: if JEM2's
# genome_specific count then rises to JEM1's, the asymmetry was sampling, not biology.
$PY $A/bin/conservation.py --clusters $OUT/pool_cluster.tsv --genomes $CACHE/genomes.tsv \
    --extra-genomes $W/jem_genomes.tsv --focus JEM1MUCOR,JEM2MUCOR \
    --drop-species "Mucor circinelloides" \
    --out $OUT/jem_conservation_no_conspecifics.tsv

# join onto the candidate table, so one file carries tier, fold reason and reach
$PY - "$OUT/jem/candidates.tsv" "$OUT/jem_conservation.tsv" "$OUT/jem_candidates_annotated.tsv" <<'PYEOF'
import csv, sys
cand, cons, out = sys.argv[1:4]
with open(cons) as fh:
    c = {r["protein_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
add = ["conservation", "cluster_members", "cluster_genomes", "cluster_species",
       "cluster_genera", "cluster_orders"]
with open(cand) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    rows = list(rd)
    cols = rd.fieldnames + add
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        r.update({k: c.get(r["protein_id"], {}).get(k, "") for k in add})
        w.writerow(r)
print(f"{len(rows)} rows -> {out}", file=sys.stderr)
PYEOF

echo "=== fold-set candidates by reach ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["fold_reason"]!="" {print $h["LOCUSTAG"]"\t"$h["conservation"]"\t"$h["fold_reason"]}' \
    $OUT/jem_candidates_annotated.tsv | sort | uniq -c | sort -rn
echo "done $(date +%T)"
