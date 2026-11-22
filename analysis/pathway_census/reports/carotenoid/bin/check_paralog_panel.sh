#!/usr/bin/env bash
# Is each candidate representative closer to the step references than to a panel of reviewed fungal paralogs?
# Used for phytoene_synthase against squalene synthase (UniProt EC 2.5.1.21). Needs internet (UniProt REST) and
# ncbi-blast. No BASH_SOURCE. Usage: check_paralog_panel.sh <step5 dir> <refs fasta> <out dir> <paralog EC> <panel name>
#   example: check_paralog_panel.sh results/step5/carotenoid/phytoene_synthase refs/carotenoid/phytoene_synthase.faa \
#            results/step5/carotenoid/phytoene_synthase/paralog_panel 2.5.1.21 ERG9
set -euo pipefail
S5=${1:?step5 dir}; REFS=${2:?refs fasta}; OUT=${3:?out dir}; EC=${4:?EC number}; PANEL=${5:?panel name}
module load ncbi-blast/2.17.0+ >/dev/null 2>&1 || true
mkdir -p "$OUT"
Q="(ec:${EC}) AND (taxonomy_id:4751) AND (reviewed:true)"
curl -sf -m 120 -G "https://rest.uniprot.org/uniprotkb/stream" --data-urlencode "query=${Q}" --data-urlencode "format=fasta" -o "$OUT/panel_raw.faa"
echo "query: ${Q}; date: $(date +%F); entries: $(grep -c '>' "$OUT/panel_raw.faa")" > "$OUT/panel_query.txt"
sed "s/^>sp|\([A-Z0-9]*\)|.*/>${PANEL}_\1/" "$OUT/panel_raw.faa" > "$OUT/panel.faa"
awk '/^>/{split($1,a," ");print ">REF_" substr(a[1],2);next}{print}' "$REFS" > "$OUT/refs.faa"
cat "$OUT/panel.faa" "$OUT/refs.faa" > "$OUT/db.faa"
makeblastdb -in "$OUT/db.faa" -dbtype prot -out "$OUT/db" >/dev/null
blastp -query "$S5/reps.faa" -db "$OUT/db" -evalue 1e-3 -max_target_seqs 50 -outfmt "6 qseqid sseqid pident length bitscore" \
  -num_threads "${SLURM_CPUS_ON_NODE:-4}" > "$OUT/hits.tsv"
/usr/bin/python3.12 -I - "$OUT" "$PANEL" <<'PY'
import collections, statistics, sys
out, panel = sys.argv[1], sys.argv[2]
best = collections.defaultdict(lambda: {panel: 0.0, "REF": 0.0})
for line in open(f"{out}/hits.tsv"):
    q, s, pid, ln, bit = line.rstrip("\n").split("\t")
    k = s.split("_")[0]
    best[q][k] = max(best[q][k], float(bit))
n = len(best)
closer_ref = sum(v["REF"] > v[panel] for v in best.values())
closer_panel = sum(v[panel] > v["REF"] for v in best.values())
d = [v["REF"] - v[panel] for v in best.values()]
with open(f"{out}/summary.tsv", "w") as fh:
    fh.write("key\tvalue\n")
    for k, v in (("representatives_with_a_hit", n), ("closer_to_references", closer_ref), (f"closer_to_{panel}_panel", closer_panel),
                 ("margin_bits_min", min(d)), ("margin_bits_median", statistics.median(d)), ("margin_bits_max", max(d)),
                 (f"median_best_{panel}_bits", statistics.median(v[panel] for v in best.values())),
                 ("median_best_reference_bits", statistics.median(v["REF"] for v in best.values()))):
        fh.write(f"{k}\t{v}\n")
print(open(f"{out}/summary.tsv").read())
PY
