#!/usr/bin/env python3
"""How much does keeping more than 3 MIBiG hits per region add? Input: dir with <id>.known.tsv.gz and
<id>.known_pairs.tsv.gz written with --topn 10. Output: counts and table sizes at TOPN 3, 5, 10 and mixed (known 5, pairs 3)."""
import csv, glob, gzip, io, os, sys, collections
d = sys.argv[1]
ids = sorted(os.path.basename(f).split(".")[0] for f in glob.glob(os.path.join(d, "*.known.tsv.gz")))
def sim(x):
    try: return int(float(x))
    except: return -1
tot = collections.Counter(); ex = []
size = collections.defaultdict(int)
def gz_size(header, rows):
    b = io.BytesIO()
    with gzip.GzipFile(fileobj=b, mode="wb", mtime=0) as g:
        g.write(("\t".join(header) + "\n").encode())
        for r in rows: g.write(("\t".join(r) + "\n").encode())
    return len(b.getvalue())
for g in ids:
    with gzip.open(os.path.join(d, g + ".known.tsv.gz"), "rt") as fh:
        rd = csv.reader(fh, delimiter="\t"); kh = next(rd); krows = list(rd)
    with gzip.open(os.path.join(d, g + ".known_pairs.tsv.gz"), "rt") as fh:
        rd = csv.reader(fh, delimiter="\t"); ph = next(rd); prows = list(rd)
    ki = {c: i for i, c in enumerate(kh)}; pi = {c: i for i, c in enumerate(ph)}
    byreg = collections.defaultdict(list)
    for r in krows: byreg[(r[ki["record"]], r[ki["region"]])].append(r)
    for key, hits in byreg.items():
        tot["regions_with_hit"] += 1
        def best(n):
            return max((sim(h[ki["similarity"]]) for h in hits if int(h[ki["rank"]]) <= n), default=-1)
        b3, b5, b10 = best(3), best(5), best(10)
        tot["n_hits_total"] += len(hits)
        if b10 > b3: tot["best_outside_1-3"] += 1
        if b10 > b5: tot["best_outside_1-5"] += 1
        if b5 > b3: tot["best_in_4-5_higher_than_1-3"] += 1
        if b5 >= 75 and b5 > b3:
            tot["rank4-5_ge75_and_higher"] += 1
            h = max((h for h in hits if int(h[ki["rank"]]) <= 5 and sim(h[ki["similarity"]]) == b5), key=lambda h: 0)
            ne = int(h[ki["n_entry_proteins"]]) if h[ki["n_entry_proteins"]] not in ("", "None") else 0
            tot["  of which n_entry_proteins==1"] += ne == 1
            tot["  of which n_entry_proteins>=3"] += ne >= 3
            ex.append((g, key, b3, b5, ne, h[ki["mibig_accession"]], h[ki["hit_description"]][:40]))
        if b10 >= 75 and b3 < 75: tot["status_changes_weak_to_known_if_top10"] += 1
        if b5 >= 75 and b3 < 75: tot["status_changes_weak_to_known_if_top5"] += 1
    for n in (3, 5, 10):
        size[f"known_top{n}"] += gz_size(kh, [r for r in krows if int(r[ki["rank"]]) <= n])
        size[f"pairs_top{n}"] += gz_size(ph, [r for r in prows if int(r[pi["rank"]]) <= n])
        tot[f"pair_rows_top{n}"] += sum(1 for r in prows if int(r[pi["rank"]]) <= n)
        tot[f"known_rows_top{n}"] += sum(1 for r in krows if int(r[ki["rank"]]) <= n)
n = len(ids)
print(n, "genomes:", ", ".join(ids))
for k, v in sorted(tot.items()): print(f"{k}: {v}")
rg = tot["regions_with_hit"]
for k in ("best_outside_1-3", "best_outside_1-5", "best_in_4-5_higher_than_1-3", "rank4-5_ge75_and_higher", "status_changes_weak_to_known_if_top5", "status_changes_weak_to_known_if_top10"):
    print(f"fraction {k}: {tot[k]}/{rg} = {100*tot[k]/rg:.2f}%")
print("gz TSV bytes per genome (mean):")
for k, v in sorted(size.items()): print(f"  {k}: {v/n:.0f}")
print("mixed (known top5 + pairs top3):", round((size["known_top5"] + size["pairs_top3"]) / n), "per genome; top3 both:", round((size["known_top3"] + size["pairs_top3"]) / n), "; top5 both:", round((size["known_top5"] + size["pairs_top5"]) / n))
print("examples (genome, region, best1-3, best1-5, n_entry, mibig, desc):")
for e in ex[:12]: print(" ", e)
