#!/usr/bin/env python3.12
"""
G3 Phase A step 1b: pull promoter sequences for the intervals from make_promoter_coords.py.

Uses the bgzip/faidx index on each genome's scaffolds file, one samtools call per genome with
every interval for that genome. Minus-strand promoters are reverse-complemented so every sequence
reads 5'->3' toward the gene.

Writes one FASTA per (set) and reports every genome it could not open and every interval that
came back short or empty, rather than dropping them silently.
"""
import argparse, os, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--coords", required=True)
ap.add_argument("--genomes-tsv", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--samtools", default="samtools")
ap.add_argument("--min-len", type=int, default=300)
a = ap.parse_args()

paths = {}
with open(a.genomes_tsv) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            lt, folder, prot = f[0], f[1], f[2]
            paths[lt] = os.path.join(os.path.dirname(prot), f"{folder}.scaffolds.fa.gz")

jobs = defaultdict(list)
with open(a.coords) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 8:
            jobs[f[1]].append(f)   # by LOCUSTAG

os.makedirs(a.outdir, exist_ok=True)
outs = {s: open(os.path.join(a.outdir, f"{s}.fa"), "w") for s in ("target", "background")}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

n_ok = 0
bad_genomes, short = [], 0
for i, (lt, recs) in enumerate(sorted(jobs.items()), 1):
    p = paths.get(lt)
    if not p or not os.path.exists(p) or not os.path.exists(p + ".fai"):
        bad_genomes.append(lt); continue
    regions = [f"{r[2]}:{r[3]}-{r[4]}" for r in recs]
    try:
        res = subprocess.run([a.samtools, "faidx", p] + regions,
                             capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        bad_genomes.append(f"{lt}:timeout"); continue
    if res.returncode != 0:
        bad_genomes.append(f"{lt}:rc={res.returncode}"); continue
    seqs, cur, buf = [], None, []
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            if cur is not None: seqs.append("".join(buf))
            cur, buf = line, []
        elif line:
            buf.append(line.strip())
    if cur is not None: seqs.append("".join(buf))
    if len(seqs) != len(recs):
        bad_genomes.append(f"{lt}:got {len(seqs)} of {len(recs)}"); continue
    for r, s in zip(recs, seqs):
        if len(s) < a.min_len:
            short += 1; continue
        if int(r[5]) < 0:
            s = rc(s)
        outs[r[6]].write(f">{r[0]}|{r[1]}|{r[6]}\n{s}\n")
        n_ok += 1
    if i % 2000 == 0:
        print(f"  {i}/{len(jobs)} genomes, {n_ok} promoters", file=sys.stderr)

for o in outs.values(): o.close()
print(f"\npromoters written: {n_ok}", file=sys.stderr)
for s in ("target", "background"):
    p = os.path.join(a.outdir, f"{s}.fa")
    n = sum(1 for l in open(p) if l.startswith(">"))
    print(f"  {s}: {n}", file=sys.stderr)
print(f"genomes unreadable: {len(bad_genomes)}", file=sys.stderr)
print(f"intervals returning < {a.min_len} bp: {short}", file=sys.stderr)
if bad_genomes:
    open(os.path.join(a.outdir, "MISSING_genomes.txt"), "w").write("\n".join(bad_genomes) + "\n")
