#!/usr/bin/env python3.12
"""Subset promoter FASTA by lineage and clade, using the representative mapping."""
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--map", required=True)
ap.add_argument("--lineage", required=True)
ap.add_argument("--clade", default=None)
ap.add_argument("--any-clade", action="store_true")
ap.add_argument("--out", required=True)
a = ap.parse_args()

keep_lt = set()
keep_gene = set()
with open(a.map) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    i = {k: j for j, k in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[i["LINEAGE_GROUP"]] != a.lineage:
            continue
        keep_lt.add(f[i["gene_id"]][:8])
        if a.any_clade or f[i["clade"]] == a.clade:
            keep_gene.add(f[i["gene_id"]])

out = open(a.out, "w"); w = False
for line in open(a.fasta):
    if line.startswith(">"):
        parts = line[1:].strip().split("|")
        gid, lt = parts[0], (parts[1] if len(parts) > 1 else parts[0][:8])
        # background sequences are selected by genome (lineage), targets by gene (clade)
        w = (lt in keep_lt) if a.any_clade else (gid in keep_gene)
    if w:
        out.write(line)
out.close()
