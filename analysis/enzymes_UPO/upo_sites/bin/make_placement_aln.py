#!/usr/bin/env python3.12
"""
Build a reference + query alignment on PF01328 MATCH STATES only, for phylogenetic placement.

Match-state columns are fixed by the profile, so reference and query alignments are directly
comparable and chunking is safe. Columns that are gaps in most references are dropped, because
RAxML's EPA behaves badly on near-empty columns.
"""
import argparse, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from upo_sites import run_hmmalign, parse_stockholm, match_column_indices, read_fasta

ap = argparse.ArgumentParser()
ap.add_argument("--hmm", required=True)
ap.add_argument("--refs", required=True)
ap.add_argument("--queries", default=None)
ap.add_argument("--out-ref", required=True)
ap.add_argument("--out-all", default=None)
ap.add_argument("--max-ref-gap", type=float, default=0.5,
                help="drop match columns that are gaps in more than this fraction of REFERENCES")
a = ap.parse_args()

work = tempfile.mkdtemp(prefix="placealn.")
cat = os.path.join(work, "all.fa")
nref = 0
with open(cat, "w") as out:
    for n, s in read_fasta(a.refs):
        out.write(f">R_{n.split()[0]}\n{s}\n"); nref += 1
    nq = 0
    if a.queries:
        for n, s in read_fasta(a.queries):
            out.write(f">Q_{n.split()[0]}\n{s}\n"); nq += 1
print(f"{nref} references, {nq} queries", file=sys.stderr)

sto = os.path.join(work, "all.sto")
run_hmmalign(a.hmm, cat, sto)
seqs, rf = parse_stockholm(sto)
mcols = match_column_indices(rf)
print(f"{len(mcols)} match-state columns", file=sys.stderr)

names = list(seqs)
mat = {n: "".join(seqs[n][i] for i in mcols).upper().replace(".", "-") for n in names}
refs = [n for n in names if n.startswith("R_")]
keep = []
for j in range(len(mcols)):
    g = sum(1 for n in refs if mat[n][j] == "-")
    if g / max(len(refs), 1) <= a.max_ref_gap:
        keep.append(j)
print(f"{len(keep)} columns kept after gap filter", file=sys.stderr)

def write(path, sel):
    with open(path, "w") as fh:
        for n in sel:
            fh.write(f">{n[2:]}\n{''.join(mat[n][j] for j in keep)}\n")

write(a.out_ref, refs)
if a.out_all:
    write(a.out_all, refs + [n for n in names if n.startswith("Q_")])
import shutil; shutil.rmtree(work, ignore_errors=True)
print("done", file=sys.stderr)
