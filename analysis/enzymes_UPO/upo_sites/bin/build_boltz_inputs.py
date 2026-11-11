#!/usr/bin/env python3.12
"""
Build Boltz-2 inputs (a3m MSA + YAML with heme) for the 19 calibration enzymes.

WHY THIS EXISTS
---------------
Only 5 of the 19 calibration enzymes have crystal structures, so the channel-width versus
naphthalene test rests on n = 5 — and on n = 4 once TteUPO is set aside as flat across all five
assays. All 19 have naphthalene measurements in Figure S4. Modelling the other 14 takes the test
to n = 19.

The 5 that already have crystal structures are modelled too, deliberately. Comparing predicted
against experimental channel width on those 5 is the validation that decides whether a predicted
structure can be measured this way at all. If it fails there, the other 14 must not be used.

MSAs FROM REAL HOMOLOGS
-----------------------
Rather than depend on an MSA server the compute nodes may not reach, the MSA is built from this
project's own data: dereplicated BFD PF01328 representatives aligned to the Pfam profile with
hmmalign. Having 97,520 homologs on disk is what makes this possible.

A trap worth recording: the first version of this built the a3m from PF01328 MATCH STATES only,
which silently truncated HspUPO from 261 residues to 181 — dropping real protein that the
structure prediction needs. The a3m must be QUERY-LENGTH: keep every column where the query has a
residue, insert columns included. The assertion below enforces that and will fail loudly rather
than hand Boltz a truncated protein.
"""
from __future__ import annotations
import argparse, os, random, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from upo_sites import parse_stockholm

ap = argparse.ArgumentParser()
ap.add_argument("--refs", required=True, help="data/clade_backbone/raw.fa (the 19)")
ap.add_argument("--homologs", required=True, help="dereplicated PF01328 FASTA for the MSA pool")
ap.add_argument("--hmm", required=True)
ap.add_argument("--n-homologs", type=int, default=2000)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--outdir", required=True, help="data/boltz")
ap.add_argument("--ligand-ccd", default="HEM")
a = ap.parse_args()


def read_fa(p, strip_bar=False):
    d, n = {}, None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            n = l[1:].split()[0]
            if strip_bar:
                n = n.split("|")[0]
            d[n] = ""
        else:
            d[n] += l.strip()
    return d


ref = read_fa(a.refs, strip_bar=True)
pool_all = read_fa(a.homologs)
print(f"{len(ref)} calibration enzymes, {len(pool_all)} homologs available", file=sys.stderr)

rng = random.Random(a.seed)
keys = sorted(pool_all)
pool = {k: pool_all[k] for k in (keys if len(keys) <= a.n_homologs
                                 else rng.sample(keys, a.n_homologs))}
print(f"using {len(pool)} homologs", file=sys.stderr)

msadir = os.path.join(a.outdir, "msa"); os.makedirs(msadir, exist_ok=True)
yamldir = os.path.join(a.outdir, "yaml"); os.makedirs(yamldir, exist_ok=True)
work = os.path.join(a.outdir, ".work"); os.makedirs(work, exist_ok=True)

cat = os.path.join(work, "all.fa")
with open(cat, "w") as fh:
    for k, v in ref.items():
        fh.write(f">Q_{k}\n{v}\n")
    for k, v in pool.items():
        fh.write(f">H_{k}\n{v}\n")

sto = os.path.join(work, "msa.sto")
subprocess.run(["hmmalign", "--amino", "-o", sto, a.hmm, cat], check=True)
seqs, _ = parse_stockholm(sto)
print(f"aligned {len(seqs)} sequences", file=sys.stderr)

for k, full in ref.items():
    qk = next(x for x in seqs if x.split("|")[0] == f"Q_{k}")
    qa = seqs[qk]
    keep = [i for i, c in enumerate(qa) if c.isalpha()]       # query-length, inserts included
    q = "".join(qa[i] for i in keep).upper()
    assert q == full.upper(), (
        f"{k}: a3m query is {len(q)} aa but the real sequence is {len(full)} aa. "
        "The MSA must be query-length; do not restrict to match states.")
    lines = [f">{k}", q]
    for h in seqs:
        if not h.startswith("H_"):
            continue
        s = seqs[h]
        lines.append(">" + h[2:].split("|")[0])
        lines.append("".join(s[i] if s[i].isalpha() else "-" for i in keep).upper())
    open(os.path.join(msadir, f"{k}.a3m"), "w").write("\n".join(lines) + "\n")

    yml = (f"version: 1\nsequences:\n  - protein:\n      id: A\n"
           f"      sequence: {full}\n"
           f"      msa: {os.path.abspath(msadir)}/{k}.a3m\n"
           f"  - ligand:\n      id: B\n      ccd: {a.ligand_ccd}\n")
    open(os.path.join(yamldir, f"{k}.yaml"), "w").write(yml)

print(f"wrote {len(ref)} a3m + yaml inputs under {a.outdir}", file=sys.stderr)
print("every a3m query row verified equal to the real sequence", file=sys.stderr)
