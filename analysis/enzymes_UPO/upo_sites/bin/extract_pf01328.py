#!/usr/bin/env python3.12
"""Pull PF01328 protein sequences out of the BFD seqdb into chunked FASTA files.

Uses the bgzip/faidx random-access index that ships with seqdb, via `samtools faidx`, with one
call per genome carrying every wanted id. That is far faster than decompressing whole proteomes:
measured on this dataset, streaming every file projected to ~2 h, indexed access takes minutes.

Reports every LOCUSTAG it could not open and every protein id it could not retrieve, rather than
silently dropping them.
"""
import argparse, hashlib, os, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--ids", required=True, help="TSV: LOCUSTAG<TAB>protein_id, with header")
ap.add_argument("--seqdb", required=True, help=".../seqdb/proteins")
ap.add_argument("--outdir", required=True)
ap.add_argument("--per-chunk", type=int, default=5000)
ap.add_argument("--samtools", default="samtools")
a = ap.parse_args()

want = defaultdict(list)
with open(a.ids) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            want[f[0]].append(f[1])
total = sum(len(v) for v in want.values())
print(f"want {total} proteins from {len(want)} genomes", file=sys.stderr)

os.makedirs(a.outdir, exist_ok=True)
state = {"i": 0, "n": 0, "fh": None}
def rotate():
    if state["fh"]:
        state["fh"].close()
    state["i"] += 1; state["n"] = 0
    state["fh"] = open(os.path.join(a.outdir, f"chunk_{state['i']:04d}.fa"), "w")
rotate()

def emit(name: str, seq: str):
    state["fh"].write(f">{name}\n{seq}\n")
    state["n"] += 1
    if state["n"] >= a.per_chunk:
        rotate()

found = 0
bad_files, missing_ids = [], []
for n, (lt, ids) in enumerate(sorted(want.items()), 1):
    sub = hashlib.sha1(lt.encode()).hexdigest()[:2]
    path = os.path.join(a.seqdb, sub, f"{lt}.proteins.fa.gz")
    if not (os.path.exists(path) and os.path.exists(path + ".fai")):
        bad_files.append(lt); continue
    try:
        r = subprocess.run([a.samtools, "faidx", path] + ids,
                           capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        bad_files.append(f"{lt}:timeout"); continue
    if r.returncode != 0:
        bad_files.append(f"{lt}:rc={r.returncode}:{r.stderr.strip()[:120]}"); continue
    got, cur, buf = set(), None, []
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            if cur:
                emit(cur, "".join(buf)); found += 1; got.add(cur)
            cur, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    if cur:
        emit(cur, "".join(buf)); found += 1; got.add(cur)
    missing_ids.extend(sorted(set(ids) - got))
    if n % 2000 == 0:
        print(f"  {n}/{len(want)} genomes, {found} proteins", file=sys.stderr)
state["fh"].close()

print(f"\nwrote {found}/{total} proteins into {state['i']} chunks in {a.outdir}", file=sys.stderr)
print(f"genome files missing/unreadable/failed: {len(bad_files)}", file=sys.stderr)
print(f"protein ids not retrieved: {len(missing_ids)}", file=sys.stderr)
if bad_files:
    open(os.path.join(a.outdir, "MISSING_genomes.txt"), "w").write("\n".join(bad_files) + "\n")
if missing_ids:
    open(os.path.join(a.outdir, "MISSING_proteins.txt"), "w").write("\n".join(missing_ids) + "\n")
