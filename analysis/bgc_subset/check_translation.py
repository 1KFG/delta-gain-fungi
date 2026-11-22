#!/usr/bin/env python3
"""Compare antiSMASH CDS translations (from a raw JSON) with the BFD protein FASTA of the same genome.
Usage: check_translation.py antismash.json.gz proteins.fa.gz   -- prints matched/mismatched counts and a few mismatches."""
import gzip, json, sys
js, fa = sys.argv[1:3]
P = {}; name = None
for l in gzip.open(fa, "rt"):
    if l.startswith(">"): name = l[1:].split()[0]; P[name] = []
    else: P[name].append(l.strip())
P = {k: "".join(v).rstrip("*") for k, v in P.items()}
d = json.load(gzip.open(js, "rt"))
n = ok = miss = 0; bad = []
for rec in d["records"]:
    for f in rec["features"]:
        if f["type"] != "CDS": continue
        q = f["qualifiers"]; gid = (q.get("locus_tag") or q.get("gene") or [""])[0]
        tr = (q.get("translation") or [""])[0] if isinstance(q.get("translation"), list) else q.get("translation", "")
        cand = [k for k in (gid + "-T1", gid) if k in P]
        if not cand: miss += 1; continue
        n += 1
        if P[cand[0]] == tr.rstrip("*"): ok += 1
        else: bad.append((gid, len(P[cand[0]]), len(tr)))
print(f"CDS compared {n}; identical {ok}; different {n - ok}; no BFD protein found {miss}")
for b in bad[:5]: print("  different:", b)
