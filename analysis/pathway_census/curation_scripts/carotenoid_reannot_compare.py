#!/usr/bin/python3.12
"""Compare carotenoid candidate rows of the 5 re-annotated genomes: old full search vs re-search.

Usage: carotenoid_reannot_compare.py OLD_CAND_DIR NEW_CAND_DIR ALIASES SEQDB_OLD_DIR OUT_DIR
A candidate row is (protein, step). Protein ids are not comparable between annotations
(new numbering), so a row is matched by exact protein SEQUENCE: old sequence from the
2026-11-06 seqdb snapshot, new sequence from the current predict file (path in ALIASES).
Per (genome, step): unchanged = sequence is a candidate in both; lost = old only; gained = new only.
"""
import csv, glob, gzip, hashlib, os, subprocess, sys
from collections import defaultdict

old_dir, new_dir, alias_f, seqdb, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
al = list(csv.DictReader(open(alias_f), delimiter="\t"))
old2new = {r["locustag"]: r["alias_prefix"] for r in al}
new2old = {v: k for k, v in old2new.items()}
info = {r["locustag"]: r for r in al}

def cand_rows(path):
    p = subprocess.run(["zstd", "-dc", path], capture_output=True, text=True, check=True)
    return list(csv.DictReader([l for l in p.stdout.splitlines() if not l.startswith("#")], delimiter="\t"))

def fasta(path):
    d, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith(">"):
                if name: d[name] = "".join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    if name: d[name] = "".join(buf)
    return d

old_seq = {g: fasta(f"{seqdb}/{hashlib.sha1(g.encode()).hexdigest()[:2]}/{g}.proteins.fa.gz") for g in old2new}
new_seq = {g: fasta(info[g]["proteins_path"]) for g in old2new}
for g in old2new:
    print(g, "old proteins", len(old_seq[g]), "new proteins", len(new_seq[g]),
          "identical sequences shared", len(set(old_seq[g].values()) & set(new_seq[g].values())))

old = defaultdict(dict); new = defaultdict(dict)  # (g, step) -> {seq: (protein_id, score)}
for f in sorted(glob.glob(f"{old_dir}/candidates.chunk*.zst")):
    for r in cand_rows(f):
        g = r["protein_id"].split("_")[0]
        if g in old2new:
            old[(g, r["step_id"])][old_seq[g][r["protein_id"]]] = (r["protein_id"], float(r["score"]))
for f in sorted(glob.glob(f"{new_dir}/candidates.chunk*.zst")):
    for r in cand_rows(f):
        g = new2old[r["protein_id"].split("_")[0]]
        new[(g, r["step_id"])][new_seq[g][r["protein_id"]]] = (r["protein_id"], float(r["score"]))

steps = sorted({s for _, s in list(old) + list(new)})
with open(f"{out}/rows_per_genome_step.tsv", "w") as o, open(f"{out}/changed_rows.tsv", "w") as c:
    o.write("old_prefix\tnew_prefix\tspecies\tstep\told_rows\tnew_rows\tunchanged\tlost\tgained\tpresence_old\tpresence_new\tpresence_changes\n")
    c.write("old_prefix\tnew_prefix\tstep\tstatus\told_protein\tnew_protein\told_score\tnew_score\n")
    for g in old2new:
        for s in steps:
            a, b = old.get((g, s), {}), new.get((g, s), {})
            if not a and not b: continue
            same, lost, gain = set(a) & set(b), set(a) - set(b), set(b) - set(a)
            po, pn = int(bool(a)), int(bool(b))
            o.write(f"{g}\t{old2new[g]}\t{info[g]['species']}\t{s}\t{len(a)}\t{len(b)}\t{len(same)}\t{len(lost)}\t{len(gain)}\t{po}\t{pn}\t{int(po!=pn)}\n")
            for q in lost: c.write(f"{g}\t{old2new[g]}\t{s}\tlost\t{a[q][0]}\t\t{a[q][1]}\t\n")
            for q in gain: c.write(f"{g}\t{old2new[g]}\t{s}\tgained\t\t{b[q][0]}\t\t{b[q][1]}\n")
            for q in same:
                if a[q][1] != b[q][1]:
                    c.write(f"{g}\t{old2new[g]}\t{s}\tscore_differs\t{a[q][0]}\t{b[q][0]}\t{a[q][1]}\t{b[q][1]}\n")
