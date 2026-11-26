#!/bin/bash
# Joint SDR tree inputs (t4hn_reductase + t3hn_reductase). Run from analysis/pathway_census after make_leaves.sh. Needs module hmmer/3.4.
I=results/step5/other_inputs; T=$I/tree_inputs/sdr_joint; L=$I/leaves; H=results/hmms/melanin_dhn/t4hn_reductase.PROVISIONAL.hmm
mkdir -p $T results/search/sdr_joint
/usr/bin/python3.12 - <<'PY'
import glob, subprocess, os
out = "results/search/sdr_joint"
n = {}
for p in sorted(glob.glob("results/search/melanin_myco_full/candidates.chunk*.tsv.zst")):
    lines = subprocess.run(["zstd", "-dc", p], capture_output=True, check=True).stdout.decode().splitlines()
    hdr = [l for l in lines if not l.startswith("#")][0]
    keep = []
    for l in lines:
        if l.startswith("#") or l == hdr:
            continue
        f = l.split("\t")
        if f[-1] in ("t4hn_reductase", "t3hn_reductase"):
            n[f[-1]] = n.get(f[-1], 0) + 1
            f[-1] = "t4hn_reductase"        # relabelled: one candidate set for the joint tree
            keep.append("\t".join(f))
    txt = "# PROVISIONAL\n" + hdr + "\n" + "".join(k + "\n" for k in keep)
    subprocess.run(["zstd", "-q", "-f", "-o", os.path.join(out, os.path.basename(p))], input=txt.encode(), check=True)
print("rows relabelled:", n)
PY
cat $L/t4hn_reductase.named.faa $L/t3hn_reductase.named.faa > $T/named_union_raw.faa
cat $L/t4hn_reductase.outgroups.faa $L/t3hn_reductase.outgroups.faa > $T/og_union_raw.faa
cat $I/decoys/t4hn_reductase.best.faa $I/decoys/t3hn_reductase.best.faa > $T/decoys_union_raw.faa
/usr/bin/python3.12 - <<'PY'
import csv
I = "results/step5/other_inputs"; T = I + "/tree_inputs/sdr_joint"
def dedupe(src, dst):
    seen = set(); keep = False
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                n = line[1:].split()[0]; keep = n not in seen; seen.add(n)
            if keep:
                fo.write(line)
    return len(seen)
floor = {"t4hn_reductase": 356.0, "t3hn_reductase": 550.4}   # a leaf at or above either floor on either HMM is dropped
bad = set()
for s, f in floor.items():
    for r in csv.DictReader(open(f"{I}/screen/{s}.family_screen.tsv"), delimiter="\t"):
        if r["best_score"] and float(r["best_score"]) >= f:
            bad.add(r["protein"].split("|")[1] if "|" in r["protein"] else r["protein"])
open(T + "/floor_excluded.txt", "w").write("\n".join(sorted(bad)) + "\n")
for a, b in (("named_union_raw", "named_union"), ("og_union_raw", "og_union"), ("decoys_union_raw", "decoys_union")):
    print(a, dedupe(f"{T}/{a}.faa", f"{T}/{b}.faa"))
PY
/usr/bin/python3.12 bin/cut_leaves_to_domain.py --fasta $T/decoys_union.faa --hmm $H --out-fasta $T/decoys.faa --out-table $T/decoys_cut.tsv | tail -1
grep '>' $T/decoys.faa | sed 's/>//' | awk '{print $1}' > $T/decoy_ids.txt
/usr/bin/python3.12 - <<'PY'
T = "results/step5/other_inputs/tree_inputs/sdr_joint"
drop = {l.strip() for l in open(T + "/decoy_ids.txt")} | {l.strip() for l in open(T + "/floor_excluded.txt")}
def filt(src, dst):
    keep = False; n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                nm = line[1:].split()[0]; a = nm.split("|")[1] if "|" in nm else nm
                keep = a not in drop and nm not in drop; n += keep
            if keep:
                fo.write(line)
    print(dst, n)
filt(T + "/named_union.faa", T + "/named_f.faa"); filt(T + "/og_union.faa", T + "/og_f.faa")
PY
/usr/bin/python3.12 bin/cut_leaves_to_domain.py --fasta $T/named_f.faa --hmm $H --out-fasta $T/members_cut.faa --out-table $T/members_cut.tsv | tail -1
/usr/bin/python3.12 bin/cut_leaves_to_domain.py --fasta $T/og_f.faa --hmm $H --out-fasta $T/outgroups.faa --out-table $T/outgroups_cut.tsv | tail -1
cat $T/members_cut.faa refs/melanin_dhn/t3hn_reductase.faa > $T/members.faa
cat $T/members_cut.faa refs/melanin_dhn/t4hn_reductase.faa > $T/members_for_t3hn.faa
grep -c '>' $T/members.faa $T/members_for_t3hn.faa $T/outgroups.faa $T/decoys.faa
