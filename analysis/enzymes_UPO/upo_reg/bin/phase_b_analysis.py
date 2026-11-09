#!/usr/bin/env python3.12
"""
G3 Phase B: do UPOs respond to oxidative stress or iron, and are they co-regulated with the
CCAAT-binding complex?

Design and pre-registered interpretation rules are in G3_PHASE_B_DESIGN.md, committed before
this ran. Read that first. The short version of the limits:

  - The reference genome Y9 has 13 PF01328 proteins and exactly ONE clade IV by placement, so the
    clade-IV-specific question is an observation on 1 gene, never a test. No p-value is reported
    for clade IV as a clade.
  - Reads are strain Z7, transcripts are strain Y9. Mapping rate is the QC and is reported first.
  - NaCl and CuSO4 are specificity controls. A response shared with them is generic.

TESTS
-----
1. CLASS LEVEL (powered). Are the 13 UPO genes shifted in log2 fold change relative to all other
   genes, under each treatment? Mann-Whitney U, two-sided, against the genome-wide distribution.
2. CLADE IV (descriptive only). The single clade IV gene's fold change, shown beside the other 12.
3. CO-REGULATION. Spearman correlation of each UPO's expression profile with the CCAAT-complex
   subunits across all 15 libraries. A correlation across 15 samples is a number; whether it
   means co-regulation is a separate question, and shared response to the same treatments will
   produce correlation without any regulatory link.
"""
from __future__ import annotations
import argparse, glob, json, os, sys
from collections import defaultdict
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--quant-dir", required=True)
ap.add_argument("--manifest", required=True, help="runs.tsv: run, condition, reads, project")
ap.add_argument("--upo-list", required=True, help="file of UPO gene ids in this genome")
ap.add_argument("--clade-tsv", required=True, help="per-protein table (scored) for clade labels")
ap.add_argument("--ccaat", nargs="*", default=[], help="CCAAT-complex gene ids")
ap.add_argument("--min-maprate", type=float, default=50.0)
ap.add_argument("--control", default="control")
ap.add_argument("--report", default=None)
a = ap.parse_args()

OUT = []
def say(s=""):
    print(s, file=sys.stderr)
    OUT.append(s)

# ---- manifest ----
cond = {}
for line in open(a.manifest):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 2:
        cond[f[0]] = f[1]

# ---- mapping-rate QC ----
say("# G3 Phase B results\n")
say("## Mapping-rate QC (strain Z7 reads against strain Y9 transcripts)\n")
say("| run | condition | % mapped |")
say("|---|---|---|")
rates, keep = {}, []
for r in sorted(cond):
    mj = os.path.join(a.quant_dir, r, "aux_info", "meta_info.json")
    if not os.path.exists(mj):
        say(f"| {r} | {cond[r]} | MISSING |")
        continue
    m = json.load(open(mj)).get("percent_mapped")
    rates[r] = m
    say(f"| {r} | {cond[r]} | {m:.1f} |")
    if m is not None and m >= a.min_maprate:
        keep.append(r)
if not rates:
    sys.exit("no quantifications found")
say(f"\nMedian mapping rate {np.median(list(rates.values())):.1f} %. "
    f"{len(keep)} of {len(rates)} runs pass the {a.min_maprate:.0f} % floor.\n")
if len(keep) < 10:
    say("**Too few runs pass the mapping floor. Cross-strain quantification has failed and the "
        "comparison is void.** Reporting this rather than working around it.")
    if a.report:
        open(a.report, "w").write("\n".join(OUT) + "\n")
    sys.exit(0)

# ---- expression matrix ----
tpm = defaultdict(dict)
for r in keep:
    p = os.path.join(a.quant_dir, r, "quant.sf")
    with open(p) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tpm[f[0].split()[0]][r] = float(f[3])
genes = sorted(tpm)
say(f"{len(genes)} transcripts quantified across {len(keep)} libraries.\n")


def norm_gene(t):
    """salmon transcript id -> gene id (strip -T1/-mRNA suffixes)."""
    g = t.split()[0]
    for suf in ("-T1", "-T2", "-T3"):
        if g.endswith(suf):
            return g[: -len(suf)]
    return g


upos = {l.strip() for l in open(a.upo_list) if l.strip()}
upos = {norm_gene(u) for u in upos}
ccaat = {norm_gene(c) for c in a.ccaat}

# clade labels
clade = {}
import subprocess
txt = subprocess.run(["zstd", "-dcq", a.clade_tsv], capture_output=True, text=True).stdout
hdr = None
for line in txt.splitlines():
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if hdr is None:
        hdr = f; ix = {k: j for j, k in enumerate(hdr)}; continue
    g = norm_gene(f[ix["protein_id"]])
    if g in upos:
        clade[g] = (f[ix["clade"]], f[ix["clade_trust"]], f[ix["t1_pocket_matches"]])

gmap = {}
for t in genes:
    gmap.setdefault(norm_gene(t), []).append(t)
mat = {}
for g, ts in gmap.items():
    mat[g] = np.array([sum(tpm[t].get(r, 0.0) for t in ts) for r in keep])
found = [g for g in upos if g in mat]
say(f"{len(found)} of {len(upos)} UPO genes present in the quantification.\n")

conds = [cond[r] for r in keep]
uc = sorted(set(conds))
ctrl_idx = [i for i, c in enumerate(conds) if c == a.control]
if not ctrl_idx:
    sys.exit(f"no control libraries labelled '{a.control}'")

L = {g: np.log2(v + 1.0) for g, v in mat.items()}
allg = [g for g in mat if mat[g].sum() > 0]


def mwu(x, y):
    """Two-sided Mann-Whitney U with normal approximation and tie correction."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n1, n2 = len(x), len(y)
    if n1 < 3 or n2 < 3:
        return float("nan"), float("nan")
    allv = np.concatenate([x, y])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv), float)
    sv = allv[order]
    i = 0
    tie = 0.0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        rr = 0.5 * (i + j) + 1.0
        ranks[order[i:j + 1]] = rr
        t = j - i + 1
        if t > 1:
            tie += t ** 3 - t
        i = j + 1
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    N = n1 + n2
    sd = np.sqrt(n1 * n2 / 12.0 * ((N + 1) - tie / (N * (N - 1))))
    if sd == 0:
        return U1, float("nan")
    z = (U1 - mu) / sd
    from math import erfc, sqrt
    return U1, erfc(abs(z) / sqrt(2))


say("## Test 1: are the UPOs as a class shifted, versus all other genes?\n")
say("Mann-Whitney U on per-gene log2 fold change, UPO set against every other expressed gene.\n")
say("| treatment | UPO median log2FC | genome median log2FC | p (two-sided) |")
say("|---|---|---|---|")
fc_store = {}
for c in uc:
    if c == a.control:
        continue
    ti = [i for i, x in enumerate(conds) if x == c]
    fc = {}
    for g in allg:
        fc[g] = float(np.mean(L[g][ti]) - np.mean(L[g][ctrl_idx]))
    fc_store[c] = fc
    u = [fc[g] for g in found if g in fc]
    o = [fc[g] for g in allg if g not in upos]
    _, p = mwu(u, o)
    say(f"| {c} | {np.median(u):+.3f} | {np.median(o):+.3f} | {p:.4g} |")

say("\n**Pre-registered rule:** a response under H2O2 that is also present under NaCl is a "
    "generic stress response. A response under FeSO4 also present under CuSO4 is a generic "
    "metal response.\n")

say("## Test 2: the individual UPO genes (descriptive)\n")
hdrs = ["gene", "clade", "trust", "pocket"] + [c for c in uc if c != a.control]
say("| " + " | ".join(hdrs) + " |")
say("|" + "---|" * len(hdrs))
for g in sorted(found, key=lambda x: -(fc_store.get("H2O2", {}).get(x, 0))):
    cl, tr, pk = clade.get(g, ("?", "?", "?"))
    row = [g, cl, tr, pk] + [f"{fc_store[c].get(g, float('nan')):+.2f}"
                             for c in uc if c != a.control]
    say("| " + " | ".join(row) + " |")
say("\nThe clade IV row is ONE gene. It is shown, not tested. No claim about clade IV as a "
    "clade can rest on it.\n")

say("## Test 3: co-regulation with the CCAAT-binding complex\n")
cc = [g for g in ccaat if g in mat]
if not cc:
    say("No CCAAT-complex genes found in the quantification; test not run.\n")
else:
    def spearman(u, v):
        ru = np.argsort(np.argsort(u)).astype(float)
        rv = np.argsort(np.argsort(v)).astype(float)
        ru -= ru.mean(); rv -= rv.mean()
        d = np.sqrt((ru ** 2).sum() * (rv ** 2).sum())
        return float((ru * rv).sum() / d) if d else float("nan")

    say(f"{len(cc)} CCAAT-complex genes. Spearman across {len(keep)} libraries.\n")
    say("| UPO gene | clade | mean rho vs CCAAT | mean rho vs random genes |")
    say("|---|---|---|---|")
    rng = np.random.default_rng(1337)
    bg = rng.choice([g for g in allg if g not in upos and g not in ccaat],
                    size=min(300, len(allg)), replace=False)
    rows = []
    for g in sorted(found):
        if mat[g].sum() == 0:
            continue
        rc = float(np.mean([spearman(L[g], L[c]) for c in cc]))
        rb = float(np.mean([spearman(L[g], L[b]) for b in bg]))
        cl = clade.get(g, ("?",))[0]
        rows.append((g, cl, rc, rb))
        say(f"| {g} | {cl} | {rc:+.3f} | {rb:+.3f} |")
    if rows:
        u_rc = [r[2] for r in rows]; u_rb = [r[3] for r in rows]
        say(f"\nUPO mean rho vs CCAAT {np.mean(u_rc):+.3f}, vs random {np.mean(u_rb):+.3f}.")
        say("\nCorrelation across 15 libraries that differ by treatment is produced by any shared "
            "response to those treatments. It is not evidence of a regulatory link, and the "
            "difference from the random-gene baseline is the only part worth reading.\n")

if a.report:
    open(a.report, "w").write("\n".join(OUT) + "\n")
    print(f"\nwrote {a.report}", file=sys.stderr)
