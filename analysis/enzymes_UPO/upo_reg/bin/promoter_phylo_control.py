#!/usr/bin/env python3.12
"""
G3 caveat 5: is the clade IV CCAAT enrichment descent rather than regulation?

THE PROBLEM
-----------
Clade IV proteins sit at 66.5 % median identity. Promoters of related genes can resemble each
other simply because the genes share an ancestor, with no implication that the shared motif is
used for anything. Cross-lineage replication argues against pure descent, because Sordariomycetes
and Dothideomycetes are different classes, but G3_FINDINGS.md records that no explicit test was
run. This runs it.

THE TEST
--------
Among clade IV promoters, split the proteins by whether FIMO finds the CCAAT motif.

If the motif is present by DESCENT, the carriers are a clade: they should be markedly more
similar to each other than a random subset of the same size. If the motif is distributed across
the clade IV phylogeny without regard to relatedness, descent does not explain it.

Statistic: mean pairwise protein identity among carriers, minus the mean over all clade IV pairs.
Null: permute the carrier labels, holding the number of carriers fixed, and recompute. This
preserves the identity matrix exactly and asks only whether the labels sit where relatedness
predicts.

A one-sided p is reported for carriers being MORE similar than chance, which is the direction
descent predicts.

WHAT A NEGATIVE RESULT DOES AND DOES NOT MEAN
---------------------------------------------
Failing to reject the null does not prove the motif is functional. It removes one specific
alternative explanation. The motif could still be neutral, or enriched for reasons unrelated to
both descent and regulation. Only expression data decides use, and that is Phase B.

Protein identity is used as the distance. It is a proxy for the phylogeny, not the phylogeny.
Across a 66.5 %-identity clade it orders relatedness adequately for this question, but a pair of
long branches can read as distant when the topology says otherwise.
"""
from __future__ import annotations
import argparse, os, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--fg", nargs="+", required=True, help="clade IV foreground promoter FASTAs")
ap.add_argument("--motif-dir", nargs="+", required=True, help="STREME dirs holding the motifs")
ap.add_argument("--aln", required=True, help="protein alignment (hmmalign match states) for reps")
ap.add_argument("--clade-tsv", required=True, help="rep_clade_lineage.tsv")
ap.add_argument("--clade", default="IV")
ap.add_argument("--motif-id", nargs="*", default=None,
                help="explicit STREME motif ids to scan with FIMO, e.g. 3-CCAATRAGAG. These are "
                     "real PWMs; a bare CCAAT 5-mer occurs by chance in ~94 %% of promoters and "
                     "cannot discriminate anything.")
ap.add_argument("--motif-re", default="CCAAT",
                help="fallback literal motif if FIMO finds no STREME motif file")
ap.add_argument("--permutations", type=int, default=10000)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--out", default=None)
a = ap.parse_args()


def read_fa(p):
    d, n = {}, None
    op = subprocess.Popen(["zstd", "-dcq", p], stdout=subprocess.PIPE, text=True).stdout \
        if p.endswith(".zst") else open(p)
    for l in op:
        l = l.rstrip("\n")
        if l.startswith(">"):
            n = l[1:].split()[0]
            d[n] = ""
        elif n is not None:
            d[n] += l.strip()
    return d


# ---- clade IV membership -------------------------------------------------
clade = {}
with open(a.clade_tsv) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    i = {k: j for j, k in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > i["clade"]:
            clade[f[i["gene_id"]]] = (f[i["clade"]], f[i["protein_id"]],
                                      f[i.get("LINEAGE_GROUP", 0)])
members = {g for g, v in clade.items() if v[0] == a.clade}
print(f"{len(members)} clade {a.clade} representatives", file=sys.stderr)

# ---- promoters -----------------------------------------------------------
prom = {}
for p in a.fg:
    for k, v in read_fa(p).items():
        g = k.split("|")[0]
        if g in members:
            prom[g] = v
print(f"{len(prom)} clade {a.clade} promoters found in the foreground FASTAs", file=sys.stderr)
if not prom:
    sys.exit("no promoters matched; check --fg and --clade-tsv id formats")

work = os.path.join(os.path.dirname(os.path.abspath(a.out or ".")), ".phylo_ctrl")
os.makedirs(work, exist_ok=True)
pf = os.path.join(work, "cladeIV.fa")
with open(pf, "w") as fh:
    for g, s in sorted(prom.items()):
        fh.write(f">{g}\n{s}\n")

# ---- motif presence ------------------------------------------------------
carriers = set()
used = None
mfile = None
for d in a.motif_dir:
    c = os.path.join(d, "streme.xml")
    if os.path.exists(c):
        mfile = c
        break
if mfile:
    fo = os.path.join(work, "fimo")
    cmd = ["fimo", "--oc", fo, "--thresh", "1e-4"]
    for mid in (a.motif_id or []):
        cmd += ["--motif", mid]
    cmd += [mfile, pf]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        print("fimo not on PATH; falling back to the literal two-strand scan", file=sys.stderr)
        r = None
    tsv = os.path.join(fo, "fimo.tsv")
    if os.path.exists(tsv):
        hits = {}
        with open(tsv) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ix = {k: j for j, k in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or line.startswith("#") or not line.strip():
                    continue
                # When explicit motif ids were requested, FIMO already restricted the scan,
                # so every row counts. Otherwise keep the CCAAT-bearing ones.
                mid = f[ix.get("motif_id", 0)]
                keep = bool(a.motif_id) or "CCAAT" in mid.upper() or "ATTGG" in mid.upper()
                if keep:
                    hits.setdefault(f[ix["sequence_name"]], 0)
                    hits[f[ix["sequence_name"]]] += 1
        if hits:
            carriers = set(hits)
            used = f"FIMO on {mfile}, p < 1e-4, motifs whose id or match contains CCAAT"
if len(carriers) < 20 and not a.motif_id:
    # A CCAAT box is bound on either strand: the Eurotiomycetes STREME motif reports it as
    # ATTGGG, the reverse complement. Scanning one strand undercounts it badly - the FIMO pass
    # above found 3 of 760 that way. Scan both strands.
    rc = str.maketrans("ACGT", "TGCA")
    m = a.motif_re.upper()
    mrc = m.translate(rc)[::-1]
    carriers = {g for g, s in prom.items() if m in s.upper() or mrc in s.upper()}
    used = f"literal match of '{m}' or its reverse complement '{mrc}', both strands"
print(f"motif presence by: {used}", file=sys.stderr)
print(f"{len(carriers)} of {len(prom)} clade {a.clade} promoters carry it "
      f"({100*len(carriers)/len(prom):.1f} %)", file=sys.stderr)

# ---- protein identity ----------------------------------------------------
aln = read_fa(a.aln)
pid2gene = {v[1]: g for g, v in clade.items()}
rows, names = [], []
for k, s in aln.items():
    base = k.split("|")[0]
    g = pid2gene.get(base) or (base if base in prom else None)
    if g in prom:
        rows.append(s.upper())
        names.append(g)
print(f"{len(names)} of {len(prom)} have an alignment row", file=sys.stderr)
if len(names) < 20:
    sys.exit("too few aligned clade IV proteins to test")

L = min(len(r) for r in rows)
M = np.array([[c for c in r[:L]] for r in rows])
gap = (M == "-") | (M == ".")
n = len(names)
ident = np.zeros((n, n))
for i in range(n):
    same = (M[i] == M) & ~gap[i] & ~gap
    comp = ~gap[i] & ~gap
    cs = comp.sum(1)
    ident[i] = np.where(cs > 0, same.sum(1) / np.maximum(cs, 1), np.nan)
np.fill_diagonal(ident, np.nan)
print(f"median pairwise identity across clade {a.clade}: "
      f"{100*np.nanmedian(ident):.1f} %", file=sys.stderr)

lab = np.array([1 if g in carriers else 0 for g in names])
k = int(lab.sum())
print(f"carriers among aligned: {k} of {n}", file=sys.stderr)
if k < 5 or n - k < 5:
    sys.exit(f"too few in one group (carriers {k}, non-carriers {n-k}) to permute")


def mean_within(mask):
    idx = np.where(mask)[0]
    sub = ident[np.ix_(idx, idx)]
    return float(np.nanmean(sub))


obs_in = mean_within(lab == 1)
overall = float(np.nanmean(ident))
obs = obs_in - overall
rng = np.random.default_rng(a.seed)
null = np.empty(a.permutations)
base = np.zeros(n, dtype=bool)
for t in range(a.permutations):
    p = rng.permutation(n)[:k]
    m = base.copy(); m[p] = True
    null[t] = mean_within(m) - overall
p_more = (np.sum(null >= obs) + 1) / (a.permutations + 1)

print("\n=== PROMOTER-PHYLOGENY CONTROL ===", file=sys.stderr)
print(f"mean pairwise identity, carriers only : {100*obs_in:.2f} %", file=sys.stderr)
print(f"mean pairwise identity, all clade {a.clade}   : {100*overall:.2f} %", file=sys.stderr)
print(f"difference                            : {100*obs:+.2f} points", file=sys.stderr)
print(f"null difference, mean of {a.permutations} permutations: "
      f"{100*null.mean():+.3f} points (sd {100*null.std():.3f})", file=sys.stderr)
print(f"one-sided p (carriers MORE similar)   : {p_more:.4f}", file=sys.stderr)
z = (obs - null.mean()) / null.std() if null.std() else float("nan")
print(f"z = {z:+.2f}", file=sys.stderr)
print("\nDescent predicts carriers are markedly more similar than chance.", file=sys.stderr)
if p_more > 0.05:
    print("NOT SUPPORTED here: carriers are not significantly more related than a random",
          file=sys.stderr)
    print("subset of the same size. Descent does not explain the motif's distribution.",
          file=sys.stderr)
else:
    print("SUPPORTED: carriers ARE more closely related than chance. The enrichment is",
          file=sys.stderr)
    print("confounded with ancestry and must not be read as convergent regulation.",
          file=sys.stderr)

if a.out:
    with open(a.out, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"clade\t{a.clade}\n")
        fh.write(f"motif_presence_method\t{used}\n")
        fh.write(f"n_aligned\t{n}\ncarriers\t{k}\n")
        fh.write(f"mean_identity_carriers\t{obs_in:.6f}\n")
        fh.write(f"mean_identity_all\t{overall:.6f}\n")
        fh.write(f"difference_points\t{100*obs:.4f}\n")
        fh.write(f"null_mean_points\t{100*null.mean():.4f}\n")
        fh.write(f"null_sd_points\t{100*null.std():.4f}\n")
        fh.write(f"p_one_sided_more_similar\t{p_more:.6f}\nz\t{z:.4f}\n")
        fh.write(f"permutations\t{a.permutations}\n")
    print(f"\nwrote {a.out}", file=sys.stderr)
