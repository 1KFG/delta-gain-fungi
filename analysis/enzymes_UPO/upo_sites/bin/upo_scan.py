#!/usr/bin/env python3.12
"""
upo_scan: one command from a proteome FASTA to a scored, domain-bounded UPO site table.

WHAT IT DOES
------------
  1. finds PF01328 proteins and reads the active-site residues   (upo_sites.py)
  2. scores the aliphatic pocket against the active consensus    (the T1 rule)
  3. states, per protein, whether it lies inside the range of active-site divergence that the
     19 characterized enzymes actually span                      (the applicability domain)

Step 3 is not decoration. Without it the score is a number with no stated scope, and 69.6 % of
BFD lies outside that scope.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
No phylogenetic placement. Placement needs the reference tree and RAxML EPA, it declines 78.9 %
of inputs, and it is a separate, slower stage. Run `make_placement_aln.py` + `jplace_to_clade.py`
if a clade call is needed. Without placement there is no `clade` column here, and the score's
central limit - that it detects clade-V-LIKE pockets - cannot be checked per protein. The header
of the output says so.

READ BEFORE USING THE SCORE COLUMN
-----------------------------------
The score is cross-validated (LOO R2 +0.496, AUC 0.956, permutation p 0.0005, n = 19) but 9 of
the 10 active calibration enzymes are clade V. With clade V removed the AUC is 0.556, which is
chance, on ONE remaining active. The score detects clade-V-like aliphatic pockets. It is NOT
shown to find actives outside clade V, and it MISSES the one that exists (GmaUPO-II, 49.3 %
active, 1 of 5 matches). A high score is a hypothesis worth assaying. A LOW SCORE IS NOT
EVIDENCE OF INACTIVITY. Substrate is 4-propylguaiacol conversion only; this says nothing about
which rings an enzyme cleaves.

MODEL_VALIDATION.md sections D, E and F hold the measurements behind every sentence above.
"""
from __future__ import annotations
import argparse, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
POCKET5 = (73, 77, 81, 183, 187)
DIST_POS = (80, 84, 108, 110, 172, 173, 175, 176, 179, 180)
CONS = "ILLAI"
MAXD = 0.50          # largest nearest-neighbour distance WITHIN the 19 (applicability.py)

ap = argparse.ArgumentParser(
    description="Proteome FASTA -> scored, domain-bounded UPO active-site table.")
ap.add_argument("--proteome", required=True, nargs="+", help="FASTA (plain, .gz or .zst)")
ap.add_argument("--hmm", required=True, help="PF01328.hmm")
ap.add_argument("--backbone", required=True, help="reference alignment (.afa/.afa.zst)")
ap.add_argument("--s2", required=True, help="active_site_S2.tsv, the 19 characterized enzymes")
ap.add_argument("--out", required=True, help="output TSV (.zst recommended)")
ap.add_argument("--anchor", default="HspUPO")
ap.add_argument("--cpu", type=int, default=4)
ap.add_argument("--cut", type=int, default=3, help="pocket matches at or above which to flag")
ap.add_argument("--sites-tsv", default=None,
                help="reuse an existing upo_sites.py table instead of recomputing")
ap.add_argument("--keep-sites", default=None, help="also write the raw site layer here")
a = ap.parse_args()


def log(m):
    print(f"[upo_scan] {m}", file=sys.stderr)


# ---- 1. site layer ----
sites = a.sites_tsv
tmp = None
if not sites:
    tmp = a.keep_sites or os.path.join(
        os.environ.get("SCRATCH", os.environ.get("TMPDIR", "/tmp")),
        f"upo_sites_{os.getpid()}.tsv.zst")
    cmd = [sys.executable, "-P", os.path.join(HERE, "upo_sites.py"),
           "--proteome", *a.proteome, "--hmm", a.hmm, "--backbone", a.backbone,
           "--anchor", a.anchor, "--out", tmp, "--cpu", str(a.cpu)]
    log("step 1/3: finding PF01328 proteins and reading active-site residues")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(f"upo_sites.py failed with exit {r.returncode}")
    sites = tmp
else:
    log(f"step 1/3: reusing {sites}")


def read_rows(p):
    if p.endswith(".zst"):
        it = subprocess.Popen(["zstd", "-dcq", p], stdout=subprocess.PIPE, text=True).stdout
    elif p.endswith(".gz"):
        import gzip
        it = gzip.open(p, "rt")
    else:
        it = open(p)
    hdr = None
    for line in it:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            yield hdr, None
            continue
        yield hdr, f


# ---- 2. calibration profiles for the domain ----
cal = []
hdr = None
for line in open(a.s2):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        ix = {k: j for j, k in enumerate(hdr)}
        continue
    if f[0] == "LOCATION":
        continue
    cal.append([f[ix[f"p{p}"]].strip().upper() if f"p{p}" in ix else "" for p in DIST_POS])
log(f"step 2/3: {len(cal)} calibration profiles over positions {DIST_POS}")
if not cal:
    sys.exit("no calibration profiles read from --s2")


def dist(p, q):
    num = den = 0
    for x, y in zip(p, q):
        if not x or x in "-." or not y or y in "-.":
            continue
        den += 1
        num += (x != y)
    return (num / den) if den else None


# ---- 3. score + domain ----
log("step 3/3: scoring the pocket and bounding the applicability domain")
out_rows = []
cache = {}
n = flagged = outside = 0
for hdr, f in read_rows(sites):
    if f is None:
        ix = {k: j for j, k in enumerate(hdr)}
        pcol = {}
        for p in DIST_POS:
            hit = [j for j, k in enumerate(hdr) if k.split("_")[0] == f"p{p}"]
            if hit:
                pcol[p] = hit[0]
        miss = [p for p in DIST_POS if p not in pcol]
        if miss:
            sys.exit(f"site table lacks active-site columns for positions {miss}")
        if "pocket_residues" not in ix:
            sys.exit("site table lacks a pocket_residues column")
        out_rows.append(hdr + ["t1_pocket_matches", "t1_pocket_flag",
                               "calib_distance", "calib_domain"])
        continue
    n += 1
    pr = f[ix["pocket_residues"]].strip().upper()
    if len(pr) != 5:
        m, flag = -1, "no_pocket"
    else:
        m = sum(pr[k] == CONS[k] for k in range(5))
        flag = "clade_V_like" if m >= a.cut else "below_cut"
        if m >= a.cut:
            flagged += 1
    aa = tuple(f[pcol[p]].strip().upper() if pcol[p] < len(f) else "" for p in DIST_POS)
    if aa not in cache:
        ds = [d for d in (dist(aa, c) for c in cal) if d is not None]
        cache[aa] = min(ds) if ds else None
    d = cache[aa]
    if d is None:
        dom, ds = "unknown", ""
    else:
        dom, ds = ("inside" if d <= MAXD else "outside"), f"{d:.3f}"
        if dom == "outside":
            outside += 1
    out_rows.append(f + [str(m), flag, ds, dom])

head = [
    "# upo_scan output: PF01328 active-site residues, T1 pocket score, applicability domain.",
    f"# t1_pocket_matches: matches (of 5) to the ACTIVE consensus {CONS} at HspUPO "
    f"p73,p77,p81,p183,p187.",
    f"# t1_pocket_flag: 'clade_V_like' at >= {a.cut} matches.",
    "# CROSS-VALIDATED on 19 enzymes: LOO R2 +0.496, LOO AUC 0.956, permutation p 0.0005.",
    "# LIMIT: 9 of 10 active calibration enzymes are clade V. With clade V removed the AUC is",
    "#   0.556 (chance) on ONE active. This detects clade-V-like pockets. It is NOT shown to find",
    "#   actives outside clade V, and it MISSES the one that exists (GmaUPO-II, 49.3 % active,",
    "#   1/5 matches). A LOW SCORE IS NOT EVIDENCE OF INACTIVITY.",
    "# Substrate: 4-propylguaiacol conversion only. Says nothing about ring cleavage.",
    f"# calib_distance: fraction of differing residues at positions {DIST_POS} to the nearest of",
    f"#   the {len(cal)} characterized enzymes. These positions are DISJOINT from the scored five.",
    f"# calib_domain: 'inside' at distance <= {MAXD:.2f}, the largest nearest-neighbour distance",
    "#   WITHIN the calibration set. This is a COVERAGE statement, not a confidence: distance was",
    "#   tested against model error and did NOT predict it (p = 0.88). See MODEL_VALIDATION.md F.",
    "# NO CLADE CALL: this tool does not run phylogenetic placement, so the score's central limit",
    "#   (that it detects clade-V-LIKE pockets) cannot be checked per protein here.",
]
buf = "\n".join(head + ["\t".join(r) for r in out_rows]) + "\n"
if a.out.endswith(".zst"):
    p = subprocess.Popen(["zstd", "-qf", "-o", a.out], stdin=subprocess.PIPE, text=True)
    p.communicate(buf)
else:
    open(a.out, "w").write(buf)

log(f"{n} PF01328 proteins; {flagged} score >= {a.cut}/5; "
    f"{outside} ({100*outside/max(n,1):.1f} %) outside the calibration domain")
log(f"wrote {a.out}")
if tmp and not a.keep_sites and os.path.exists(tmp):
    os.remove(tmp)
