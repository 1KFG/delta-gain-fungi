#!/usr/bin/env python3.12
"""
Applicability domain for the T1 pocket score: WHICH proteins is the score entitled to speak about?

THE PROBLEM THIS ADDRESSES
--------------------------
The T1 score currently ships with ONE global caveat: "9 of 10 active calibration enzymes are
clade V, out-of-clade AUC is 0.556, a low score is not evidence of inactivity." That caveat is
true, but it is the same sentence for all 97,520 proteins. It cannot distinguish a protein whose
active site closely resembles a characterized enzyme from one that resembles nothing ever
measured. Those two cases deserve different levels of trust and currently get the same one.

A model is not made reliable by a better point estimate. It is made reliable by knowing where it
applies. That is a per-protein property and it is measurable without any new biochemistry.

THE DISTANCE, AND WHY IT IS INDEPENDENT OF THE SCORE
-----------------------------------------------------
Distance is the fraction of differing residues at active-site positions that the score does NOT
use. The score uses POCKET5 = (73, 77, 81, 183, 187). This distance uses the other active-site
positions shared between the calibration table (Figure S2) and the BFD site layer:

    80, 84, 108, 110, 172, 173, 175, 176, 179, 180

No position is used by both. A protein can therefore sit far from the calibration set while
scoring 5/5, and that combination is exactly the case the global caveat hides.

THE TEST (this is what decides whether the distance is worth emitting)
-----------------------------------------------------------------------
Within the 19 calibration enzymes, leave-one-out:
  1. distance from the held-out enzyme to its NEAREST neighbour among the other 18;
  2. the pocket model's absolute LOO prediction error for that enzyme (from score_t1.py).
Spearman correlation between the two, with a permutation test.

If distance predicts error, then distance is a usable reliability statement and BFD proteins far
from all 19 can be flagged as outside the calibration domain. If it does not predict error, that
is reported as a null result and the flag is NOT dressed up as a confidence score. Either way the
fraction of BFD lying outside the range of distances ever seen among the 19 is a fact worth
knowing, because it bounds how far the score has been pushed.

LIMITS
------
n = 19. A Spearman on 19 points is weak and the permutation p is reported so the reader can see
how often noise does this well. 10 positions give a coarse distance: possible values are
multiples of 0.1. Ties are frequent and are handled by the tie-corrected rank statistic.
Only 4 of the 19 enzymes have a public full-length sequence, so a whole-sequence identity
version of this test cannot be run. This is a limit of the published data, not a choice.
"""
from __future__ import annotations
import argparse, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--s2", required=True, help="active_site_S2.tsv")
ap.add_argument("--loo", required=True, help="t1_score_loo.tsv from score_t1.py")
ap.add_argument("--bfd-all", default=None, help="upo_sites_bfd_all.tsv.zst (site layer)")
ap.add_argument("--bfd-scored", default=None, help="upo_sites_bfd_scored.tsv.zst")
ap.add_argument("--out", default=None, help="per-protein domain table (.tsv.zst)")
ap.add_argument("--permutations", type=int, default=10000)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--report", default=None)
a = ap.parse_args()

POCKET5 = (73, 77, 81, 183, 187)
# active-site positions in BOTH the S2 calibration table and the BFD site layer, minus POCKET5
DIST_POS = (80, 84, 108, 110, 172, 173, 175, 176, 179, 180)

OUT = []
def say(s=""):
    print(s, file=sys.stderr)
    OUT.append(s)


def read_tsv(p):
    rows, hdr = [], None
    for line in open(p):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))
    return rows


# ---------------- calibration profiles ----------------
s2 = [r for r in read_tsv(a.s2) if r["enzyme"] != "LOCATION"]
loo = {r["enzyme"]: r for r in read_tsv(a.loo)}

enz, prof, clade = [], [], []
for r in s2:
    if r["enzyme"] not in loo:
        continue
    aa = [r.get(f"p{p}", "").strip().upper() for p in DIST_POS]
    enz.append(r["enzyme"]); prof.append(aa); clade.append(r["clade"])
n = len(enz)
say(f"# Applicability domain of the T1 pocket score\n")
say(f"{n} calibration enzymes with both an active-site profile and a LOO prediction.")
say(f"Distance positions (disjoint from the scored POCKET5): {', '.join(str(p) for p in DIST_POS)}\n")


def dist(p, q):
    """Fraction of differing residues over positions where BOTH are observed."""
    num = den = 0
    for x, y in zip(p, q):
        if not x or x in "-." or not y or y in "-.":
            continue
        den += 1
        num += (x != y)
    return (num / den) if den else float("nan")


# nearest-neighbour distance, leave-one-out
nnd = np.array([min(dist(prof[i], prof[j]) for j in range(n) if j != i) for i in range(n)])
err = np.array([abs(float(loo[e]["activity_mean"]) - float(loo[e]["loo_pred_pocket"])) for e in enz])
act = np.array([int(loo[e]["active"]) for e in enz])


def spearman(u, v):
    def rk(x):
        o = np.argsort(x, kind="mergesort")
        r = np.empty(len(x), float)
        sx = np.asarray(x, float)[o]
        i = 0
        while i < len(sx):
            j = i
            while j + 1 < len(sx) and sx[j + 1] == sx[i]:
                j += 1
            r[o[i:j + 1]] = 0.5 * (i + j) + 1.0
            i = j + 1
        return r
    ru, rv = rk(u), rk(v)
    ru -= ru.mean(); rv -= rv.mean()
    d = np.sqrt((ru ** 2).sum() * (rv ** 2).sum())
    return float((ru * rv).sum() / d) if d else float("nan")


rho = spearman(nnd, err)
rng = np.random.default_rng(a.seed)
cnt = sum(1 for _ in range(a.permutations)
          if spearman(nnd, rng.permutation(err)) >= rho)
p_emp = (cnt + 1) / (a.permutations + 1)

say("## Does distance from the calibration set predict the model's error?\n")
say(f"Nearest-neighbour distance among the {n}: "
    f"median {np.median(nnd):.2f}, range {nnd.min():.2f}-{nnd.max():.2f}.\n")

# --- the confound, which has to be dealt with before the raw correlation means anything ---
# Distance is almost perfectly confounded with clade. Every enzyme at distance <= 0.10 is clade V,
# and clade V holds 9 of the 10 actives, whose activities are far from the overall mean. ABSOLUTE
# error is larger wherever the response is far from the mean, for ANY model including the null.
# So a raw Spearman(distance, |error|) measures the activity distribution, not the model.
# The question that is actually about the model is whether its SKILL over the baseline decays
# with distance: skill = |null error| - |pocket error|, positive when the pocket model helps.
err_null = np.array([abs(float(loo[e]["activity_mean"]) - float(loo[e]["loo_pred_null"]))
                     for e in enz])
skill = err_null - err
v_lo = [i for i in range(n) if nnd[i] <= 0.10]
v_hi = [i for i in range(n) if nnd[i] >= 0.30]
say(f"Distance is confounded with clade: {sum(1 for i in v_lo if clade[i]=='V')} of "
    f"{len(v_lo)} enzymes at distance <= 0.10 are clade V, against "
    f"{sum(1 for i in v_hi if clade[i]=='V')} of {len(v_hi)} at >= 0.30. Clade V holds "
    f"{int(act[[i for i in range(n) if clade[i]=='V']].sum())} of the {int(act.sum())} actives, "
    f"whose activities sit far from the mean. Absolute error is larger there for ANY model, "
    f"including the null. The raw correlation below therefore measures the activity "
    f"distribution, not the model.\n")
say(f"Raw Spearman(distance, |LOO error|) = {rho:+.3f}, one-sided permutation p = {p_emp:.4f} "
    f"({a.permutations} shuffles). **Confounded - do not read this as a result.**\n")

rho_s = spearman(nnd, -skill)
cnt_s = sum(1 for _ in range(a.permutations)
            if spearman(nnd, rng.permutation(-skill)) >= rho_s)
p_skill = (cnt_s + 1) / (a.permutations + 1)
say("### The unconfounded version: does the model's SKILL decay with distance?\n")
say("Skill = |null-model error| - |pocket-model error|. Positive means the pocket model beat the\n"
    "baseline on that enzyme. Both models face the same response, so the magnitude of the\n"
    "response cancels.\n")
say(f"Mean skill, distance <= 0.10: {np.mean(skill[v_lo]):+.1f} "
    f"(n = {len(v_lo)}); distance >= 0.30: {np.mean(skill[v_hi]):+.1f} (n = {len(v_hi)}).\n")
say(f"Spearman(distance, -skill) = {rho_s:+.3f}, one-sided permutation p = {p_skill:.4f}.\n")
rho, p_emp = rho_s, p_skill
if p_emp <= 0.05 and rho > 0:
    say("**Distance predicts error.** Proteins far from every characterized enzyme are measurably\n"
        "harder for the score, so distance is a usable per-protein reliability statement.\n")
    verdict = "predictive"
else:
    say("**Distance does NOT predict error on this calibration set.** This is a null result and is\n"
        "reported as one. The domain flag below is therefore a statement about COVERAGE - whether a\n"
        "protein's active site resembles anything ever measured - and NOT a calibrated confidence.\n"
        "It must not be presented as a probability of being right.\n")
    verdict = "coverage_only"

say("| enzyme | clade | active | NN distance | \\|LOO error\\| |")
say("|---|---|---|---|---|")
for i in np.argsort(-nnd):
    say(f"| {enz[i]} | {clade[i]} | {act[i]} | {nnd[i]:.2f} | {err[i]:.1f} |")
say("")

maxd = float(np.nanmax(nnd))
say(f"The largest nearest-neighbour distance inside the calibration set is {maxd:.2f}. "
    f"A BFD protein further than that from all {n} enzymes sits outside the range of active-site "
    f"divergence the score was ever calibrated over.\n")

if not (a.bfd_all and a.bfd_scored):
    if a.report:
        open(a.report, "w").write("\n".join(OUT) + "\n")
    sys.exit(0)

# ---------------- apply to BFD ----------------
say("## How much of BFD lies outside that range?\n")
cols = [f"p{p}" for p in DIST_POS]
proc = subprocess.Popen(["zstd", "-dcq", a.bfd_all], stdout=subprocess.PIPE, text=True)
hdr, ix = None, None
res = {}
for line in proc.stdout:
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        ix = {}
        for c in cols:
            hit = [j for j, k in enumerate(hdr) if k.split("_")[0] == c]
            if hit:
                ix[c] = hit[0]
        missing = [c for c in cols if c not in ix]
        if missing:
            sys.exit(f"BFD site layer lacks columns {missing}")
        pid_i = hdr.index("protein_id")
        continue
    aa = [f[ix[c]].strip().upper() if ix[c] < len(f) else "" for c in cols]
    res[f[pid_i]] = aa
proc.wait()
say(f"{len(res)} proteins read from the site layer.\n")

P = [p for p in prof]
def nn_to_calib(aa):
    ds = [dist(aa, q) for q in P]
    ds = [d for d in ds if d == d]
    return min(ds) if ds else float("nan")

cache = {}
rows = []
proc = subprocess.Popen(["zstd", "-dcq", a.bfd_scored], stdout=subprocess.PIPE, text=True)
hdr2 = None
nout = 0
inside = outside = unknown = 0
out_by_clade = {}
out_hi = 0
hi_total = 0
for line in proc.stdout:
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if hdr2 is None:
        hdr2 = f
        j = {k: i for i, k in enumerate(hdr2)}
        rows.append(hdr2 + ["calib_distance", "calib_domain"])
        continue
    pid = f[j["protein_id"]]
    aa = res.get(pid)
    if aa is None:
        d = float("nan")
    else:
        key = tuple(aa)
        if key not in cache:
            cache[key] = nn_to_calib(aa)
        d = cache[key]
    if d != d:
        dom = "unknown"; unknown += 1
    elif d <= maxd:
        dom = "inside"; inside += 1
    else:
        dom = "outside"; outside += 1
        out_by_clade[f[j["clade"]]] = out_by_clade.get(f[j["clade"]], 0) + 1
    try:
        m = int(f[j["t1_pocket_matches"]])
    except ValueError:
        m = -1
    if m >= 3:
        hi_total += 1
        if dom == "outside":
            out_hi += 1
    rows.append(f + ["" if d != d else f"{d:.3f}", dom])
    nout += 1
proc.wait()

tot = inside + outside + unknown
say(f"| domain | n | % |")
say("|---|---|---|")
for lab, v in (("inside", inside), ("outside", outside), ("unknown (no profile)", unknown)):
    say(f"| {lab} | {v} | {100*v/tot:.1f} |")
say("")
say(f"**{out_hi} of the {hi_total} proteins scoring >= 3/5 ({100*out_hi/max(hi_total,1):.1f} %) sit "
    f"OUTSIDE the calibration range.** For those, a high pocket score is an extrapolation on both "
    f"axes at once: the clade is usually not V, and the rest of the active site resembles nothing "
    f"that has been measured.\n")
if out_by_clade:
    say("Outside-domain proteins by clade call:\n")
    say("| clade | n |")
    say("|---|---|")
    for k, v in sorted(out_by_clade.items(), key=lambda x: -x[1]):
        say(f"| {k} | {v} |")
    say("")

say("### How to read the two columns this adds\n")
say("`calib_distance` is the fraction of differing residues, over the 10 non-scored active-site\n"
    "positions, to the nearest of the 19 characterized enzymes.\n")
say(f"`calib_domain` is `inside` when that distance is within the largest distance seen WITHIN the\n"
    f"calibration set ({maxd:.2f}), `outside` otherwise.\n")
if verdict == "coverage_only":
    say("`outside` does NOT mean a protein is inactive, and `inside` does NOT mean the score is\n"
        "right about it. Distance did not predict model error on 19 enzymes. The column says only\n"
        "whether the protein's active site resembles something that has been measured. Use it to\n"
        "decide what to assay, not to decide what is active.\n")

if a.out:
    buf = ["# calib_distance / calib_domain: applicability domain of the T1 pocket score.",
           f"# Distance = fraction of differing residues at positions {DIST_POS},",
           "#   which are DISJOINT from the scored POCKET5 (73,77,81,183,187).",
           f"# 'inside' = within {maxd:.2f}, the largest nearest-neighbour distance among the 19.",
           f"# Spearman(distance, |LOO error|) = {rho:+.3f}, permutation p = {p_emp:.4f}.",
           "# Distance is a COVERAGE statement, not a calibrated confidence." if
           verdict == "coverage_only" else
           "# Distance predicted model error on the calibration set.",
           ] + ["\t".join(r) for r in rows]
    pz = subprocess.Popen(["zstd", "-qf", "-o", a.out], stdin=subprocess.PIPE, text=True)
    pz.communicate("\n".join(buf) + "\n")
    print(f"wrote {nout} rows to {a.out}", file=sys.stderr)

if a.report:
    open(a.report, "w").write("\n".join(OUT) + "\n")
    print(f"wrote {a.report}", file=sys.stderr)
