#!/usr/bin/env python3.12
"""
Build and HONESTLY cross-validate a T1 score: can active-site residues predict peroxygenase
activity on 4-propylguaiacol?

WHY THIS EXISTS
---------------
The final table deliberately emitted no activity score. That was a judgement made ahead of the
measurement. The user asked for a score. The correct answer is to build one, cross-validate it
OUT OF SAMPLE, and report what it is worth - not to refuse in advance and not to ship an
in-sample correlation as if it were predictive.

WHAT IS MEASURED
----------------
Response: mean peroxygenase conversion fraction over the three pH values in Figure S8 (n = 19).

Four models, each scored by LEAVE-ONE-OUT cross-validation. Everything a model learns - the
active-residue consensus, the clade means, the regression coefficients - is refit on the 18
training enzymes only. The held-out enzyme never contributes to anything used to predict it.

  null     predict the training mean for everyone
  clade    predict the held-out enzyme's clade mean, from training enzymes only
  pocket   regress on the 5 alignment-reliable aliphatic pocket positions
  both     clade mean plus pocket

THE QUESTION THAT DECIDES IT
----------------------------
Does `pocket` beat `null`, and does `both` beat `clade`? If the pocket adds nothing over knowing
the clade, then a pocket-based score is a clade label wearing a number, and it must not be sold
as structural insight. In-sample r is not evidence here: with n = 19 and 5 free positions,
in-sample fit is close to guaranteed.

LIMITS, read before using any number from this
----------------------------------------------
- n = 19. A leave-one-out estimate on 19 points has a wide confidence interval. It is reported
  with a permutation test so the reader can see how often noise does this well.
- One substrate, one assay, one laboratory, one publication.
- Figure S8 values are printed data labels transcribed by eye, not author-supplied numbers.
- Any score this produces applies to 4-propylguaiacol peroxygenase conversion. It is NOT a
  general activity score and says nothing about which rings an enzyme cleaves.
"""
from __future__ import annotations
import argparse, itertools, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--s2", required=True)
ap.add_argument("--t1", required=True)
ap.add_argument("--pocket", default="73,77,81,183,187")
ap.add_argument("--active-threshold", type=float, default=10.0)
ap.add_argument("--permutations", type=int, default=10000)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--out", default=None)
a = ap.parse_args()


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


s2 = [r for r in read_tsv(a.s2) if r["enzyme"] != "LOCATION"]
t1 = read_tsv(a.t1)
pocket_cols = [f"p{p}" for p in a.pocket.split(",")]

act = {}
for r in t1:
    v = [float(r[k]) for k in ("perox_pH3.5", "perox_pH5.5", "perox_pH8.0")]
    act[r["enzyme"]] = float(np.mean(v))

enz, clade, res, y = [], [], [], []
for r in s2:
    if r["enzyme"] not in act:
        continue
    aa = [r[c].strip() for c in pocket_cols]
    if any(x in ("", "-") for x in aa):
        print(f"  skip {r['enzyme']}: missing pocket residue", file=sys.stderr)
        continue
    enz.append(r["enzyme"]); clade.append(r["clade"]); res.append(aa); y.append(act[r["enzyme"]])

y = np.array(y); n = len(y)
print(f"{n} enzymes with both pocket residues and T1 activity", file=sys.stderr)
print(f"activity: mean {y.mean():.1f}, median {np.median(y):.1f}, "
      f"range {y.min():.1f}-{y.max():.1f}", file=sys.stderr)
binary = (y > a.active_threshold).astype(int)
print(f"active (>{a.active_threshold}%): {binary.sum()} of {n}", file=sys.stderr)


def consensus_feature(train_idx, test_idx):
    """Residue-match count against the ACTIVE consensus learned from training enzymes only."""
    # Uses the CURRENT global y, so that under permutation the active set permutes too.
    # If it used the fixed `binary`, the null model would keep real information and the
    # permutation p-value would be anticonservative.
    tr_act = [i for i in train_idx if y[i] > a.active_threshold]
    if not tr_act:
        return np.zeros(len(test_idx))
    cons = []
    for k in range(len(pocket_cols)):
        col = [res[i][k] for i in tr_act]
        cons.append(max(set(col), key=col.count))
    return np.array([sum(res[i][k] == cons[k] for k in range(len(cons))) for i in test_idx],
                    dtype=float)


def fit_predict(train_idx, test_idx, model):
    tr_y = y[train_idx]
    if model == "null":
        return np.full(len(test_idx), tr_y.mean())
    if model == "clade":
        gm = {}
        for c in set(clade):
            v = [y[i] for i in train_idx if clade[i] == c]
            gm[c] = float(np.mean(v)) if v else float(tr_y.mean())
        return np.array([gm.get(clade[i], tr_y.mean()) for i in test_idx])
    xf_tr = consensus_feature(train_idx, train_idx)
    xf_te = consensus_feature(train_idx, test_idx)
    if model == "pocket":
        X_tr = np.c_[np.ones(len(train_idx)), xf_tr]
        X_te = np.c_[np.ones(len(test_idx)), xf_te]
    else:  # both
        gm = {}
        for c in set(clade):
            v = [y[i] for i in train_idx if clade[i] == c]
            gm[c] = float(np.mean(v)) if v else float(tr_y.mean())
        g_tr = np.array([gm[clade[i]] for i in train_idx])
        g_te = np.array([gm.get(clade[i], tr_y.mean()) for i in test_idx])
        X_tr = np.c_[np.ones(len(train_idx)), xf_tr, g_tr]
        X_te = np.c_[np.ones(len(test_idx)), xf_te, g_te]
    beta, *_ = np.linalg.lstsq(X_tr, tr_y, rcond=None)
    return X_te @ beta


def loo(model, yv=None):
    global y
    keep = y
    y = keep if yv is None else yv
    try:
        pred = np.empty(n)
        for i in range(n):
            tr = [j for j in range(n) if j != i]
            pred[i] = fit_predict(tr, [i], model)[0]
    finally:
        y = keep
    return pred


def r2(yt, yp):
    ss = ((yt - yp) ** 2).sum()
    return 1.0 - ss / ((yt - yt.mean()) ** 2).sum()


def spearman(u, v):
    ru = np.argsort(np.argsort(u)).astype(float)
    rv = np.argsort(np.argsort(v)).astype(float)
    ru -= ru.mean(); rv -= rv.mean()
    d = np.sqrt((ru ** 2).sum() * (rv ** 2).sum())
    return float((ru * rv).sum() / d) if d else float("nan")


def auc(lab, sc):
    pos = sc[lab == 1]; neg = sc[lab == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    w = sum((1.0 if p > q else 0.5 if p == q else 0.0) for p in pos for q in neg)
    return w / (len(pos) * len(neg))


models = ["null", "clade", "pocket", "both"]
results = {}
print("\n=== LEAVE-ONE-OUT, out of sample ===", file=sys.stderr)
print(f"{'model':8} {'LOO R2':>9} {'LOO rho':>9} {'LOO AUC':>9} {'MAE':>7}", file=sys.stderr)
for m in models:
    p = loo(m)
    results[m] = dict(r2=r2(y, p), rho=spearman(y, p), auc=auc(binary, p),
                      mae=float(np.abs(y - p).mean()), pred=p.tolist())
    r = results[m]
    print(f"{m:8} {r['r2']:9.3f} {r['rho']:9.3f} {r['auc']:9.3f} {r['mae']:7.1f}", file=sys.stderr)

# in-sample, to show the gap
print("\n=== IN-SAMPLE (shown only to expose the overfit gap) ===", file=sys.stderr)
for m in models:
    p = fit_predict(list(range(n)), list(range(n)), m)
    print(f"{m:8} in-sample R2 {r2(y, p):6.3f}   LOO R2 {results[m]['r2']:6.3f}", file=sys.stderr)

# permutation test: how often does a shuffled response do this well?
rng = np.random.default_rng(a.seed)
print(f"\n=== PERMUTATION TEST ({a.permutations} shuffles of the response) ===", file=sys.stderr)
perm = {}
for m in ("pocket", "both"):
    obs = results[m]["r2"]
    cnt = 0
    for _ in range(a.permutations):
        yp = rng.permutation(y)
        if r2(yp, loo(m, yp)) >= obs:
            cnt += 1
    p_emp = (cnt + 1) / (a.permutations + 1)
    perm[m] = p_emp
    print(f"{m:8} LOO R2 = {obs:+.3f}   p = {p_emp:.4f}", file=sys.stderr)

print("\n=== VERDICT ===", file=sys.stderr)
gain_null = results["pocket"]["r2"] - results["null"]["r2"]
gain_clade = results["both"]["r2"] - results["clade"]["r2"]
print(f"pocket over null:  LOO R2 {gain_null:+.3f}", file=sys.stderr)
print(f"both over clade:   LOO R2 {gain_clade:+.3f}", file=sys.stderr)

if a.out:
    with open(a.out, "w") as fh:
        fh.write("enzyme\tclade\tactivity_mean\tactive\t" +
                 "\t".join(f"loo_pred_{m}" for m in models) + "\n")
        for i in range(n):
            fh.write(f"{enz[i]}\t{clade[i]}\t{y[i]:.1f}\t{binary[i]}\t" +
                     "\t".join(f"{results[m]['pred'][i]:.2f}" for m in models) + "\n")
    print(f"\nwrote {a.out}", file=sys.stderr)

# ---------------------------------------------------------------------------
# THE DECISIVE STRATIFICATION
# ---------------------------------------------------------------------------
# 9 of the 10 active enzymes are clade V. A model that only learned "clade V is active" would
# still post a high overall AUC. The question that matters for applying a score to BFD is whether
# the pocket can find an active enzyme OUTSIDE the active clade, because that is the only case
# where a residue-based score adds anything to a clade label.
#
# So: recompute AUC with all clade V enzymes removed. What is left is 1 active (GmaUPO-II) against
# the inactive non-V enzymes. The resulting AUC is the honest measure of out-of-clade detection,
# and its effective sample size is the number of out-of-clade actives, which is ONE.
print("\n=== AUC WITH CLADE V REMOVED (out-of-clade detection) ===", file=sys.stderr)
nv = [i for i in range(n) if clade[i] != "V"]
nv_lab = binary[nv]
print(f"non-V enzymes: {len(nv)}, of which active: {nv_lab.sum()} "
      f"({', '.join(enz[i] for i in nv if binary[i] == 1)})", file=sys.stderr)
for m in models:
    p = np.array(results[m]["pred"])[nv]
    print(f"{m:8} AUC(non-V) {auc(nv_lab, p):6.3f}   vs AUC(all) {results[m]['auc']:6.3f}",
          file=sys.stderr)
print("\nEffective sample size for the out-of-clade question is the number of actives outside",
      file=sys.stderr)
print(f"clade V: {int(nv_lab.sum())}. No AUC computed on {int(nv_lab.sum())} positive case is",
      file=sys.stderr)
print("evidence of anything. It is reported to show the limit, not to claim a result.",
      file=sys.stderr)
