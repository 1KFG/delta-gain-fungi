#!/usr/bin/env python3.12
"""
Can STRUCTURE predict T1 activity where residue identity cannot?

WHY THIS IS THE ONE ROUTE LEFT
------------------------------
The residue-based pocket score works (LOO AUC 0.956) but is confounded with clade: 9 of its 10
active calibration enzymes are clade V, and with clade V removed its AUC is 0.556 on one active.
A residue consensus learned almost entirely inside one clade is a clade marker by another name.

Geometry is a different kind of feature. A heme access channel wide enough to admit
4-propylguaiacol is a physical requirement that does not care which clade an enzyme belongs to.
If channel geometry predicts activity, and keeps predicting it with clade V removed, then the
model would generalize out of clade for a mechanical reason rather than a phylogenetic one. That
is the specific gap the project has been unable to close.

The 19 Boltz-2 structures exist for exactly this question (`CHANNEL_FINDINGS.md` records that the
channel work was stuck at n = 4 until they were built). This script is the first use of them
against T1.

WHAT IS MEASURED
----------------
Features, from `heme_channel.py` over each Boltz model:
  bottleneck  narrowest diameter along the corridor from the heme iron
  d1..d9      widest clearance at 1, 3, 5, 7 and 9 A from the iron (the corridor profile)

Models, all scored by LEAVE-ONE-OUT with every coefficient refit on the 18 training enzymes:
  null        training mean
  clade       training clade mean
  pocket      the residue score (the incumbent)
  channel     bottleneck alone
  profile     the corridor profile
  pocket+chan pocket plus bottleneck

THE TEST THAT DECIDES IT
------------------------
AUC with clade V REMOVED. Everything else is secondary. If a structural feature cannot beat 0.556
out of clade, it has not solved the problem that matters, however well it does overall.

LIMITS, read before using any number here
-----------------------------------------
- n = 19, of which 10 active and only 1 active outside clade V. The out-of-clade AUC still rests
  on ONE positive case. A better number there is not proof; it is a weaker form of the same
  problem, and is reported as such.
- The structures are PREDICTED (Boltz-2), not experimental. Confidence is high (complex pLDDT
  about 0.91-0.93) but a predicted side-chain rotamer can move a bottleneck by an angstrom.
- Channel geometry is measured on an apo model with no substrate bound. Real channels flex.
- Same substrate caveat as everywhere else: 4-propylguaiacol only.
"""
from __future__ import annotations
import argparse, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--channel", required=True, help="boltz19_channel_profile.tsv")
ap.add_argument("--s2", required=True)
ap.add_argument("--t1", required=True)
ap.add_argument("--pocket", default="73,77,81,183,187")
ap.add_argument("--active-threshold", type=float, default=10.0)
ap.add_argument("--permutations", type=int, default=10000)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--report", default=None)
ap.add_argument("--out", default=None)
a = ap.parse_args()

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


chan = {}
for r in read_tsv(a.channel):
    # the Boltz run labelled AaeUPO's model by its evolved variant's parentage
    nm = r["structure"].replace("_parent_AaeUPO", "")
    chan[nm] = r

s2 = [r for r in read_tsv(a.s2) if r["enzyme"] != "LOCATION"]
t1 = read_tsv(a.t1)
pocket_cols = [f"p{p}" for p in a.pocket.split(",")]
act = {r["enzyme"]: float(np.mean([float(r[k]) for k in
       ("perox_pH3.5", "perox_pH5.5", "perox_pH8.0")])) for r in t1}

FEATS = ["bottleneck", "d1", "d3", "d5", "d7", "d9"]
enz, clade, res, y, X = [], [], [], [], []
for r in s2:
    e = r["enzyme"]
    if e not in act or e not in chan:
        continue
    aa = [r[c].strip() for c in pocket_cols]
    if any(x in ("", "-") for x in aa):
        continue
    try:
        xv = [float(chan[e][k]) for k in FEATS]
    except (KeyError, ValueError):
        continue
    enz.append(e); clade.append(r["clade"]); res.append(aa); y.append(act[e]); X.append(xv)

y = np.array(y); X = np.array(X); n = len(y)
binary = (y > a.active_threshold).astype(int)
say("# Structure versus residues for T1 activity\n")
say(f"{n} enzymes with a Boltz structure, pocket residues and a T1 value. "
    f"Active (> {a.active_threshold:.0f} %): {binary.sum()}.\n")
missing = sorted({r['enzyme'] for r in s2 if r['enzyme'] in act} - set(enz))
if missing:
    say(f"Excluded for want of a structure or residue: {', '.join(missing)}\n")

FI = {k: i for i, k in enumerate(FEATS)}


def consensus_feature(train_idx, test_idx):
    tr_act = [i for i in train_idx if y[i] > a.active_threshold]
    if not tr_act:
        return np.zeros(len(test_idx))
    cons = []
    for k in range(len(pocket_cols)):
        col = [res[i][k] for i in tr_act]
        cons.append(max(set(col), key=col.count))
    return np.array([sum(res[i][k] == cons[k] for k in range(len(cons))) for i in test_idx],
                    dtype=float)


def design(idx, model, train_idx):
    one = np.ones(len(idx))
    if model == "channel":
        return np.c_[one, X[idx][:, FI["bottleneck"]]]
    if model == "profile":
        return np.c_[one, X[idx][:, [FI["d1"], FI["d3"], FI["d5"]]]]
    if model == "pocket":
        return np.c_[one, consensus_feature(train_idx, idx)]
    if model == "pocket+chan":
        return np.c_[one, consensus_feature(train_idx, idx), X[idx][:, FI["bottleneck"]]]
    raise ValueError(model)


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
    Xtr = design(train_idx, model, train_idx)
    Xte = design(test_idx, model, train_idx)
    beta, *_ = np.linalg.lstsq(Xtr, tr_y, rcond=None)
    return Xte @ beta


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
    return 1.0 - ((yt - yp) ** 2).sum() / ((yt - yt.mean()) ** 2).sum()


def auc(lab, sc):
    pos, neg = sc[lab == 1], sc[lab == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    return sum(1.0 if p > q else 0.5 if p == q else 0.0
               for p in pos for q in neg) / (len(pos) * len(neg))


models = ["null", "clade", "pocket", "channel", "profile", "pocket+chan"]
preds = {m: loo(m) for m in models}

say("## Leave-one-out, out of sample\n")
say("| model | LOO R2 | LOO AUC (all) | **LOO AUC, clade V removed** |")
say("|---|---|---|---|")
nv = [i for i in range(n) if clade[i] != "V"]
nv_lab = binary[nv]
for m in models:
    p = preds[m]
    say(f"| {m} | {r2(y, p):+.3f} | {auc(binary, p):.3f} | **{auc(nv_lab, p[nv]):.3f}** |")
say("")
say(f"Non-clade-V enzymes: {len(nv)}, of which active: {int(nv_lab.sum())} "
    f"({', '.join(enz[i] for i in nv if binary[i] == 1)}).\n")

rng = np.random.default_rng(a.seed)
say("## Permutation test on the structural models\n")
say("| model | LOO R2 | p |")
say("|---|---|---|")
for m in ("channel", "profile", "pocket+chan"):
    obs = r2(y, preds[m])
    cnt = 0
    for _ in range(a.permutations):
        yp = rng.permutation(y)
        if r2(yp, loo(m, yp)) >= obs:
            cnt += 1
    say(f"| {m} | {obs:+.3f} | {(cnt + 1) / (a.permutations + 1):.4f} |")
say("")

say("## Per-enzyme channel geometry against activity\n")
say("| enzyme | clade | active | activity | bottleneck | d3 | d5 |")
say("|---|---|---|---|---|---|---|")
for i in np.argsort(-y):
    say(f"| {enz[i]} | {clade[i]} | {binary[i]} | {y[i]:.1f} | "
        f"{X[i][FI['bottleneck']]:.2f} | {X[i][FI['d3']]:.1f} | {X[i][FI['d5']]:.1f} |")
say("")

b = X[:, FI["bottleneck"]]
say(f"Bottleneck, active enzymes: mean {b[binary == 1].mean():.2f} A "
    f"(range {b[binary == 1].min():.2f}-{b[binary == 1].max():.2f}). "
    f"Inactive: mean {b[binary == 0].mean():.2f} A "
    f"(range {b[binary == 0].min():.2f}-{b[binary == 0].max():.2f}).\n")

best_nv = max(models, key=lambda m: (auc(nv_lab, preds[m][nv]) if
                                     auc(nv_lab, preds[m][nv]) == auc(nv_lab, preds[m][nv])
                                     else -1))
say("## Verdict\n")
say(f"Out-of-clade AUC by model: " +
    ", ".join(f"{m} {auc(nv_lab, preds[m][nv]):.3f}" for m in models) + ".\n")
say(f"The out-of-clade question rests on {int(nv_lab.sum())} active enzyme outside clade V. "
    f"No AUC computed on {int(nv_lab.sum())} positive case is evidence, whichever model wins it. "
    f"A structural feature that beat the residue score here would be a HYPOTHESIS worth testing "
    f"on new enzymes, not a validated out-of-clade model.\n")

if a.out:
    with open(a.out, "w") as fh:
        fh.write("enzyme\tclade\tactivity\tactive\tbottleneck\t" +
                 "\t".join(f"loo_{m}" for m in models) + "\n")
        for i in range(n):
            fh.write(f"{enz[i]}\t{clade[i]}\t{y[i]:.1f}\t{binary[i]}\t"
                     f"{X[i][FI['bottleneck']]:.2f}\t" +
                     "\t".join(f"{preds[m][i]:.2f}" for m in models) + "\n")
    print(f"wrote {a.out}", file=sys.stderr)
if a.report:
    open(a.report, "w").write("\n".join(OUT) + "\n")
    print(f"wrote {a.report}", file=sys.stderr)
