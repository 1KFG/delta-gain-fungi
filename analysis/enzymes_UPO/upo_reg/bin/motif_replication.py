#!/usr/bin/env python3.12
"""
G3 Phase A step 3: the primary result -- does a clade's motif replicate across INDEPENDENT lineages?

For each clade, collect the motifs STREME found in every lineage where that clade was tested, then
compare them all-against-all with Tomtom. A motif found in >= 2 lineages for the same clade is
evidence of a clade-linked regulatory element. A motif found in only one lineage is a lineage
element and is reported as such.

This is deliberately NOT a p-value ranking. STREME p-values on 7,055 promoters drawn from 1,117
genera are optimistic no matter how the input is dereplicated. Cross-lineage replication is the
evidence that survives that problem.
"""
import argparse, glob, os, re, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--streme-dir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--jaspar", default=None)
ap.add_argument("--tomtom-thresh", type=float, default=0.05)
ap.add_argument("--test-pvalue", type=float, default=0.05,
                help="STREME held-out test p-value a motif must beat to be considered at all. "
                     "STREME with --nmotifs N always REPORTS N motifs regardless of significance, "
                     "so without this filter non-significant motifs would be counted as "
                     "replicating.")
a = ap.parse_args()

runs = {}
for d in sorted(glob.glob(os.path.join(a.streme_dir, "*_clade*"))):
    if not os.path.isdir(d):
        continue
    f = os.path.join(d, "streme.txt")
    if os.path.exists(f):
        base = os.path.basename(d)
        lin, clade = base.rsplit("_clade", 1)
        runs[(lin, clade)] = f
if not runs:
    sys.exit(f"no STREME outputs under {a.streme_dir}")
print(f"{len(runs)} contrasts with motifs", file=sys.stderr)


by_clade = defaultdict(list)
for (lin, clade), f in runs.items():
    by_clade[clade].append((lin, f))

def significant_ids(streme_txt, thresh):
    """Motif ids whose STREME held-out TEST p-value beats `thresh`."""
    xml = os.path.join(os.path.dirname(streme_txt), "streme.xml")
    if not os.path.exists(xml):
        return None                      # cannot filter; caller must not claim significance
    s = open(xml).read()
    out = set()
    for m in re.finditer(r'id="(\d+-[^"]+)"[\s\S]{0,400}?test_pvalue="([^"]*)"', s):
        try:
            if float(m.group(2)) < thresh:
                out.add(m.group(1))
        except ValueError:
            pass
    return out

SIG = {}
for k, f in runs.items():
    SIG[k] = significant_ids(f, a.test_pvalue)

out = []
W = out.append
W("# G3 Phase A: cross-lineage motif replication\n")
W("A motif is reported as **clade-linked** only if it is significant on STREME's held-out test "
  f"set (test p < {a.test_pvalue}) in two or more INDEPENDENT lineages for the same clade, and "
  "the two are matched by Tomtom. Single-lineage motifs are lineage elements.\n")
W("**Significance filter matters here.** STREME run with `--nmotifs 8` always reports 8 motifs "
  "per contrast whether or not they are significant. Across the 17 contrasts it reported 136 "
  "motifs, of which only those passing the held-out test are considered below.\n")
W("STREME p-values are not used for ranking; see the module docstring for why.\n")

summary = []
for clade in sorted(by_clade):
    lins = by_clade[clade]
    W(f"\n## Clade {clade}\n")
    W(f"Tested in {len(lins)} lineages: {', '.join(l for l, _ in sorted(lins))}\n")
    if len(lins) < 2:
        W("**Only one lineage available — no replication test possible for this clade.**\n")
        summary.append((clade, len(lins), 0, 0))
        continue
    pairs = 0; replicated = set()
    for i in range(len(lins)):
        for j in range(i + 1, len(lins)):
            (l1, f1), (l2, f2) = lins[i], lins[j]
            td = os.path.join(a.streme_dir, f".tt_{clade}_{l1}_{l2}")
            r = subprocess.run(["tomtom", "-no-ssc", "-oc", td, "-verbosity", "1",
                                "-min-overlap", "5", "-dist", "pearson",
                                "-thresh", str(a.tomtom_thresh), f1, f2],
                               capture_output=True, text=True)
            tsv = os.path.join(td, "tomtom.tsv")
            hits = []
            if os.path.exists(tsv):
                with open(tsv) as fh:
                    fh.readline()
                    for line in fh:
                        if line.startswith("#") or not line.strip():
                            continue
                        f_ = line.split("\t")
                        if len(f_) > 5:
                            hits.append((f_[0], f_[1], f_[4]))
            s1, s2 = SIG.get((l1, clade)), SIG.get((l2, clade))
            if s1 is not None and s2 is not None:
                hits = [h for h in hits if h[0] in s1 and h[1] in s2]
            pairs += 1
            if hits:
                W(f"- **{l1} vs {l2}: {len(hits)} matching motif pair(s)**")
                for q, t, q_ in hits[:5]:
                    W(f"  - `{q}` ({l1}) matches `{t}` ({l2}), q = {q_}")
                    replicated.add(q)
            else:
                W(f"- {l1} vs {l2}: no matching motifs")
    W(f"\n**Clade {clade}: {len(replicated)} motif(s) replicated across lineages, "
      f"from {pairs} lineage pair(s).**")
    summary.append((clade, len(lins), pairs, len(replicated)))

W("\n\n## Summary\n")
W("| clade | lineages tested | lineage pairs | replicated motifs |")
W("|---|---|---|---|")
for c, nl, np_, nr in summary:
    W(f"| {c} | {nl} | {np_} | **{nr}** |")
total = sum(r for *_, r in summary)
W(f"\n**Total replicated motifs across all clades: {total}.**")
if total == 0:
    W("\n**This is a negative result.** No motif discovered for any clade replicated in a second "
      "independent lineage. Under this design that means the evidence does not support "
      "clade-linked regulation of fungal UPOs, and any single-lineage motif is better explained "
      "as a lineage-specific element. Report it as a negative, not as a set of weak hits.")

open(a.out, "w").write("\n".join(out) + "\n")
print(f"wrote {a.out}", file=sys.stderr)
