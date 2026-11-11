#!/usr/bin/env python3.12
"""
Turn the 5,115-protein shortlist into an assay-ready panel.

WHY A PANEL AND NOT A RANKING
-----------------------------
The score is a clade-V-like-pocket detector. Its measured AUC is 0.956 overall but 0.556 once
clade V is removed, on one out-of-clade active. Ranking by score alone would therefore order the
panel by how clade-V-like each protein is, and would fill a plate with near-duplicates of enzymes
that are already characterized. That tests nothing.

The panel is instead built for INFORMATION. Every enzyme assayed should be able to change what we
believe. That means spreading the panel across distinct pocket chemistries, clades and lineages,
and deliberately over-weighting the cases the current model is WORST at, because those are the
ones whose results carry information:

  - proteins OUTSIDE clade V, where out-of-clade AUC is 0.556 and the model is untested;
  - early-diverging phyla, which no characterized enzyme represents at all;
  - pocket residue combinations not present in any of the 19 characterized enzymes.

A panel that confirms the model teaches little. A panel that can break it is worth running.

SELECTION
---------
Greedy max-min diversity. Each protein is a feature vector of (pocket 5-mer, clade, phylum,
lineage). Starting from the highest-priority seed, repeatedly add the candidate that is most
dissimilar to everything already chosen, with priority as a tiebreak. This fills the plate with
spread rather than with the top of a ranked list.

Priority score, stated so it can be argued with:
  +3  clade call is not V            (where the model is untested)
  +3  early-diverging phylum         (no characterized representative exists)
  +2  pocket 5-mer unseen in the 19  (new chemistry, not a re-test)
  +1  clade_trust is reliable        (we know what we are assaying)
  +1  5/5 pocket matches             (strongest prior, worth including some)

DELIBERATE CONTROLS
-------------------
The panel reserves slots for NEGATIVE controls: proteins with the catalytic core but a LOW pocket
score. The model says these are inactive. The one out-of-clade active in the calibration set
(GmaUPO-II) also had a low score, so this is the cell where the model is most likely wrong, and
it cannot be tested without assaying some. A panel of predicted positives alone cannot measure a
false-negative rate.
"""
from __future__ import annotations
import argparse, subprocess, sys
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--scored", required=True, help="upo_sites_bfd_scored.tsv.zst")
ap.add_argument("--s2", required=True, help="active_site_S2.tsv, for the known pocket 5-mers")
ap.add_argument("--size", type=int, default=96, help="panel size (one plate)")
ap.add_argument("--negatives", type=int, default=12,
                help="slots reserved for predicted-INACTIVE controls")
ap.add_argument("--out", required=True)
ap.add_argument("--report", default=None)
a = ap.parse_args()

EARLY = {"Chytridiomycota", "Mucoromycota", "Basidiobolomycota", "Glomeromycota",
         "Zoopagomycota", "Blastocladiomycota", "Neocallimastigomycota", "Olpidiomycota",
         "Entomophthoromycota", "Kickxellomycota", "Monoblepharomycota", "Cryptomycota",
         "Microsporidia", "Aphelidiomycota", "Sanchytriomycota"}

# pocket 5-mers already represented among the 19 characterized enzymes
known = set()
cols = ["p73", "p77", "p81", "p183", "p187"]
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
    aa = [f[ix[c]].strip() for c in cols]
    if all(x and x != "-" for x in aa):
        known.add("".join(aa).upper())
print(f"{len(known)} distinct pocket 5-mers among the characterized enzymes: "
      f"{sorted(known)}", file=sys.stderr)

txt = subprocess.run(["zstd", "-dcq", a.scored], capture_output=True, text=True).stdout
hdr = None
cands, negs = [], []
for line in txt.splitlines():
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if hdr is None:
        hdr = f
        ix = {k: j for j, k in enumerate(hdr)}
        continue
    if f[ix["catalytic_core"]] != "yes":
        continue
    m = int(f[ix["t1_pocket_matches"]])
    pr = f[ix["pocket_residues"]].strip().upper()
    cl = f[ix["clade"]]
    ph = f[ix["PHYLUM"]]
    rec = dict(pid=f[ix["protein_id"]], sample=f[ix["sample"]], clade=cl,
               trust=f[ix["clade_trust"]], matches=m, pocket=pr, phylum=ph,
               lineage=f[ix["LINEAGE_GROUP"]], species=f[ix["SPECIES"]],
               ident=f[ix["identity_to_nearest_ref"]])
    if m >= 3:
        pr_score = 0
        if cl != "V":
            pr_score += 3
        if ph in EARLY:
            pr_score += 3
        if pr not in known:
            pr_score += 2
        if f[ix["clade_trust"]].startswith("reliable"):
            pr_score += 1
        if m == 5:
            pr_score += 1
        rec["priority"] = pr_score
        cands.append(rec)
    elif m <= 1:
        rec["priority"] = (3 if cl != "V" else 0) + (3 if ph in EARLY else 0)
        negs.append(rec)

print(f"{len(cands)} positive candidates, {len(negs)} low-score control candidates",
      file=sys.stderr)


def pick(pool, k):
    """Greedy max-min diversity over (pocket, clade, phylum, lineage), priority as tiebreak."""
    if not pool:
        return []
    pool = sorted(pool, key=lambda r: (-r["priority"], r["pid"]))
    chosen = [pool[0]]
    rest = pool[1:]
    while len(chosen) < k and rest:
        best, bi = None, None
        for i, r in enumerate(rest):
            d = min(sum((r["pocket"] != c["pocket"], r["clade"] != c["clade"],
                         r["phylum"] != c["phylum"], r["lineage"] != c["lineage"]))
                    for c in chosen)
            key = (d, r["priority"])
            if best is None or key > best:
                best, bi = key, i
        chosen.append(rest.pop(bi))
    return chosen


n_pos = a.size - a.negatives
panel = pick(cands, n_pos)
controls = pick(negs, a.negatives)
for r in panel:
    r["role"] = "predicted_active"
for r in controls:
    r["role"] = "predicted_inactive_control"
allp = panel + controls

cols_out = ["rank", "role", "protein_id", "sample", "species", "PHYLUM", "LINEAGE_GROUP",
            "clade", "clade_trust", "identity_to_nearest_ref", "t1_pocket_matches",
            "pocket_residues", "pocket_is_novel", "priority"]
with open(a.out, "w") as fh:
    fh.write("# Assay panel from the BFD UPO shortlist.\n")
    fh.write("# Built for INFORMATION, not for confirmation: weighted toward the cases where the\n")
    fh.write("#   score is untested (outside clade V, out-of-clade AUC 0.556) and toward pocket\n")
    fh.write("#   combinations absent from the 19 characterized enzymes.\n")
    fh.write(f"# {len(panel)} predicted active + {len(controls)} predicted INACTIVE controls.\n")
    fh.write("# The controls are not filler. Without them no false-negative rate can be measured,\n")
    fh.write("#   and GmaUPO-II shows that an active enzyme CAN carry a low pocket score.\n")
    fh.write("\t".join(cols_out) + "\n")
    for i, r in enumerate(allp, 1):
        fh.write("\t".join([str(i), r["role"], r["pid"], r["sample"], r["species"], r["phylum"],
                            r["lineage"], r["clade"], r["trust"], r["ident"],
                            str(r["matches"]), r["pocket"],
                            "yes" if r["pocket"] not in known else "no",
                            str(r["priority"])]) + "\n")
print(f"wrote {a.out}", file=sys.stderr)

rep = [f"# Assay panel: {len(allp)} proteins\n",
       f"{len(panel)} predicted active, {len(controls)} predicted inactive controls.\n",
       "## Predicted-active slots\n", "| property | n | % of slots |", "|---|---|---|"]
nv = sum(1 for r in panel if r["clade"] != "V")
ed = sum(1 for r in panel if r["phylum"] in EARLY)
nov = sum(1 for r in panel if r["pocket"] not in known)
rel = sum(1 for r in panel if r["trust"].startswith("reliable"))
for lab, v in [("outside clade V", nv), ("early-diverging phylum", ed),
               ("pocket 5-mer not in the 19", nov), ("reliable clade call", rel)]:
    rep.append(f"| {lab} | {v} | {100*v/len(panel):.0f} |")
rep += ["\n## Spread\n",
        f"- distinct pocket 5-mers: {len({r['pocket'] for r in panel})}",
        f"- distinct species: {len({r['species'] for r in panel})}",
        f"- distinct lineages: {len({r['lineage'] for r in panel})}", ""]
rep += ["## Clade composition of the active slots\n", "| clade | n |", "|---|---|"]
for k, v in Counter(r["clade"] for r in panel).most_common():
    rep.append(f"| {k} | {v} |")
rep += ["\n## Phylum composition\n", "| phylum | n |", "|---|---|"]
for k, v in Counter(r["phylum"] for r in panel).most_common():
    rep.append(f"| {k} | {v} |")
rep += ["\n## Why the controls matter\n",
        "The model predicts these are inactive. GmaUPO-II is a characterized enzyme that is 49.3 %",
        "active on 4-propylguaiacol and carries only 1 of 5 pocket matches, so the predicted-",
        "inactive cell is exactly where the model is known to fail. Assaying only predicted",
        "positives would measure precision and leave the false-negative rate unknown.\n"]
if a.report:
    open(a.report, "w").write("\n".join(rep) + "\n")
    print(f"wrote {a.report}", file=sys.stderr)
print("\n".join(rep), file=sys.stderr)
