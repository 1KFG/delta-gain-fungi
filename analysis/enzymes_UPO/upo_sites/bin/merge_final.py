#!/usr/bin/env python3.12
"""
Build the final per-protein table for every PF01328 protein in BFD.

Merges the three layers, each carrying its own validated reliability, so a consumer cannot
accidentally treat a bounded call as a validated one:

  1. Active-site residues and catalytic-core status   -- VALIDATED (95/95 vs the published
     structure-based matrix; catalytic positions rank 6/12/15 of 186 match states by conservation)
  2. Clade by phylogenetic placement, with likelihood weight -- BOUNDED (100 % precision above
     70 % identity, 88 % at 55-70 %, 38 % at 45-55 %; BFD median is 43.8 %)
  3. Lineage and species

It deliberately emits NO activity score. See MODEL_VALIDATION.md section D.

The `clade_trust` column is the one to filter on. `reliable` means the placement sits above the
identity at which the caller was actually validated.
"""
from __future__ import annotations
import argparse, glob, os, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True, help="upo_sites_bfd_all.tsv.zst")
ap.add_argument("--placement-dir", required=True, help="data/placement/full")
ap.add_argument("--reps-placement", default=None,
                help="placement_clades.tsv for the 7,055 representatives (fallback/cross-check)")
ap.add_argument("--validated-identity", type=float, default=55.0)
ap.add_argument("--out", required=True)
ap.add_argument("--report", default=None)
a = ap.parse_args()


def read_zst(p):
    return subprocess.run(["zstd", "-dcq", p], capture_output=True, text=True).stdout


place = {}
files = sorted(glob.glob(os.path.join(a.placement_dir, "*.clades.tsv")))
for f in files:
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {k: j for j, k in enumerate(hdr)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) < 3:
                continue
            place[r[i["query"]]] = (r[i["clade"]], float(r[i["placement_weight"]]),
                                    r[i["note"]] if "note" in i else "")
print(f"placement: {len(place)} proteins from {len(files)} chunks", file=sys.stderr)

if a.reps_placement and os.path.exists(a.reps_placement):
    n_add = 0
    with open(a.reps_placement) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {k: j for j, k in enumerate(hdr)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) >= 3 and r[i["query"]] not in place:
                place[r[i["query"]]] = (r[i["clade"]], float(r[i["placement_weight"]]), "")
                n_add += 1
    if n_add:
        print(f"  + {n_add} from representative placement", file=sys.stderr)

txt = read_zst(a.sites)
lines = txt.splitlines()
hdr = None
out_rows = []
from collections import Counter
tally = Counter(); trust_tally = Counter()

for line in lines:
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if hdr is None:
        hdr = f; ix = {k: j for j, k in enumerate(hdr)}
        continue
    pid = f[ix["protein_id"]]
    core = (f[ix.get("p39_Cys_proximal", 0)] == "C" and f[ix.get("p109_E_of_EHD", 0)] == "E"
            and f[ix.get("p111_D_of_EHD", 0)] == "D" and f[ix.get("p180_Acid", 0)] == "E")
    nn_clade = f[ix["clade"]]
    nn_id = float(f[ix["clade_identity"]])
    pl = place.get(pid)
    if pl:
        pclade, pw, _ = pl
        method = "placement"
    else:
        pclade, pw, method = "not_placed", 0.0, "none"
    if pclade in ("unassigned", "not_placed"):
        trust = "unassigned"
    elif nn_id >= a.validated_identity:
        trust = "reliable"
    else:
        trust = "extrapolated"
    if pclade == "IV":
        trust += "+provisional_cladeIV"
    tally[pclade] += 1; trust_tally[trust] += 1
    out_rows.append((pid, f[ix["sample"]], f[ix["length"]],
                     "yes" if core else "no",
                     pclade, f"{pw:.3f}", method, trust,
                     nn_clade, f"{nn_id:.1f}",
                     f[ix["pocket_score"]], f[ix["pocket_residues"]],
                     f[ix.get("PHYLUM", 0)], f[ix.get("LINEAGE_GROUP", 0)],
                     f[ix.get("SPECIES", 0)]))

cols = ["protein_id", "sample", "length", "catalytic_core", "clade", "placement_weight",
        "clade_method", "clade_trust", "nn_clade", "identity_to_nearest_ref",
        "pocket_score_of_5", "pocket_residues", "PHYLUM", "LINEAGE_GROUP", "SPECIES"]
buf = ["# Final per-protein table for all PF01328 proteins in BFD.",
       "# Site residues and catalytic_core are VALIDATED (95/95 vs the published structure-based",
       "#   matrix; catalytic positions rank 6/12/15 of 186 match states by conservation).",
       "# clade is BOUNDED: filter on clade_trust. 'reliable' = identity to nearest reference is at",
       f"#   or above {a.validated_identity} %, where the caller measured 88-100 % precision.",
       "#   'extrapolated' = below it, where precision measured 38-50 %. Treat as a hypothesis.",
       "# Clade IV always carries +provisional_cladeIV: it has ONE characterized member.",
       "# NO ACTIVITY SCORE is emitted; see MODEL_VALIDATION.md section D for why.",
       "\t".join(cols)]
buf += ["\t".join(r) for r in out_rows]

p = subprocess.Popen(["zstd", "-qf", "-o", a.out], stdin=subprocess.PIPE, text=True)
p.communicate("\n".join(buf) + "\n")
print(f"wrote {len(out_rows)} rows to {a.out}", file=sys.stderr)

rep = ["# Final BFD UPO table: composition\n",
       f"**{len(out_rows)} PF01328 proteins.**\n",
       "## Clade (phylogenetic placement)\n",
       "| clade | n |", "|---|---|"]
for k, v in tally.most_common():
    rep.append(f"| {k} | {v} |")
rep += ["\n## Trust tier — the column to filter on\n", "| tier | n | % |", "|---|---|---|"]
for k, v in trust_tally.most_common():
    rep.append(f"| `{k}` | {v} | {100*v/len(out_rows):.1f} |")
if a.report:
    open(a.report, "w").write("\n".join(rep) + "\n")
    print(f"wrote {a.report}", file=sys.stderr)
print("\n".join(rep[:40]), file=sys.stderr)
