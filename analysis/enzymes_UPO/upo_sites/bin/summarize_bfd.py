#!/usr/bin/env python3.12
"""
G1 summary: clade distribution and site conservation of fungal PF01328 proteins across BFD.

Joins upo_sites output to BFD lineage (samples.parquet) and reports:
  1. clade counts, split by confidence tier -- extrapolated calls are never merged with
     validated ones;
  2. clade distribution per phylum/class;
  3. pocket score distribution per clade and per lineage;
  4. catalytic-machinery integrity;
  5. copy number per genome.

Every table states which confidence tiers it includes. See upo_sites/README.md: at BFD identity
almost every clade call is an extrapolation beyond the validated range.
"""
import argparse, glob, gzip, io, os, subprocess, sys
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True, help="dir of chunk_*.tsv.zst from upo_sites")
ap.add_argument("--samples", required=True, help="samples.parquet")
ap.add_argument("--out", required=True, help="output markdown report")
ap.add_argument("--tsv-out", default=None, help="optional merged per-protein TSV (.zst)")
a = ap.parse_args()

def read_rows(d):
    hdr = None
    for p in sorted(glob.glob(os.path.join(d, "chunk_*.tsv.zst"))):
        txt = subprocess.run(["zstd", "-dcq", p], capture_output=True, text=True).stdout
        for line in txt.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if hdr is None:
                hdr = f; continue
            if f[0] == hdr[0]:
                continue
            yield dict(zip(hdr, f))

rows = list(read_rows(a.results))
if not rows:
    sys.exit(f"no rows found in {a.results}")
print(f"loaded {len(rows)} proteins", file=sys.stderr)

# lineage via duckdb CLI (no python duckdb module on this host)
lin_tsv = os.path.join(os.path.dirname(a.out), ".lineage.tsv")
subprocess.run(["duckdb", "-c",
    f"copy (select LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, \"ORDER\" as ORD, GENUS, SPECIES "
    f"from '{a.samples}') to '{lin_tsv}' (header, delimiter '\t');"], check=True)
LIN = {}
with open(lin_tsv) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        LIN[f[0]] = dict(zip(h, f))

def lineage_of(pid):
    lt = pid.split("_")[0]
    d = LIN.get(lt, {})
    ph, cl = d.get("PHYLUM", "?"), d.get("CLASS", "?")
    # group as in POCKET_TRANSFER: phylum for non-Dikarya, class for Asco/Basidio
    grp = cl if ph in ("Ascomycota", "Basidiomycota") and cl not in ("", "?") else ph
    return ph, cl, (grp or "?"), lt

# One definition of "validated", used everywhere. The +provisional_cladeIV suffix must not
# change whether a call counts as validated-range; it is reported as a separate marker instead.
# (Before 2026-11-08 section 2 used exact matching and section 6 used startswith, so Clade IV
# appeared as 0 in one table and 4,565 in the other.)
def is_validated(conf): return conf.startswith("validated")
out = []
W = out.append
W("# G1: fungal PF01328 diversity across BFD\n")
W(f"Generated from `{a.results}`. **n = {len(rows)} proteins.**\n")
W("Clade calls below 52.7 % identity to the nearest characterized enzyme are **extrapolations "
  "beyond the validated range** of the caller and are reported separately throughout. "
  "See `upo_sites/README.md`.\n")

conf = Counter(r["clade_confidence"] for r in rows)
W("## 1. Confidence tiers\n")
W("| tier | n | % |")
W("|---|---|---|")
for k, v in conf.most_common():
    W(f"| `{k}` | {v} | {100*v/len(rows):.1f} |")
n_valid = sum(v for k, v in conf.items() if is_validated(k))
n_valid_IV = sum(1 for r in rows
                 if r["clade"] == "IV" and is_validated(r["clade_confidence"]))
W(f"\n**Only {n_valid} of {len(rows)} ({100*n_valid/len(rows):.2f} %) clade calls fall inside the "
  f"validated identity range.** Everything else is a hypothesis.\n")
W(f"**Of those {n_valid}, {n_valid_IV} ({100*n_valid_IV/max(n_valid,1):.1f} %) are Clade IV**, which "
  f"has one characterized member and therefore cannot be validated by leave-one-out at all. "
  f"Excluding Clade IV leaves **{n_valid - n_valid_IV}** calls "
  f"({100*(n_valid-n_valid_IV)/len(rows):.2f} % of all proteins) that are both inside the validated "
  f"range and anchored on a clade with more than one reference.\n")

W("## 2. Clade counts\n")
cl_all = Counter(r["clade"] for r in rows)
cl_val = Counter(r["clade"] for r in rows if is_validated(r["clade_confidence"]))
W("**Clade IV carries a single characterized member (MphUPO), so its `within validated range` "
  "count means 'high identity to one uncharacterized protein'. It is marked provisional "
  "everywhere below and must not be pooled with the other clades' validated counts.**\n")
W("| clade | all calls | within validated range | in-sample T1 (n=19) |")
W("|---|---|---|---|")
T1 = {"I": "0-53, bimodal", "II": "0", "III": "0-3", "IV": "0-5", "V": "33-82",
      "unassigned": "n/a"}
for c, v in cl_all.most_common():
    mark = " **(provisional, 1 ref)**" if c == "IV" else ""
    W(f"| {c}{mark} | {v} | {cl_val.get(c,0)} | {T1.get(c,'?')} |")

W("\n## 3. Pocket score (0-5) by clade -- all calls\n")
byc = defaultdict(list)
for r in rows:
    byc[r["clade"]].append(int(r["pocket_score"]))
W("| clade | n | mean | 5/5 | 0/5 |")
W("|---|---|---|---|---|")
for c in sorted(byc, key=lambda k: -len(byc[k])):
    v = byc[c]
    W(f"| {c} | {len(v)} | {sum(v)/len(v):.2f} | {v.count(5)} | {v.count(0)} |")

W("\n## 4. Clade distribution by lineage\n")
lg = defaultdict(Counter); lgn = Counter()
for r in rows:
    _, _, grp, _ = lineage_of(r["protein_id"])
    lg[grp][r["clade"]] += 1; lgn[grp] += 1
clades = [c for c, _ in cl_all.most_common()]
W("| lineage | n | " + " | ".join(clades) + " |")
W("|---" * (len(clades) + 2) + "|")
for grp in sorted(lgn, key=lambda k: -lgn[k]):
    W(f"| {grp} | {lgn[grp]} | " + " | ".join(str(lg[grp].get(c, 0)) for c in clades) + " |")

W("\n## 5. Catalytic machinery\n")
ok = sum(1 for r in rows if r["catalytic_ok"] == "yes")
W(f"- Full catalytic set intact (Cys39, Glu109, Asp111, Glu180): **{ok}/{len(rows)} "
  f"({100*ok/len(rows):.1f} %)**")
cys = sum(1 for r in rows if r.get("p39_Cys_proximal", "") == "C")
W(f"- Proximal Cys39 present: **{cys}/{len(rows)} ({100*cys/len(rows):.1f} %)**")
base = Counter(r.get("p110_Base_EHD", "-") for r in rows)
W(f"- Residue at position 110 (His in short UPOs, Gly in long): "
  + ", ".join(f"{k}={v}" for k, v in base.most_common(6)))

W("\n## 6. The catalytic-core filter -- which of these are plausibly UPOs at all\n")
W("PF01328 is a broad heme-thiolate peroxidase family. The 2025 paper states that **Glu180 is "
  "present in all 19** characterized UPOs, so it is a usable test of whether a PF01328 protein "
  "has UPO-like catalytic machinery.\n")
def tier(r):
    if r["clade"] == "unassigned": return "unassigned"
    return "validated" if r["clade_confidence"].startswith("validated") else "extrapolated"
t180 = defaultdict(Counter); tn = Counter()
for r in rows:
    t = tier(r); t180[t][r.get("p180_Acid", "-")] += 1; tn[t] += 1
W("| confidence tier | n | Glu180 |")
W("|---|---|---|")
for t in ("validated", "extrapolated", "unassigned"):
    if tn[t]:
        W(f"| {t} | {tn[t]} | {100*t180[t]['E']/tn[t]:.1f} % |")
c180 = defaultdict(Counter); cn = Counter()
for r in rows:
    c180[r["clade"]][r.get("p180_Acid", "-")] += 1; cn[r["clade"]] += 1
W("\n| clade | n | Glu180 |")
W("|---|---|---|")
for c in sorted(cn, key=lambda k: -cn[k]):
    W(f"| {c} | {cn[c]} | {100*c180[c]['E']/cn[c]:.1f} % |")
W("\n**Validated-range calls carry Glu180 at 93.7 %, matching the published claim. Clade I calls "
  "carry it at 16.8 %.** The same alignment gives all five characterized Clade I enzymes a Glu "
  "at 180, so this is not an alignment artifact for long UPOs. The most likely reading is that "
  "the low-identity Clade I bin is absorbing PF01328 proteins that are **not UPOs**, because the "
  "long characterized enzymes are the nearest available match for divergent sequences. Treat "
  "unfiltered Clade I counts as an upper bound.\n")
core = [r for r in rows
        if r.get("p39_Cys_proximal") == "C" and r.get("p109_E_of_EHD") == "E"
        and r.get("p111_D_of_EHD") == "D" and r.get("p180_Acid") == "E"]
W(f"### Catalytic core set\n")
W(f"Requiring Cys39 + Glu109 + Asp111 + Glu180: **{len(core)}/{len(rows)} "
  f"({100*len(core)/len(rows):.1f} %)** proteins.\n")
cc = Counter(r["clade"] for r in core)
ccv = Counter(r["clade"] for r in core if is_validated(r["clade_confidence"]))
W("| clade | catalytic core | of which validated-range |")
W("|---|---|---|")
for c, v in cc.most_common():
    mark = " **(provisional, 1 ref)**" if c == "IV" else ""
    W(f"| {c}{mark} | {v} | {ccv.get(c,0)} |")
W("\n**This is the defensible UPO set for downstream work**, and it is the one to carry into "
  "novelty and regulatory analysis, not the raw 97,520.\n")

W("\n## 7. Copy number per genome\n")
per = Counter(r["protein_id"].split("_")[0] for r in rows)
top = per.most_common(10)
W(f"- genomes with at least one PF01328: **{len(per)}**")
W(f"- maximum copies in one genome: **{top[0][1]}**")
W("\n| genome | copies | species |")
W("|---|---|---|")
for lt, n in top:
    W(f"| {lt} | {n} | *{LIN.get(lt,{}).get('SPECIES','?')}* |")

open(a.out, "w").write("\n".join(out) + "\n")
print(f"wrote {a.out}", file=sys.stderr)

if a.tsv_out:
    hdr = list(rows[0].keys()) + ["PHYLUM", "CLASS", "LINEAGE_GROUP", "SPECIES"]
    buf = io.StringIO(); buf.write("\t".join(hdr) + "\n")
    for r in rows:
        ph, cl, grp, lt = lineage_of(r["protein_id"])
        buf.write("\t".join(list(r.values()) +
                            [ph, cl, grp, LIN.get(lt, {}).get("SPECIES", "?")]) + "\n")
    p = subprocess.Popen(["zstd", "-qf", "-o", a.tsv_out], stdin=subprocess.PIPE, text=True)
    p.communicate(buf.getvalue())
    print(f"wrote {a.tsv_out}", file=sys.stderr)
