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

VALID = {"validated_range", "validated_low_margin"}
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
n_valid = sum(v for k, v in conf.items() if k in VALID)
W(f"\n**Only {n_valid} of {len(rows)} ({100*n_valid/len(rows):.2f} %) clade calls fall inside the "
  f"validated identity range.** Everything else is a hypothesis.\n")

W("## 2. Clade counts\n")
cl_all = Counter(r["clade"] for r in rows)
cl_val = Counter(r["clade"] for r in rows if r["clade_confidence"] in VALID)
W("| clade | all calls | within validated range | in-sample T1 (n=19) |")
W("|---|---|---|---|")
T1 = {"I": "0-53, bimodal", "II": "0", "III": "0-3", "IV": "0-5", "V": "33-82",
      "unassigned": "n/a"}
for c, v in cl_all.most_common():
    W(f"| {c} | {v} | {cl_val.get(c,0)} | {T1.get(c,'?')} |")

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

W("\n## 6. Copy number per genome\n")
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
