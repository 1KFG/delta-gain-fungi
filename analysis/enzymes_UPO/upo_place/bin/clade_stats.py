#!/usr/bin/env python3
"""Per-clade statistics under placement clade calls, for the G1 write-up and for revisiting the
numbers in COORDINATION.md / SYNTHESIS_REVIEW.md. Prints markdown tables.

  clade_stats.py --table upo_place_bfd_all.tsv.zst \
      [--nr novelty/reps090_vs_nr.tsv.zst] [--loo placement_ref/loo.tsv] > CLADE_STATS.md

Identity columns (`clade_identity`, `nearest_ref`) are the old upo_sites values: identity to the
nearest of the 19 characterized enzymes over PF01328 match columns. They are reported here because
the earlier documents argued from them; they no longer decide any call.

Novelty (`--nr`): DIAMOND blastp of the 90 %-dereplicated catalytic-core representatives against
nr; best hit = first row per query (DIAMOND orders by e-value); "novel" = best-hit identity < 70 %,
the cut used in G1_FINDINGS.md Finding 5. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import io
import statistics
import subprocess
from collections import Counter, defaultdict

CLADES = ["Clade I", "Clade II", "Clade III", "Clade IV", "Clade V", "unassigned"]


def lines(path):
    if path.endswith(".zst"):
        return io.StringIO(subprocess.run(["zstd", "-dcq", path], capture_output=True, text=True,
                                          check=True).stdout)
    return open(path, newline="")


def is_core(r):
    return (r.get("p39_Cys_proximal") == "C" and r.get("p109_E_of_EHD") == "E"
            and r.get("p111_D_of_EHD") == "D" and r.get("p180_Acid") == "E")


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def md(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


def pct(a, b):
    return f"{100 * a / b:.1f} %" if b else "-"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", required=True)
    ap.add_argument("--nr")
    ap.add_argument("--loo")
    a = ap.parse_args(argv)
    rows = list(csv.DictReader((l for l in lines(a.table) if not l.startswith("#")),
                               delimiter="\t"))
    by = defaultdict(list)
    for r in rows:
        by[r["place_our_clade"]].append(r)
    n = len(rows)
    out = [f"# Clade statistics under placement calls (n = {n})\n"]

    out.append("## A. Per clade: old identity metric, length, Glu180, catalytic core\n")
    t = []
    for c in CLADES:
        rs = by[c]
        ids = [v for v in (fnum(r["clade_identity"]) for r in rs) if v is not None]
        lens = [int(r["length"]) for r in rs]
        core = [r for r in rs if is_core(r)]
        hi = sum(1 for v in ids if v >= 52.7)
        t.append([c, len(rs), f"{statistics.mean(ids):.1f}" if ids else "-",
                  f"{statistics.median(ids):.1f}" if ids else "-", hi,
                  f"{statistics.mean(lens):.1f}", pct(sum(r.get('p180_Acid') == 'E' for r in rs),
                                                       len(rs)),
                  len(core), sum(1 for r in core if fnum(r["place_our_lwr"]) >= 0.95)])
    out.append(md(["clade", "n", "mean id to nearest of 19", "median id", "n id >= 52.7",
                   "mean length", "Glu180", "catalytic core", "core, LWR >= 0.95"], t))

    out.append("## B. Glu180 by call confidence\n")
    tiers = defaultdict(list)
    for r in rows:
        if r["place_our_clade"] != "unassigned":
            k = "assigned LWR >= 0.95" if fnum(r["place_our_lwr"]) >= 0.95 else "assigned 0.5-0.95"
        else:
            k = "unassigned"
        tiers[k].append(r)
    out.append(md(["call", "n", "Glu180"],
                  [[k, len(v), pct(sum(r.get("p180_Acid") == "E" for r in v), len(v))]
                   for k, v in tiers.items()]))
    old_val = [r for r in rows if r["clade_confidence"].startswith("validated")]
    out.append(f"Old validated-range calls (identity >= 52.7): n = {len(old_val)}, Glu180 "
               f"{pct(sum(r.get('p180_Acid') == 'E' for r in old_val), len(old_val))}.\n")

    out.append("## C. Clade IV\n")
    iv = by["Clade IV"]
    lg = Counter(r.get("LINEAGE_GROUP") or "?" for r in iv)
    lgc = Counter(r.get("LINEAGE_GROUP") or "?" for r in iv if is_core(r))
    out.append(md(["lineage", "Clade IV", "catalytic core"],
                  [[k, v, lgc[k]] for k, v in lg.most_common(12)]))
    nref = Counter(r["nearest_ref"] for r in iv)
    out.append("Old nearest reference of Clade IV proteins: "
               + ", ".join(f"{k} {v}" for k, v in nref.most_common(6)) + "\n")
    oldc = Counter(r["clade"] or "unassigned" for r in iv)
    out.append("Old clade call of Clade IV proteins: "
               + ", ".join(f"{k} {v}" for k, v in oldc.most_common()) + "\n")
    allm = Counter(r["nearest_ref"] for r in rows)
    out.append(f"MphUPO is the old nearest reference for {allm['MphUPO']} proteins overall.\n")

    out.append("## D. Lineage x clade, Clade I by phylum\n")
    ph = defaultdict(Counter)
    for r in by["Clade I"]:
        ph[r.get("PHYLUM") or "?"][r.get("p180_Acid") == "E"] += 1
    out.append(md(["phylum", "Clade I", "Glu180"],
                  [[k, v[True] + v[False], pct(v[True], v[True] + v[False])]
                   for k, v in sorted(ph.items(), key=lambda x: -sum(x[1].values()))]))
    lin = Counter(r.get("LINEAGE_GROUP") or "?" for r in by["Clade I"])
    out.append("Clade I by lineage: " + ", ".join(f"{k} {v}" for k, v in lin.most_common(8))
               + "\n")

    out.append("## E. Outside Clades I-V: Figure 1 lineage x catalytic core\n")
    un = by["unassigned"]
    f1 = Counter(r["place_fig1_clade"] for r in un)
    f1c = Counter(r["place_fig1_clade"] for r in un if is_core(r))
    out.append(md(["Figure 1 clade (unassigned in I-V)", "n", "catalytic core"],
                  [[k, v, f1c[k]] for k, v in f1.most_common()]))
    sl = Counter(r.get("LINEAGE_GROUP") or "?" for r in un)
    out.append("Unassigned by lineage: " + ", ".join(f"{k} {v}" for k, v in sl.most_common(10))
               + "\n")

    if a.nr:
        out.append("## F. Novelty vs nr by placement clade (90 % reps of the catalytic core)\n")
        best = {}
        for l in lines(a.nr):
            f = l.rstrip("\n").split("\t")
            if f and f[0] not in best:
                best[f[0]] = float(f[2])
        info = {r["protein_id"].split()[0]: r for r in rows}
        nc, nf = defaultdict(lambda: [0, 0]), defaultdict(lambda: [0, 0])
        miss = 0
        for q, pid in best.items():
            r = info.get(q)
            if r is None:
                miss += 1
                continue
            for d, k in ((nc, r["place_our_clade"]), (nf, r["place_fig1_clade"])):
                d[k][1] += 1
                d[k][0] += pid < 70
        tot = [sum(v[0] for v in nc.values()), sum(v[1] for v in nc.values())]
        out.append(f"Representatives with a hit: {len(best)} (not in table: {miss}); novel "
                   f"(< 70 %) overall {tot[0]}/{tot[1]} = {pct(*tot)}.\n")
        out.append(md(["clade", "novel (< 70 %)", "n"],
                      [[k, pct(*v), v[1]] for k, v in sorted(nc.items(), key=lambda x: -x[1][1])]))
        out.append(md(["Figure 1 clade", "novel (< 70 %)", "n"],
                      [[k, pct(*v), v[1]] for k, v in sorted(nf.items(), key=lambda x: -x[1][1])]))

    if a.loo:
        out.append("## G. Leave-one-out by expected clade\n")
        res = defaultdict(lambda: [0, 0])
        with open(a.loo, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                res[r["our_expected"]][1] += 1
                res[r["our_expected"]][0] += r["our_ok"] == "yes"
        out.append(md(["expected clade", "correct", "n"],
                      [[k, v[0], v[1]] for k, v in sorted(res.items())]))
    print("\n".join(out))


if __name__ == "__main__":
    main()
