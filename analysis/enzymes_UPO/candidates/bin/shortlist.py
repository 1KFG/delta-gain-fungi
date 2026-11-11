#!/usr/bin/env python3
"""Expression shortlist: UPO candidates from lineages with no (or one) characterized enzyme.

  shortlist.py --placed upo_place_bfd_all.tsv.zst --scored upo_sites_bfd_scored.tsv.zst \
      --geneqc geneqc.tsv --asm asm.tsv --clusters70 clusters_070.tsv.zst \
      --clusters90 clusters_090.tsv.zst --nr reps090_vs_nr.tsv.zst \
      --panel assay_panel_96.tsv --fasta all_pf01328.faa --out candidates/out

The selection rules below were fixed before the script was run (CANDIDATES.md, "Rules").
"Novel" here means LINEAGE novelty: a member of a supported clade in which no enzyme has been
characterized, or only one. It is not sequence novelty; every catalytic-core representative has a
public relative in nr (G1_FINDINGS Finding 5).

Eligibility, every pool:
  - placement best edge not flagged `long_pendant`;
  - complete gene model in BFD (`is_partial` false, start and stop codon present);
  - length within the 5th-95th percentile of the pool's eligible proteins (drops fragments and
    fusions);
  - pools A and B: catalytic core (Cys39, Glu109, Asp111, Glu180) and placement LWR >= 0.95 for
    the pool's clade; pool C: Cys39, Glu109, Asp111 (Glu180 is the open question) and LWR >= 0.95.

Pools and quotas: see POOLS. Selection inside a pool: one protein per 70 % identity cluster
(catalytic-core clusters from the G1 dereplication; proteins outside that set are their own
cluster), at most one per genus, clusters taken largest first, round-robin across lineages
(LINEAGE_GROUP) so no class dominates. Within a cluster the representative is the member from a
RefSeq (GCF_) genome if there is one, then the highest LWR, then the protein ID.
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import subprocess
import sys
from collections import Counter, defaultdict

POOLS = [
    # (pool, label, quota, predicate on a row, which LWR column, needs Glu180)
    ("A", "Clade IV (one characterized member, MphUPO)", 10,
     lambda r: r["place_our_clade"] == "Clade IV", "place_our_lwr", True),
    ("B1", "Ascomycetous 1st (no characterized member)", 5,
     lambda r: r["place_fig1_clade"] == "Ascomycetous 1st", "place_fig1_lwr", True),
    ("B2", "Basidiomycetous 2nd (no characterized member)", 5,
     lambda r: r["place_fig1_clade"] == "Basidiomycetous 2nd", "place_fig1_lwr", True),
    ("B3", "Family I (short), outside Clades I-V (no characterized member)", 5,
     lambda r: (r["place_fig1_clade"] == "Ascomycetous 2nd (Family I, short)"
                and r["place_our_clade"] == "unassigned"), "place_fig1_lwr", True),
    ("B4", "Glomeromycete clade (no characterized member)", 3,
     lambda r: r["place_fig1_clade"].startswith("Glomeromycetous"), "place_fig1_lwr", True),
    ("B5", "Mucorales basal clade (no characterized member)", 3,
     lambda r: r["place_fig1_clade"] == "Mucoromycetous basal", "place_fig1_lwr", True),
    ("C1", "Ascomycete long UPOs, Clade I (Glu180 replaced; contingent on structure check)", 3,
     lambda r: r["place_our_clade"] == "Clade I" and r.get("PHYLUM") == "Ascomycota",
     "place_our_lwr", False),
    ("C2", "Basidiomycetous 1st (Glu180 replaced; contingent on structure check)", 3,
     lambda r: r["place_fig1_clade"] == "Basidiomycetous 1st", "place_fig1_lwr", False),
]


def text(path):
    if path.endswith(".zst"):
        return io.StringIO(subprocess.run(["zstd", "-dcq", path], capture_output=True, text=True,
                                          check=True).stdout)
    return open(path, newline="")


def tsv(path, header=True):
    fh = text(path)
    lines = (l for l in fh if not l.startswith("#"))
    return list(csv.DictReader(lines, delimiter="\t")) if header else [
        l.rstrip("\n").split("\t") for l in lines]


def pid(x):
    return x.split()[0]


def core3(r):
    return (r.get("p39_Cys_proximal") == "C" and r.get("p109_E_of_EHD") == "E"
            and r.get("p111_D_of_EHD") == "D")


def core4(r):
    return core3(r) and r.get("p180_Acid") == "E"


def fnum(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def percentile(vals, q):
    v = sorted(vals)
    if not v:
        return 0
    k = (len(v) - 1) * q
    lo = int(k)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (k - lo)


def eligible(rows, pred, lwr_col, need_glu, geneqc):
    pool = []
    reasons = Counter()
    for r in rows:
        if not pred(r):
            continue
        reasons["in pool"] += 1
        if r.get("place_flags"):
            reasons["long_pendant"] += 1
            continue
        q = geneqc.get(pid(r["protein_id"]))
        if q is None or q["is_partial"] != "false" or q["has_start_codon"].upper() != "TRUE" \
                or q["has_stop_codon"].upper() != "TRUE":
            reasons["incomplete gene model"] += 1
            continue
        # need_glu True: full core (Glu180 present). False: pool C, the lineage trait itself, so
        # Cys39/Glu109/Asp111 present and Glu180 ABSENT (deviation 1 in CANDIDATES.md).
        ok = core4(r) if need_glu else (core3(r) and r.get("p180_Acid") != "E")
        if not ok:
            reasons["catalytic residues"] += 1
            continue
        if fnum(r.get(lwr_col)) < 0.95:
            reasons["LWR < 0.95"] += 1
            continue
        pool.append(r)
    lens = [int(r["length"]) for r in pool]
    lo, hi = percentile(lens, 0.05), percentile(lens, 0.95)
    kept = [r for r in pool if lo <= int(r["length"]) <= hi]
    reasons["length outside 5-95 %"] = len(pool) - len(kept)
    reasons["eligible"] = len(kept)
    return kept, reasons, (lo, hi)


def select(pool, quota, lwr_col, cluster_of, cluster_size, gcf):
    """One per 70 % cluster, at most one per genus, largest clusters first, round-robin over
    lineages. Returns the chosen rows in selection order."""
    by_cluster = defaultdict(list)
    for r in pool:
        by_cluster[cluster_of.get(pid(r["protein_id"]), pid(r["protein_id"]))].append(r)

    def rep_key(r):
        return (not gcf(r), -fnum(r.get(lwr_col)), pid(r["protein_id"]))
    reps = []
    for c, members in by_cluster.items():
        best = sorted(members, key=rep_key)[0]
        reps.append((cluster_size.get(c, len(members)), c, best))
    reps.sort(key=lambda x: (-x[0], rep_key(x[2])))
    by_lineage = defaultdict(list)
    for size, c, r in reps:
        by_lineage[r.get("LINEAGE_GROUP") or "?"].append((size, c, r))
    lineages = sorted(by_lineage, key=lambda k: -by_lineage[k][0][0])
    chosen, genera = [], set()
    while len(chosen) < quota and any(by_lineage.values()):
        for lg in lineages:
            while by_lineage[lg]:
                size, c, r = by_lineage[lg].pop(0)
                genus = (r.get("SPECIES") or "?").split()[0]
                if genus in genera:
                    continue
                genera.add(genus)
                chosen.append((size, c, r))
                break
            if len(chosen) >= quota:
                break
    return chosen


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("placed", "scored", "geneqc", "asm", "clusters70", "clusters90", "nr", "panel",
              "fasta", "out"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    rows = tsv(a.placed)
    scored = {pid(r["protein_id"]): r for r in tsv(a.scored)}
    geneqc = {r["transcript_id"]: r for r in tsv(a.geneqc)}
    asm = {r["LOCUSTAG"]: r["ASMID"] for r in tsv(a.asm)}

    def gcf(r):
        return asm.get(pid(r["protein_id"]).split("_")[0], "").startswith("GCF_")
    c70, size70 = {}, Counter()
    for rep, mem in tsv(a.clusters70, header=False):
        c70[mem] = rep
        size70[rep] += 1
    c90 = {mem: rep for rep, mem in tsv(a.clusters90, header=False)}
    best_nr = {}
    for f in tsv(a.nr, header=False):
        best_nr.setdefault(f[0], (float(f[2]), f[7] if len(f) > 7 else ""))
    panel = set()
    for r in tsv(a.panel):
        panel.add(pid(r.get("protein_id", "")))

    out_rows, summary = [], []
    for code, label, quota, pred, lwr_col, need_glu in POOLS:
        pool, reasons, (lo, hi) = eligible(rows, pred, lwr_col, need_glu, geneqc)
        chosen = select(pool, quota, lwr_col, c70, size70, gcf)
        summary.append((code, label, quota, reasons, lo, hi, len(chosen)))
        for i, (size, c, r) in enumerate(chosen, 1):
            p = pid(r["protein_id"])
            rep90 = c90.get(p, p)
            nr_id, nr_sp = best_nr.get(rep90, ("", ""))
            s = scored.get(p, {})
            out_rows.append({
                "pool": code, "rank_in_pool": i, "protein_id": p,
                "species": r.get("SPECIES", ""), "lineage": r.get("LINEAGE_GROUP", ""),
                "phylum": r.get("PHYLUM", ""), "our_clade": r["place_our_clade"],
                "fig1_clade": r["place_fig1_clade"], "lwr": r.get(lwr_col, ""),
                "pendant_length": r.get("place_pendant_length", ""), "length": r["length"],
                "genome": asm.get(p.split("_")[0], ""),
                "refseq_genome": "yes" if gcf(r) else "no",
                "cluster70_size": size, "nr_best_identity_of_90rep": nr_id,
                "nr_best_species": nr_sp, "p180": r.get("p180_Acid", ""),
                "pocket_residues": r.get("pocket_residues", ""),
                "t1_pocket_matches": s.get("t1_pocket_matches", ""),
                "in_skylark_panel": "yes" if p in panel else "no",
            })
    cols = list(out_rows[0]) if out_rows else ["pool"]
    with open(os.path.join(a.out, "shortlist.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    want = {r["protein_id"] for r in out_rows}
    with open(os.path.join(a.out, "shortlist.faa"), "w") as fo:
        name, keep = None, False
        for line in text(a.fasta):
            if line.startswith(">"):
                keep = pid(line[1:]) in want
            if keep:
                fo.write(line)
    with open(os.path.join(a.out, "selection_summary.tsv"), "w") as fh:
        keys = ["in pool", "long_pendant", "incomplete gene model", "catalytic residues",
                "LWR < 0.95", "length outside 5-95 %", "eligible"]
        fh.write("pool\tlabel\tquota\t" + "\t".join(keys) + "\tlength_lo\tlength_hi\tselected\n")
        for code, label, quota, reasons, lo, hi, n in summary:
            fh.write(f"{code}\t{label}\t{quota}\t" + "\t".join(str(reasons[k]) for k in keys)
                     + f"\t{lo:.0f}\t{hi:.0f}\t{n}\n")
    print(f"shortlist: {len(out_rows)} candidates in {len(POOLS)} pools -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
