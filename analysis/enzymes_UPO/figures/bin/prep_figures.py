#!/usr/bin/env python3
"""Build the data behind Figures 1, 3, 4 and Table T1 of the UPO clade summary.

  prep_figures.py --root analysis/enzymes_UPO --bfd <upo_place_bfd_all.tsv.zst> --out figures/data

Writes (all small TSVs, committed):
  fig1_tips.tsv        one row per reference leaf: clade labels and the Figure 1 colour group
  fig1_support.tsv     SH-aLRT / UFBoot of each named clade's node in the rooted mfp tree
  fig3_flows.tsv       old nearest-neighbour call (clade, tier) -> placement clade, counts
  fig4_precision.tsv   matched-identity outcomes by identity bin, both callers
  table_T1.tsv / .md   per-clade summary (Clades I-V, outside, and Zamocky Figure 1 lineages)

Inputs are the committed outputs of refset/, upo_place/ and the skylark session's placement
validation; the BFD table is the merged per-protein table on the shared project. Standard library
only; reuses the tree code in refset/bin.
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict

CLADES = ["Clade I", "Clade II", "Clade III", "Clade IV", "Clade V"]
BINS = ["70-100", "55-70", "45-55", "35-45", "<35"]

# Zamocky's 11 Figure 1 clades folded into 8 colour groups (the categorical cap). The three
# non-dikaryan minorities share one group; Family I and II are named by their size class.
FIG1_GROUP = {
    "Mucoromycetous basal": "Mucorales (basal)",
    'Glomeromycetous (Fig 1 "Mucoromycetous small")': "Glomeromycete",
    "Basidiomycetous 1st": "Basidiomycetous 1st",
    "Ascomycetous 1st": "Ascomycetous 1st",
    "Basidiomycetous 2nd": "Basidiomycetous 2nd",
    "Stramenopilous 1st": "Oomycete, chytrid, amoebozoan",
    "Stramenopilous 2nd": "Oomycete, chytrid, amoebozoan",
    "Chytridiomycetous": "Oomycete, chytrid, amoebozoan",
    "Amoebozoan": "Oomycete, chytrid, amoebozoan",
    "Basidiomycetous 3rd (Family II, long)": "Family II (long)",
    "Ascomycetous 2nd (Family I, short)": "Family I (short)",
}
MUCORALES = ["Z__AglauHTP__Mucoro", "Z__CcucHTP__Mucoro", "Z__MucambHTP__Mucoro",
             "Z__MuclusHTP__Mucoro", "Z__MuccirHTP__Mucoro", "Z__PpaHTP__Mucoro",
             "Z__RdelHTP__Mucoro", "Z__RstolHTP1__Mucoro", "Z__RmicroHTP1__Mucoro",
             "Z__RmicroHTP2__Mucoro"]


def tsv_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def write(path, rows, cols=None):
    cols = cols or list(rows[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def bin_of(x: float) -> str:
    return ("70-100" if x >= 70 else "55-70" if x >= 55 else "45-55" if x >= 45
            else "35-45" if x >= 35 else "<35")


# ---------- Figure 1 ----------

def fig1(members_tsv, tree_path, refset_bin):
    sys.path.insert(0, refset_bin)
    import name_clades as nc
    import tree_outline as to
    mem = tsv_rows(members_tsv)
    tips = [{"label": r["label"], "set": r["set"], "name": r["name"], "group": r["group"],
             "fig1_clade": r["fig1_clade"], "fig1_group": FIG1_GROUP.get(r["fig1_clade"], ""),
             "our_clade": r["our_clade"]} for r in mem]
    root = to.reroot(to.parse_newick(open(tree_path).read()), set(MUCORALES))
    leaves, below = nc.index(root)
    sets = defaultdict(set)
    for r in mem:
        if r["our_clade"]:
            sets[("our", r["our_clade"])].add(r["label"])
        if r["fig1_clade"]:
            sets[("fig1", r["fig1_clade"])].add(r["label"])
    sup = []
    for (system, clade), m in sorted(sets.items()):
        n = nc.mrca(leaves, m)
        sh, ufb = (n.label.split("/") + ["", ""])[:2] if n.label else ("", "")
        sup.append({"system": system, "clade": clade, "n_leaves": len(m),
                    "exact": "yes" if below[id(n)] == frozenset(m) else "no",
                    "sh_alrt": sh, "ufboot": ufb})
    return tips, sup


# ---------- Figure 3 ----------

def old_tier(conf: str) -> str:
    if conf.startswith("validated"):
        return "validated"
    if conf.startswith("below"):
        return "below floor"
    return "extrapolated"


def read_bfd(path):
    if path.endswith(".zst"):
        txt = subprocess.run(["zstd", "-dcq", path], capture_output=True, text=True,
                             check=True).stdout
        return list(csv.DictReader(io.StringIO(txt), delimiter="\t"))
    return tsv_rows(path)


def fig3(bfd):
    c = Counter()
    for r in bfd:
        oc = "unassigned" if r["clade"] in ("", "unassigned") else f"Clade {r['clade']}"
        c[(oc, old_tier(r["clade_confidence"]), r["place_our_clade"])] += 1
    return [{"old_clade": k[0], "old_tier": k[1], "new_clade": k[2], "n": v}
            for k, v in sorted(c.items())]


# ---------- Figure 4 ----------

def parse_matched_md(path, section="characterized, matched identity"):
    """Rows of the named section of MATCHED_IDENTITY.md (written by matched_identity.py)."""
    out, on = [], False
    for line in open(path):
        if line.startswith("### "):
            on = line[4:].startswith(section)
            continue
        m = re.match(r"\| (70-100|55-70|45-55|35-45|<35) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \|", line)
        if on and m:
            out.append({"bin": m.group(1), "tested": int(m.group(2)), "correct": int(m.group(3)),
                        "declined": int(m.group(4)), "wrong": int(m.group(5))})
    return out


def fig4(matched_md, skylark_tsv):
    rows = [dict(r, caller="187-leaf reference (Zamocky 2022)") for r in parse_matched_md(matched_md)]
    agg = defaultdict(Counter)
    for r in tsv_rows(skylark_tsv):
        # ceilings 70/55/45/40/35 match the 187-leaf design; 95 is near-plain leave-one-out
        if float(r["threshold"]) > 70 or r["outcome"].startswith("impossible"):
            continue
        b = bin_of(float(r["nearest_surviving_identity"]))
        agg[b]["tested"] += 1
        agg[b][{"correct": "correct", "declined": "declined"}.get(r["outcome"], "wrong")] += 1
    for b in BINS:
        if agg[b]["tested"]:
            rows.append({"caller": "19-reference (characterized only)", "bin": b,
                         "tested": agg[b]["tested"], "correct": agg[b]["correct"],
                         "declined": agg[b]["declined"], "wrong": agg[b]["wrong"]})
    return rows


# ---------- Table T1 ----------

def is_core(r):
    return (r.get("p39_Cys_proximal") == "C" and r.get("p109_E_of_EHD") == "E"
            and r.get("p111_D_of_EHD") == "D" and r.get("p180_Acid") == "E")


def t1_activity(s8_tsv):
    """Per clade: per-enzyme max peroxygenase conversion over pH, then min-max across enzymes."""
    by = defaultdict(list)
    for r in tsv_rows(s8_tsv):
        v = max(float(r[k]) for k in r if k.startswith("perox_"))
        by[f"Clade {r['clade']}"].append(v)
    return {c: f"{min(v):.0f}-{max(v):.0f} (n={len(v)})" for c, v in by.items()}


def table_t1(bfd, members, support, loo_tsv, s8_tsv):
    act = t1_activity(s8_tsv)
    sup = {(s["system"], s["clade"]): s for s in support}
    refs = defaultdict(lambda: [0, 0])          # (leaves, characterized)
    for r in members:
        for system, col in (("our", "our_clade"), ("fig1", "fig1_clade")):
            if r[col]:
                refs[(system, r[col])][0] += 1
                refs[(system, r[col])][1] += r["set"] == "backbone"
    loo = defaultdict(lambda: [0, 0])
    for r in tsv_rows(loo_tsv):
        loo[r["our_expected"]][1] += 1
        loo[r["our_expected"]][0] += r["our_ok"] == "yes"
    groups = defaultdict(list)
    for r in bfd:
        groups[("our", r["place_our_clade"])].append(r)
        groups[("fig1", r["place_fig1_clade"] or "unassigned")].append(r)
    rows = []
    order = [("our", c) for c in CLADES + ["unassigned"]] + sorted(
        k for k in groups if k[0] == "fig1" and k[1] != "unassigned") + [("fig1", "unassigned")]
    for key in order:
        rs = groups.get(key, [])
        if not rs:
            continue
        system, clade = key
        nref, nchar = refs.get(key, [0, 0])
        s = sup.get(key, {})
        lin = Counter(r.get("LINEAGE_GROUP") or "?" for r in rs).most_common(3)
        core = sum(is_core(r) for r in rs)
        rows.append({
            "system": "Clades I-V" if system == "our" else "Zamocky Figure 1",
            "clade": "outside Clades I-V" if (system, clade) == ("our", "unassigned")
                     else ("no Figure 1 clade" if clade == "unassigned" else clade),
            "ref_leaves": nref, "characterized": nchar,
            "support_sh_ufb": f"{s['sh_alrt']}/{s['ufboot']}" if s.get("sh_alrt") else "",
            "loo_correct": (f"{loo[clade][0]}/{loo[clade][1]}" if system == "our" else ""),
            "bfd_proteins": len(rs), "catalytic_core": core,
            "core_pct": f"{100 * core / len(rs):.1f}",
            "glu180_pct": f"{100 * sum(r.get('p180_Acid') == 'E' for r in rs) / len(rs):.1f}",
            "mean_length": f"{statistics.mean(int(r['length']) for r in rs):.0f}",
            "high_conf_pct": f"{100 * sum(float(r['place_our_lwr'] or 0) >= 0.95 for r in rs) / len(rs):.1f}"
                             if system == "our" and clade != "unassigned" else "",
            "top_lineages": "; ".join(f"{k} {v}" for k, v in lin),
            "t1_peroxygenase_pct": act.get(clade, "not measured") if system == "our" else "",
        })
    return rows


T1_COLS = [("clade", "clade"), ("ref_leaves", "ref. leaves"), ("characterized", "charact."),
           ("support_sh_ufb", "SH-aLRT/UFBoot"), ("loo_correct", "LOO"),
           ("bfd_proteins", "BFD proteins"), ("catalytic_core", "catalytic core"),
           ("glu180_pct", "Glu180 %"), ("mean_length", "mean aa"),
           ("high_conf_pct", "LWR>=0.95 %"), ("t1_peroxygenase_pct", "T1 peroxygenase %"),
           ("top_lineages", "top lineages")]


def t1_markdown(rows):
    out = ["# Table T1: UPO clades, references and BFD members", "",
           "Generated by `figures/bin/prep_figures.py`. Clade calls: placement on the 187-leaf "
           "reference (`upo_place`). T1 peroxygenase: per-enzyme maximum over pH on "
           "4-propylguaiacol (Barros-Reguera et al. 2025, Fig. S8), range across the clade's "
           "characterized enzymes. Catalytic core: Cys39, Glu109, Asp111, Glu180 (HspUPO "
           "numbering).", ""]
    for system in ("Clades I-V", "Zamocky Figure 1"):
        out += [f"## {system}", "", "| " + " | ".join(h for _, h in T1_COLS) + " |",
                "|" + "---|" * len(T1_COLS)]
        for r in rows:
            if r["system"] == system:
                out.append("| " + " | ".join(
                    f"{r[k]:,}" if isinstance(r[k], int) else str(r[k] or "")
                    for k, _ in T1_COLS) + " |")
        out.append("")
    return "\n".join(out) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", required=True, help="analysis/enzymes_UPO")
    ap.add_argument("--bfd", required=True, help="merged placement table (upo_place_bfd_all.tsv.zst)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    R = a.root
    tree = os.path.join(R, "data/refset_zamocky2022/tree")
    pref = os.path.join(R, "data/refset_zamocky2022/placement_ref")
    os.makedirs(a.out, exist_ok=True)
    members = tsv_rows(os.path.join(tree, "mfp.clade_members.tsv"))
    tips, sup = fig1(os.path.join(tree, "mfp.clade_members.tsv"),
                     os.path.join(tree, "mfp.treefile"), os.path.join(R, "refset/bin"))
    write(os.path.join(a.out, "fig1_tips.tsv"), tips)
    write(os.path.join(a.out, "fig1_support.tsv"), sup)
    write(os.path.join(a.out, "fig4_precision.tsv"),
          fig4(os.path.join(pref, "MATCHED_IDENTITY.md"),
               os.path.join(R, "data/placement/matched_identity_validation.tsv")),
          ["caller", "bin", "tested", "correct", "declined", "wrong"])
    bfd = read_bfd(a.bfd)
    write(os.path.join(a.out, "fig3_flows.tsv"), fig3(bfd))
    t1 = table_t1(bfd, members, sup, os.path.join(pref, "loo.tsv"),
                  os.path.join(R, "data/calib_dl/t1_peroxygenase_S8.tsv"))
    write(os.path.join(a.out, "table_T1.tsv"), t1)
    open(os.path.join(a.out, "table_T1.md"), "w").write(t1_markdown(t1))
    print(f"wrote figure data for {len(bfd)} BFD proteins to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
