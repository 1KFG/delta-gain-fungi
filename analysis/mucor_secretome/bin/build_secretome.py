#!/usr/bin/env python3
"""Mucoromycota secretome, and a tiered screen set for one focal genome (Mucor hiemalis).

  build_secretome.py --dir <sql export dir> --codons codon_check.tsv --clusters clusters.tsv \
      --fasta candidates.faa --focus F00CDB58 --out out/

Built to supply the input list for an AlphaFold3-style screen of fungal secreted proteins against
plant candidate genes. Rules are fixed here, before the run, so they can be read off the code.

Secretion call (from BFD's own predictions):
  signal peptide        SignalP probability >= 0.7
  membrane              a TMHMM helix starting more than 5 aa after the signal peptide cleavage
                        site (helices inside the signal peptide are the signal peptide itself)
  GPI-anchored          a PredGPI feature; kept, but called `gpi_anchored`, not `secreted`
  mitochondrial         TargetP mTP at probability >= 0.7; excluded
  `secreted` means: signal peptide, no membrane helix, no GPI anchor, not mitochondrial.

5' completeness: taken from the sequence (bin/check_codons.py), never from BFD's
`has_start_codon` flag, which is wrong in these genomes before Fungi_BFD fb1d7cd
(enzymes_UPO/bfd_findings/REPORT.md). A model without a real ATG start may have lost its signal
peptide, so it cannot enter the screen set.

Tiers, highest first. A protein takes the first tier it qualifies for:
  T1 chitin/chitosan   CE4, GH18, GH20, GH46, GH75 or a chitin-binding module (CBM12/14/18/50),
                       or Pfam PF01522 (polysaccharide deacetylase) / PF01476 (LysM).
                       The Mucorales cell wall is chitosan-rich and plant LysM receptors bind
                       chitosan poorly, so this is where perception and masking would act.
  T2 plant-facing      a CAZyme whose dbCAN substrate names a plant polysaccharide, or a cutinase,
                       lipase or peroxygenase/peroxidase domain: enzymes acting on plant material.
  T3 effector-like     no Pfam domain, <= 300 aa, >= 4 cysteines and >= 2 % cysteine.
  T4 other secreted    everything else.

Conservation comes from MMseqs2 clusters over the whole Mucoromycota candidate secretome: for each
focal protein, how many species, genera and orders hold a cluster member.
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_proteins import read_fasta  # noqa: E402

SIGP_MIN = 0.7
MITO_MIN = 0.7
TM_SLACK = 5          # aa after the cleavage site within which a helix is the signal peptide
EFF_MAX_LEN = 300
EFF_MIN_CYS = 4
EFF_MIN_CYS_FRAC = 0.02

CHITIN_CAZY = ("CE4", "GH18", "GH20", "GH46", "GH75", "CBM12", "CBM14", "CBM18", "CBM50")
CHITIN_PFAM = {"PF01522": "Polysacc_deac_1", "PF01476": "LysM"}
PLANT_WORDS = ("cellulose", "xylan", "pectin", "xyloglucan", "mannan", "arabinan", "galactan",
               "lignin", "cutin", "starch", "arabinoxylan", "rhamnogalacturonan", "glucuronan",
               "pectate", "cellobiose", "xylose", "polygalacturonan")
PLANT_PFAM = {"PF01083": "Cutinase", "PF01328": "Peroxidase_2 (UPO)", "PF00141": "peroxidase",
              "PF01764": "Lipase_3", "PF00657": "SGNH lipase/esterase", "PF03583": "LIP",
              "PF07859": "Abhydrolase_3"}
TIER_LABEL = {"T1": "chitin/chitosan", "T2": "plant-facing enzyme", "T3": "effector-like unknown",
              "T4": "other secreted"}


def tsv(path, optional=False):
    if optional and not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gpi_set(rows):
    """Protein IDs with a GPI anchor. PredGPI emits a `Chain` row for every protein and a
    `GPI-anchor` row only for the anchored ones, so an unfiltered table would mark everything."""
    return {r["protein_id"] for r in rows
            if r.get("feature", "GPI-anchor").strip() == "GPI-anchor"}


def fnum(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def parse_topology(topo):
    """TMHMM topology string ('i7-29o44-66i') -> [(start, end)] of the predicted helices."""
    return [(int(a), int(b)) for a, b in re.findall(r"(\d+)-(\d+)", topo or "")]


def membrane_helices(topo, n_helices, cleavage):
    """Helices that start more than TM_SLACK residues after the signal-peptide cleavage site.
    Falls back to the helix count when the topology string is missing."""
    hel = parse_topology(topo)
    if not hel:
        return max(0, int(n_helices or 0) - (1 if cleavage else 0))
    return sum(1 for s, _ in hel if s > (cleavage or 0) + TM_SLACK)


def classify(sig_p, topo, n_helices, cleavage, gpi, mito_p):
    """-> one of secreted / gpi_anchored / membrane / mitochondrial / weak_signal."""
    if sig_p < SIGP_MIN:
        return "weak_signal"
    if membrane_helices(topo, n_helices, cleavage) > 0:
        return "membrane"
    if gpi:
        return "gpi_anchored"
    if mito_p >= MITO_MIN:
        return "mitochondrial"
    return "secreted"


def cys_stats(seq):
    n = len(seq.rstrip("*"))
    c = seq.upper().count("C")
    return c, (c / n if n else 0.0), n


def is_effector_like(length, n_cys, cys_frac, has_pfam):
    return (not has_pfam and length <= EFF_MAX_LEN and n_cys >= EFF_MIN_CYS
            and cys_frac >= EFF_MIN_CYS_FRAC)


def cazy_root(fam):
    """'CE4(145-269)' or 'GH18_3' -> 'CE4' / 'GH18'."""
    return re.split(r"[(_]", (fam or "").strip())[0]


def tier_of(cazy_fams, substrates, pfam_accs, length, n_cys, cys_frac):
    roots = {cazy_root(f) for f in cazy_fams}
    hits = sorted(roots & set(CHITIN_CAZY)) + sorted(set(pfam_accs) & set(CHITIN_PFAM))
    if hits:
        return "T1", ",".join(hits)
    sub = " ".join(substrates).lower()
    words = sorted({w for w in PLANT_WORDS if w in sub})
    pf = sorted(set(pfam_accs) & set(PLANT_PFAM))
    if words or pf:
        return "T2", ",".join(words + [PLANT_PFAM[p] for p in pf])
    if is_effector_like(length, n_cys, cys_frac, bool(pfam_accs)):
        return "T3", f"{n_cys} Cys, {100 * cys_frac:.1f} %, {length} aa"
    return "T4", ""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True, help="directory of the SQL export")
    ap.add_argument("--codons", required=True)
    ap.add_argument("--clusters", required=True, help="MMseqs2 cluster TSV: representative, member")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--focus", required=True, help="LOCUSTAG of the focal genome")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    D = lambda f: os.path.join(a.dir, f)  # noqa: E731

    genomes = {r["LOCUSTAG"]: r for r in tsv(D("genomes.tsv"))}
    focus = genomes.get(a.focus)
    if focus is None:
        sys.exit(f"focal LOCUSTAG {a.focus} is not in {D('genomes.tsv')}")

    sig = {r["protein_id"]: r for r in tsv(D("signalp.tsv"))}
    tm = {r["protein_id"]: r for r in tsv(D("tmhmm.tsv"))}
    gpi = gpi_set(tsv(D("gpi.tsv"), optional=True))
    tp = {r["protein_id"]: r for r in tsv(D("targetp.tsv"), optional=True)}
    length = {r["protein_id"]: int(r["length"]) for r in tsv(D("lengths.tsv"), optional=True)}
    pfam, cazy, subs = defaultdict(list), defaultdict(list), defaultdict(list)
    for r in tsv(D("pfam.tsv")):
        pfam[r["protein_id"]].append((r["pfam_acc"], r["pfam_id"]))
    for r in tsv(D("cazy.tsv"), optional=True):
        cazy[r["protein_id"]].append(r["cazyme_fam"])
        if r.get("substrate"):
            subs[r["protein_id"]].append(r["substrate"])
    codon = {r["transcript_id"]: r for r in tsv(a.codons)}
    seqs = dict(read_fasta(a.fasta))

    cluster_of, members = {}, defaultdict(list)
    with open(a.clusters) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                cluster_of[f[1]] = f[0]
                members[f[0]].append(f[1])

    # --- classify every candidate in every Mucoromycota genome
    rows = {}
    for p, s in sig.items():
        lt = p.split("_")[0]
        if lt not in genomes:
            continue
        t = tm.get(p, {})
        cleave = int(fnum(s.get("peptide_end")))
        cat = classify(fnum(s.get("probability")), t.get("Topology"), t.get("PredHel"), cleave,
                       p in gpi, fnum(tp.get(p, {}).get("probability"))
                       if (tp.get(p, {}).get("prediction") or "").lower().startswith("mt") else 0.0)
        rows[p] = {"protein_id": p, "LOCUSTAG": lt, "category": cat,
                   "sigp_prob": s.get("probability"), "sigp_cleavage": cleave}

    secreted = {p for p, r in rows.items() if r["category"] in ("secreted", "gpi_anchored")}

    # --- conservation of each cluster, over the whole Mucoromycota candidate secretome
    def spec(lt, key):
        return (genomes.get(lt) or {}).get(key, "")
    cons = {}
    for rep, mem in members.items():
        lts = {m.split("_")[0] for m in mem if m in secreted}
        cons[rep] = {"n_members": len(mem),
                     "n_genomes": len(lts),
                     "species": sorted({spec(l, "SPECIES") for l in lts} - {""}),
                     "genera": sorted({spec(l, "GENUS") for l in lts} - {""}),
                     "orders": sorted({spec(l, "ORD") for l in lts} - {""})}

    # --- focal genome table
    out_rows = []
    for p, r in sorted(rows.items()):
        if r["LOCUSTAG"] != a.focus:
            continue
        seq = seqs.get(p, "")
        n_cys, frac, slen = cys_stats(seq)
        if not slen:
            slen = length.get(p, 0)
        pacc = [x[0] for x in pfam.get(p, [])]
        tier, why = tier_of(cazy.get(p, []), subs.get(p, []), pacc, slen, n_cys, frac)
        cl = cluster_of.get(p, p)
        c = cons.get(cl, {"n_members": 1, "n_genomes": 1, "species": [], "genera": [], "orders": []})
        cc = codon.get(p, {})
        five_ok = cc.get("seq_start") == "True"
        out_rows.append({
            "protein_id": p, "category": r["category"], "tier": tier,
            "tier_label": TIER_LABEL[tier], "tier_evidence": why,
            "length": slen, "n_cys": n_cys, "cys_pct": round(100 * frac, 2),
            "sigp_prob": r["sigp_prob"], "sigp_cleavage": r["sigp_cleavage"],
            "gpi_anchored": "yes" if p in gpi else "no",
            "tm_helices_after_sp": membrane_helices(tm.get(p, {}).get("Topology"),
                                                    tm.get(p, {}).get("PredHel"),
                                                    r["sigp_cleavage"]),
            "pfam": ";".join(f"{x[0]}:{x[1]}" for x in pfam.get(p, [])),
            "cazy": ";".join(sorted({cazy_root(f) for f in cazy.get(p, [])})),
            "cazy_substrate": ";".join(sorted(set(subs.get(p, [])))),
            "seq_start_ATG": "yes" if five_ok else ("no" if cc else "not_checked"),
            "seq_has_stop": "yes" if cc.get("seq_stop") == "True" else ("no" if cc else "not_checked"),
            "bfd_flag_has_start": cc.get("has_start_codon", ""),
            "cluster": cl, "cluster_members": c["n_members"], "cluster_genomes": c["n_genomes"],
            "cluster_species": len(c["species"]), "cluster_genera": len(c["genera"]),
            "cluster_orders": ";".join(c["orders"]),
            "conservation": ("genome_specific" if c["n_genomes"] <= 1 else
                             "genus_restricted" if len(c["genera"]) <= 1 else
                             "mucorales_wide" if c["orders"] == ["Mucorales"] else
                             "mucoromycota_wide"),
            "sequence": seq,
        })

    cols = [c for c in out_rows[0] if c != "sequence"] if out_rows else ["protein_id"]
    with open(os.path.join(a.out, "focal_secretome.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    # --- screen set: secreted (incl. GPI), real ATG start, tier order then conservation
    screen = [r for r in out_rows
              if r["category"] in ("secreted", "gpi_anchored") and r["seq_start_ATG"] == "yes"]
    screen.sort(key=lambda r: (r["tier"], -r["cluster_species"], r["protein_id"]))
    with open(os.path.join(a.out, "screen_set.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rank"] + cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for i, r in enumerate(screen, 1):
            w.writerow({"rank": i, **r})
    with open(os.path.join(a.out, "screen_set.faa"), "w") as fo:
        for r in screen:
            fo.write(f">{r['protein_id']} tier={r['tier']} {r['tier_label']}\n{r['sequence']}\n")

    # --- family presence/absence across Mucoromycota, for the focal clusters
    focal_clusters = {r["cluster"] for r in screen}
    species = sorted({g["SPECIES"] for g in genomes.values() if g["SPECIES"]})
    with open(os.path.join(a.out, "family_presence_absence.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cluster", "focal_proteins", "n_species"] + species)
        for cl in sorted(focal_clusters):
            have = {spec(m.split("_")[0], "SPECIES") for m in members.get(cl, [cl])
                    if m in secreted}
            foc = [r["protein_id"] for r in screen if r["cluster"] == cl]
            w.writerow([cl, ";".join(foc), len(have - {""})]
                       + ["1" if s in have else "0" for s in species])

    # --- summary
    cat = Counter(r["category"] for r in out_rows)
    tiers = Counter(r["tier"] for r in screen)
    consv = Counter(r["conservation"] for r in screen)
    checked = [r for r in out_rows if r["seq_start_ATG"] in ("yes", "no")]
    bad_flag = sum(1 for r in checked
                   if r["seq_start_ATG"] == "yes" and r["bfd_flag_has_start"] == "FALSE")
    L = [r["length"] for r in screen]
    md = [f"# {focus['SPECIES']} secretome and screen set",
          "",
          f"Focal genome `{focus['ASMID']}` (LOCUSTAG `{a.focus}`). Candidate pool: "
          f"{len(genomes)} Mucoromycota genomes, {len(sig)} proteins with a signal peptide at "
          f"p >= 0.5.",
          "",
          "## Focal genome",
          "",
          "| category | proteins |", "|---|---|"]
    md += [f"| {k} | {v} |" for k, v in cat.most_common()]
    md += ["", f"Sequence-level start check: {len(checked)} focal proteins checked, "
               f"{sum(1 for r in checked if r['seq_start_ATG'] == 'yes')} begin with ATG. "
               f"{bad_flag} of them are flagged `has_start_codon = FALSE` by BFD "
               f"(the flag bug; see enzymes_UPO/bfd_findings/REPORT.md).",
           "", f"## Screen set: {len(screen)} proteins", "",
           "| tier | what it is | proteins |", "|---|---|---|"]
    md += [f"| {t} | {TIER_LABEL[t]} | {tiers[t]} |" for t in ("T1", "T2", "T3", "T4") if tiers[t]]
    md += ["", "| conservation across Mucoromycota | proteins |", "|---|---|"]
    md += [f"| {k} | {v} |" for k, v in consv.most_common()]
    if L:
        L = sorted(L)
        md += ["", f"Length: min {L[0]}, median {L[len(L) // 2]}, max {L[-1]} aa."]
    md += ["", "Files: `focal_secretome.tsv` (every candidate, classified), `screen_set.tsv` and "
               "`screen_set.faa` (the proteins to fold, ranked), `family_presence_absence.tsv` "
               "(each screen-set cluster across Mucoromycota species)."]
    with open(os.path.join(a.out, "SUMMARY.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"focal {len(out_rows)} candidates, screen set {len(screen)} -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
