#!/usr/bin/env python3
"""Host target panel for the liverworts, built identically for both species.

  plant_panel.py --dir out/panels --locustags MPTAK,MPALE --out out/plant_panel

The host species of the JEM isolates is not recorded, so both Marchantia species are annotated
and panelled: M. polymorpha (MPTAK, GCA_039105155.1) and M. paleacea (MPALE, GCA_014180765.2).
They differ in the way that matters here - M. paleacea forms arbuscular mycorrhizal symbiosis
and M. polymorpha has lost it - so a class present in one and absent in the other is itself the
interesting result, and a class present in both makes the unknown host moot.

Classes are defined on Pfam accessions (see CLASSES). Receptor kinases enter the panel as their
**ectodomain**, not the full receptor: the interaction surface is outside the membrane, and a
cytoplasmic kinase domain is both the largest part of the chain and irrelevant to an effector
contact, so including it would dominate any structure prediction. The ectodomain is taken as
the mature N-terminus (after the signal-peptide cleavage site) up to the residue before the
first TMHMM helix. Proteins with a signal peptide and no helix are taken as the whole mature
chain; proteins with neither are skipped.

What this does NOT do is assert orthology between the two species' proteins. It counts and
constructs by class; which MPTAK protein corresponds to which MPALE protein is a tree question,
not a Pfam question. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_secretome as bs  # noqa: E402
from extract_proteins import read_fasta  # noqa: E402

# Pfam accessions per class. A protein takes the first class it matches, so the more specific
# classes come first.
CLASSES = [
    ("lysm_receptor", {"PF01476"},
     "LysM domain - the chitin/peptidoglycan perception module; MpLYK1/MpLYR in M. polymorpha, "
     "MpaLYR/MpaCERK1 in M. paleacea"),
    ("chitinase", {"PF00704", "PF00182", "PF03067"},
     "GH18/GH19 chitinase - releases chitin oligomers from a fungal wall"),
    ("chitin_binding", {"PF00187", "PF01607"},
     "CBM18/CBM14 chitin-binding module"),
    ("lrr_receptor", {"PF00560", "PF13855", "PF07714", "PF08263", "PF13516"},
     "LRR receptor kinase - the SERK/BAK1-type co-receptor class"),
    ("cys_protease", {"PF00112"},
     "papain-like cysteine protease - the Rcr3/Pip1 class targeted by fungal inhibitors"),
    ("protease_inhibitor", {"PF00031", "PF00197", "PF00280", "PF00079"},
     "cystatin / Kunitz / potato inhibitor I"),
    ("subtilase", {"PF00082"}, "subtilisin-like serine protease"),
    ("pr_protein", {"PF00188", "PF00407", "PF00314", "PF00967"},
     "PR-1 (CAP), Bet v 1, thaumatin, barwin"),
    ("germin", {"PF01072"}, "germin / oxalate oxidase"),
    ("lectin", {"PF01419", "PF00139"}, "jacalin-related or legume lectin"),
    ("expansin", {"PF03330", "PF01357"}, "expansin - cell-wall loosening"),
]
CLASS_OF = {acc: name for name, accs, _ in CLASSES for acc in accs}
DESC_OF = {name: d for name, _, d in CLASSES}


def classify_protein(pfam_accs):
    """-> class name, by the first entry of CLASSES that matches, or '' for none."""
    accs = set(pfam_accs)
    for name, want, _ in CLASSES:
        if accs & want:
            return name
    return ""


def construct_range(length, cleavage, topo, n_helices):
    """-> (start, end, kind) 1-based inclusive, the part of the chain to model.

    ectodomain  signal peptide and a membrane helix: mature N-terminus to just before the helix
    mature      signal peptide, no helix: the whole secreted chain
    none        no signal peptide: not in the panel (an intracellular protein is not a surface
                target for a secreted effector)"""
    start = (cleavage or 0) + 1
    if start <= 1:
        return (0, 0, "none")
    helices = bs.parse_topology(topo)
    first = None
    for s, _e in helices:
        if s > (cleavage or 0) + bs.TM_SLACK:
            first = s
            break
    if first is None and bs.fnum(n_helices) > 0 and not helices:
        return (0, 0, "none")          # helices present but no topology to place them
    if first is not None:
        end = first - 1
        if end - start + 1 < 30:
            return (0, 0, "none")      # nothing left outside the membrane worth modelling
        return (start, end, "ectodomain")
    return (start, length, "mature")


def build(tables, seqs, locustags):
    sig, tm, gpi, tp, pfam, cazy, subs, length = tables
    want = set(locustags)
    rows = []
    for p, accs in pfam.items():
        lt = p.split("_")[0]
        if lt not in want:
            continue
        cls = classify_protein(a for a, _ in accs)
        if not cls:
            continue
        s = sig.get(p, {})
        t = tm.get(p, {})
        cleave = int(bs.fnum(s.get("peptide_end"))) if s else 0
        slen = length.get(p, 0) or len(seqs.get(p, "").rstrip("*"))
        start, end, kind = construct_range(slen, cleave, t.get("Topology"), t.get("PredHel"))
        if kind == "none":
            continue
        seq = seqs.get(p, "").rstrip("*")
        rows.append({
            "protein_id": p, "species": lt, "class": cls, "construct": kind,
            "start": start, "end": end, "construct_len": end - start + 1,
            "full_len": slen, "n_helices": t.get("PredHel", "0"),
            "pfam": ";".join(sorted({a for a, _ in accs})),
            "seq": seq[start - 1:end] if seq else "",
        })
    return sorted(rows, key=lambda r: (r["class"], r["species"], r["protein_id"]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--locustags", required=True)
    ap.add_argument("--fasta", default="", help="proteome FASTA(s), comma separated")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    lts = [x for x in a.locustags.split(",") if x]

    import call_effectors as ce
    seqs = {}
    for f in [x for x in a.fasta.split(",") if x]:
        seqs.update(dict(read_fasta(f)))
    if not seqs:
        for lt in lts:
            f = os.path.join(a.dir, "faa", f"{lt}.faa")
            if os.path.exists(f):
                seqs.update(dict(read_fasta(f)))

    rows = build(ce.load_tables(a.dir), seqs, lts)
    cols = ["protein_id", "species", "class", "construct", "start", "end", "construct_len",
            "full_len", "n_helices", "pfam"]
    with open(os.path.join(a.out, "panel.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(a.out, "panel.faa"), "w") as fh:
        for r in rows:
            if r["seq"]:
                fh.write(f">{r['protein_id']}|{r['class']}|{r['construct']}\n{r['seq']}\n")

    # side by side, which is the point of annotating both species
    per = defaultdict(Counter)
    for r in rows:
        per[r["class"]][r["species"]] += 1
    with open(os.path.join(a.out, "compare.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["class"] + lts + ["description"])
        for name, _accs, desc in CLASSES:
            if name in per:
                w.writerow([name] + [per[name].get(lt, 0) for lt in lts] + [desc])
    print(f"{len(rows)} panel constructs -> {a.out}", file=sys.stderr)
    for name, _a, _d in CLASSES:
        if name in per:
            print("  " + name + ": " + ", ".join(f"{lt}={per[name].get(lt, 0)}" for lt in lts),
                  file=sys.stderr)


if __name__ == "__main__":
    main()
