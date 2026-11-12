#!/usr/bin/env python3
"""How far does each candidate effector reach across Mucoromycota?

  conservation.py --clusters clu_cluster.tsv --genomes genomes.tsv --out conservation.tsv \
      [--focus JEM1MUCOR,JEM2MUCOR] [--extra-genomes extra.tsv]

Lineage restriction is the sequence-independent effector signal: a secreted protein found only
in one strain is a different proposition from one conserved across the phylum, whatever its
sequence says. Together with the fold calls this gives the two axes that matter - "no known
fold" and "found nowhere else" - without either standing alone.

Input is the MMseqs2 easy-cluster TSV (representative<TAB>member), clustered over the candidate
secretomes of every Mucoromycota genome in BFD plus the genomes being screened, and the genome
table that maps LOCUSTAG to species, genus and order. Classes, widest first:

  mucoromycota_wide  members in more than one order
  mucorales_wide     members in more than one genus of one order
  genus_restricted   members in more than one species of one genus
  species_restricted members in more than one genome of one species
  genome_specific    every member is in the one genome

`genome_specific` is the strongest effector-like signal here and also the one most easily
produced by a gene-prediction artefact, so it is reported with the cluster size rather than
alone. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict


def read_genomes(path, *paths):
    """genomes.tsv -> {LOCUSTAG: (species, genus, order, asmid)}; later files win on conflict."""
    out = {}
    for p in (path,) + paths:
        if not p:
            continue
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                lt = r.get("LOCUSTAG") or r.get("locustag")
                if not lt:
                    continue
                out[lt] = (r.get("SPECIES", ""), r.get("GENUS", ""),
                           r.get("ORD") or r.get("ORDER", ""), r.get("ASMID", ""))
    return out


def read_clusters(handle):
    """MMseqs easy-cluster TSV -> {representative: [members]}."""
    clu = defaultdict(list)
    for line in handle:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] and f[1]:
            clu[f[0]].append(f[1])
    return clu


def classify(members, genomes):
    """-> (class, n_members, n_genomes, n_species, n_genera, n_orders)."""
    lts = [m.split("_")[0] for m in members]
    known = [genomes[lt] for lt in lts if lt in genomes]
    # Only genomes present in the genome table count towards reach. A member whose LOCUSTAG is
    # unknown contributes no lineage, so it must not contribute a genome either - otherwise a
    # single unrecognised member silently promotes a genome_specific family to
    # species_restricted, which is the direction that loses effector candidates.
    n_genomes = len({lt for lt in lts if lt in genomes})
    species = {i[0] for i in known if i[0]}
    genera = {i[1] for i in known if i[1]}
    orders = {i[2] for i in known if i[2]}
    if len(orders) > 1:
        cls = "mucoromycota_wide"
    elif len(genera) > 1:
        cls = "mucorales_wide"
    elif len(species) > 1:
        cls = "genus_restricted"
    elif n_genomes > 1:
        cls = "species_restricted"
    else:
        cls = "genome_specific"
    return cls, len(members), n_genomes, len(species), len(genera), len(orders)


def build(clusters, genomes, focus=None):
    """-> {protein_id: row} for every protein of the focus locustags (all, if focus is None)."""
    want = set(focus) if focus else None
    rows = {}
    for rep, members in clusters.items():
        allm = members if rep in members else [rep] + members
        cls, n, ng, nsp, ngen, nord = classify(allm, genomes)
        for m in allm:
            lt = m.split("_")[0]
            if want and lt not in want:
                continue
            rows[m] = {"protein_id": m, "LOCUSTAG": lt, "conservation": cls,
                       "cluster": rep, "cluster_members": n, "cluster_genomes": ng,
                       "cluster_species": nsp, "cluster_genera": ngen, "cluster_orders": nord}
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--extra-genomes", help="genome table for genomes outside BFD")
    ap.add_argument("--focus", default="", help="comma-separated LOCUSTAGs to report")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    genomes = read_genomes(a.genomes, a.extra_genomes)
    with open(a.clusters) as fh:
        clusters = read_clusters(fh)
    focus = [x for x in a.focus.split(",") if x]
    rows = build(clusters, genomes, focus or None)

    cols = ["protein_id", "LOCUSTAG", "conservation", "cluster_members", "cluster_genomes",
            "cluster_species", "cluster_genera", "cluster_orders", "cluster"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for p in sorted(rows):
            w.writerow(rows[p])
    c = Counter(r["conservation"] for r in rows.values())
    per = defaultdict(Counter)
    for r in rows.values():
        per[r["LOCUSTAG"]][r["conservation"]] += 1
    print("conservation: " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())), file=sys.stderr)
    for lt in sorted(per):
        print(f"  {lt}: " + ", ".join(f"{k}={v}" for k, v in sorted(per[lt].items())),
              file=sys.stderr)
    print(f"{len(rows)} proteins -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
