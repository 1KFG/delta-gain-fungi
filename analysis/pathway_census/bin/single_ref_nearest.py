#!/usr/bin/env python3
"""Reading aid for a step 5 run of a single-reference step. Reads one step-5 output folder; writes <out-dir>/single_ref_report/.

Files:
  nearest_labelled.tsv     each representative: nearest labelled leaf (reference, named member, decoy, outgroup), patristic
                           distance, clade call, candidates in its mmseqs cluster, genus and class of its genome.
  reference_neighbours.tsv representatives sorted by patristic distance to the reference (all of them), with the same columns.
  labelled_distances.tsv   patristic distance from the reference to every other labelled leaf.
  target_by_genus.tsv      candidates called `target` per genus and class (genome of the protein).
  calls_by_class.tsv       per class: genomes, share with a candidate (floor of the run), share with an assigned call.
  summary.json             counts (representatives per call, candidates per call, genomes).
Distances are patristic (sum of branch lengths on the path) in the ML tree ref.raxml.bestTree, LG+G4 substitutions per site.
Nothing here is a call; it describes the tree and the placement.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

_HERE = Path(__file__).parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


sp = _load("step5_place")
hbc = _load("hits_by_class_floors")
pvh = hbc.pvh


def distances_from(root, source, others):
    """[(leaf, patristic distance from `source`)] for the leaves in `others`, sorted by distance then name."""
    paths = sp.leaf_depths(root)
    return sorted(((o, sp.tree_distance(paths, source, o)) for o in others if o != source), key=lambda x: (x[1], x[0]))


def nearest_among(root, queries, labelled):
    """[(query, nearest labelled leaf, distance)]; sp.nearest_labelled, kept here for one import point."""
    return sp.nearest_labelled(root, queries, labelled)


def genus_class(genome, tax):
    t = tax.get(genome, {})
    return t.get("GENUS") or "(no label)", t.get("CLASS") or "(no label)"


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def run(out_dir, ref_id, species, genomes_tsv):
    import duckdb
    d = Path(out_dir)
    rep_dir = d / "single_ref_report"
    rep_dir.mkdir(exist_ok=True)
    upo = sp.load_upo()
    top = upo.parse_newick((d / "ref.raxml.bestTree").read_text())
    cm = sp.read_tsv(d / "clade_members.tsv")
    role = {r["label"]: r["role"] for r in cm}
    clade = {r["label"]: r["our_clade"] or "(none)" for r in cm}
    labelled = [l for l, r in role.items() if r != "representative"]
    reps = [l for l, r in role.items() if r == "representative"]
    size = {rep: len(m) for rep, m in sp.parse_clusters(d / "clusters.tsv.zst").items()}
    tax = {r[0]: {"GENUS": r[1], "CLASS": r[2]} for r in duckdb.connect().execute(
        f"select LOCUSTAG, GENUS, \"CLASS\" from read_parquet('{species}')").fetchall()}

    def gc(leaf):
        return genus_class(leaf.split("_")[0], tax)
    cols = ["representative", "nearest_labelled", "nearest_role", "distance", "clade_call", "n_candidates_in_cluster", "genus", "class"]
    rows = []
    for q, n, dist in nearest_among(top, reps, labelled):
        g, c = gc(q)
        rows.append([q, n, role[n], f"{dist:.4f}", clade[q], size.get(q, ""), g, c])
    write_tsv(rep_dir / "nearest_labelled.tsv", cols, rows)
    cols2 = ["representative", "distance_to_reference", "clade_call", "n_candidates_in_cluster", "genus", "class"]
    rows2 = [[q, f"{dist:.4f}", clade[q], size.get(q, ""), *gc(q)] for q, dist in distances_from(top, ref_id, reps)]
    write_tsv(rep_dir / "reference_neighbours.tsv", cols2, rows2)
    ld = [[l, role[l], f"{dist:.4f}"] for l, dist in distances_from(top, ref_id, labelled)]
    write_tsv(rep_dir / "labelled_distances.tsv", ["leaf", "role", "distance_to_reference"], ld)
    det = sp.read_tsv(d / "placements_detail.tsv")
    cand = sp.read_tsv(d / "candidates_step.tsv")
    tg = Counter(gc(r["protein_id"]) for r in det if r["clade_call"] == "target")
    write_tsv(rep_dir / "target_by_genus.tsv", ["genus", "class", "n_candidates_target"],
              [[g, c, n] for (g, c), n in sorted(tg.items(), key=lambda x: (-x[1], x[0]))])
    asg = sp.read_tsv(next(d.glob("assigned*.tsv")))
    hit_g = {r["protein_id"].split("_")[0] for r in cand}
    asg_g = {r["protein_id"].split("_")[0] for r in asg if r["label"] == "assigned"}
    genomes = [r["locustag"] for r in pvh.read_tsv(genomes_tsv) if r["has_proteins"].strip().lower() == "true"]
    classes = {k: v["CLASS"] for k, v in tax.items()}
    table = hbc.class_table(genomes, classes, {"candidate": hit_g, "assigned": asg_g}, 20)
    gs = set(genomes)
    with open(rep_dir / "calls_by_class.tsv", "w") as fh:
        fh.write("class\tn_genomes\tn_with_candidate\tn_with_assigned\tshare_candidate\tshare_assigned\n")
        for c, n, s in sorted(table, key=lambda r: -r[1]):
            fh.write(f"{c}\t{n}\t{round(s['candidate'] * n)}\t{round(s['assigned'] * n)}\t{s['candidate']:.4f}\t{s['assigned']:.4f}\n")
        fh.write(f"ALL\t{len(genomes)}\t{len(hit_g & gs)}\t{len(asg_g & gs)}\t{len(hit_g & gs) / len(genomes):.4f}\t{len(asg_g & gs) / len(genomes):.4f}\n")
    summ = {"representatives_by_call": dict(Counter(clade[q] for q in reps)),
            "candidates_by_call": dict(Counter(r["clade_call"] for r in det)),
            "genomes_with_candidate": len(hit_g & gs), "genomes_with_assigned": len(asg_g & gs), "genomes_searched": len(genomes)}
    (rep_dir / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ref-id", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--genomes", required=True)
    a = ap.parse_args(argv)
    print(json.dumps(run(a.out_dir, a.ref_id, a.species, a.genomes)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
