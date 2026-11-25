#!/usr/bin/env python3
"""Reading aid for one step 5 run folder (any number of references). Writes <run>/report/.

Reads clade_report.json, clade_members.tsv, ref.raxml.bestTree, clusters.tsv.zst, placements_detail.tsv, candidates_step.tsv and
assigned*.tsv in the run folder, and the genome table (results/genomes.tsv.zst: locustag, phylum, class_, genus, has_proteins).
Files:
  summary.json          leaves by role; target size and references inside; calls (target / each other call) with support;
                        assigned proteins, genomes, copies per genome; check that candidates_step rows = placements rows.
  assigned_by_class.tsv per class: genomes searched, genomes with an assigned protein, share (classes with fewer than --min-n
                        genomes are merged); also by phylum.
  nearest_labelled.tsv  per representative: nearest labelled leaf (reference, named member, decoy, outgroup), patristic distance,
                        clade call, candidates in its mmseqs cluster.
  nearest_summary.tsv   per role of the nearest labelled leaf: representatives, candidates in their clusters, median distance.
Nothing here is a call. Distances are patristic distances in the ML tree (LG+G4 substitutions per site).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("step5_place", Path(__file__).parent / "step5_place.py")
sp = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sp)


def genome_of(protein_id: str) -> str:
    return protein_id.split("_")[0]


def class_table(searched, label_of, assigned, min_n=20):
    """searched: list of genome ids; label_of: genome -> label; assigned: set of genome ids with an assigned protein.
    Labels with fewer than min_n genomes are merged into '(other, small groups)'. Rows (label, n, n_assigned, share), by n."""
    by = defaultdict(list)
    for g in searched:
        by[label_of.get(g) or "(no label)"].append(g)
    big = {k: v for k, v in by.items() if len(v) >= min_n}
    small = [g for k, v in by.items() if len(v) < min_n for g in v]
    if small:
        big["(other, small groups)"] = small
    rows = [(k, len(v), sum(g in assigned for g in v), sum(g in assigned for g in v) / len(v)) for k, v in big.items()]
    return sorted(rows, key=lambda r: (-r[1], r[0]))


def nearest_summary(near_rows, cluster_size):
    """near_rows: [(rep, nearest, distance, role)] -> rows (role, n_reps, n_candidates, median distance)."""
    by = defaultdict(list)
    for rep, _, dist, role in near_rows:
        by[role].append((dist, cluster_size.get(rep, 1)))
    return [(role, len(v), sum(c for _, c in v), round(statistics.median(d for d, _ in v), 4)) for role, v in sorted(by.items())]


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def run(run_dir, genomes_tsv, min_n=20):
    d = Path(run_dir)
    out = d / "report"
    out.mkdir(exist_ok=True)
    cm = sp.read_tsv(d / "clade_members.tsv")
    role = {r["label"]: r["role"] for r in cm}
    call = {r["label"]: r["our_clade"] or "(none)" for r in cm}
    labelled = [l for l, r in role.items() if r not in ("representative", "other")]
    reps = [l for l, r in role.items() if r == "representative"]
    rep_json = json.loads((d / "clade_report.json").read_text())
    top = sp.load_upo().parse_newick((d / "ref.raxml.bestTree").read_text())
    size = {r: len(m) for r, m in sp.parse_clusters(d / "clusters.tsv.zst").items()}
    near = [(q, n, dist, role[n]) for q, n, dist in sp.nearest_labelled(top, reps, labelled)] if labelled and reps else []
    write_tsv(out / "nearest_labelled.tsv", ["representative", "nearest_labelled", "nearest_role", "distance", "clade_call", "n_candidates_in_cluster"],
              [[q, n, r, f"{dist:.4f}", call[q], size.get(q, "")] for q, n, dist, r in near])
    write_tsv(out / "nearest_summary.tsv", ["nearest_role", "n_representatives", "n_candidates_in_their_clusters", "median_distance"],
              nearest_summary(near, size))
    det = sp.read_tsv(d / "placements_detail.tsv")
    cand = sp.read_tsv(d / "candidates_step.tsv")
    asg = [r for r in sp.read_tsv(next(d.glob("assigned*.tsv"))) if r["label"] == "assigned"]
    asg_genomes = {genome_of(r["protein_id"]) for r in asg}
    gt = sp.read_tsv(genomes_tsv)
    searched = [r["locustag"] for r in gt if r["has_proteins"].strip().lower() == "true"]
    cls = {r["locustag"]: r["class_"] for r in gt}
    phy = {r["locustag"]: r["phylum"] for r in gt}
    rows = [("class",) + r for r in class_table(searched, cls, asg_genomes, min_n)]
    rows += [("phylum",) + r for r in class_table(searched, phy, asg_genomes, 1)]
    write_tsv(out / "assigned_by_class.tsv", ["level", "group", "n_genomes", "n_with_assigned", "share"],
              [[a, b, c, e, f"{f:.4f}"] for a, b, c, e, f in rows])
    calls = Counter(r["clade_call"] or "(none)" for r in det)
    sup = sorted(float(r["support"]) for r in det if r["clade_call"] == "target")
    per = Counter(genome_of(r["protein_id"]) for r in asg)
    summ = {"step_id": rep_json.get("step_id"), "n_candidates_step_rows": len(cand), "n_placements": len(det),
            "candidates_equal_placements": len(cand) == len(det),
            "leaves_by_role": dict(Counter(role.values())),
            "refs_form_clade_free_of_labelled_leaves": rep_json.get("refs_form_clade_free_of_decoys"),
            "target_size": rep_json.get("target_size"), "smallest_rule_target_size": rep_json.get("smallest_target_size"),
            "target_refs": rep_json.get("target_refs"), "target_representatives": rep_json.get("target_representatives"),
            "n_representatives": len(reps),
            "calls": dict(calls), "support_target_min": sup[0] if sup else None, "support_target_median": sup[len(sup) // 2] if sup else None,
            "assigned_proteins": len(asg), "genomes_with_assigned": len(asg_genomes), "genomes_searched": len(searched),
            "copies_per_genome": dict(sorted(Counter(per.values()).items())),
            "long_pendant": sum(r["long_pendant"] == "True" for r in det)}
    (out / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--min-n", type=int, default=20)
    a = ap.parse_args(argv)
    print(json.dumps(run(a.run, a.genomes, a.min_n)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
