#!/usr/bin/env python3.12
"""Bootstrap support and clade-rule stability for the joint PKS tree (step 5). Standard library plus the step 5 tools.

Subcommands:
  support    Support values (raxml-ng `--support`: FBP and TBE trees on the best ML tree) at the nodes that define the
             target clades: the target side of each step, the node of all references, and the nodes of the heptaketide
             references with the decoy wA. Writes a TSV (one row per node) and a JSON.
  stability  For every bootstrap tree, apply the clade rule of `step5_joint.define_joint_targets` and count how often each
             step has a target, how often the target equals the target of the best tree, and how often the smallest side that
             holds both heptaketide references does not hold the decoy wA (the wA question; the default-rule target of that
             step also needs that side to hold no other labelled leaf). Writes a TSV and a JSON.
  subset     Rarefaction input: keep all labelled leaves and a random fraction of the representatives; write the pruned tree,
             the pruned alignment and the kept representative ids.
  boundary   The labelled leaves outside a target that lie nearest to it (patristic distance).
  rep-coverage  Cluster sizes, and how many domains have a representative within the identity and coverage rule (mmseqs search table).
  compare    Compare calls (target / decoy / unassigned) of two runs of the same steps, by direction, on the proteins with a
             row of the step. Writes a TSV and a JSON.
Nothing here is an ortholog call. Support values say how often a split of the best tree is found in the bootstrap trees.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import random
import sys
from collections import Counter
from pathlib import Path

_HERE = Path(__file__).parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


sj = _load("step5_joint")
sp = sj.sp


def parse_tree(upo, newick):
    top = upo.parse_newick(newick)
    for i, n in enumerate(sp._nodes(top)[1:]):
        n.edge = i
    return top


def sides_of_tree(upo, newick):
    """(top, set of all sides) of a Newick string."""
    top = parse_tree(upo, newick)
    return top, sp.all_sides(upo.edge_bipartitions(top))


def edge_support(upo, newick):
    """{side: support float} for every internal edge with a label; both sides of an edge map to the value.
    Edges without a label (a leaf edge, or a root edge) are absent."""
    top = parse_tree(upo, newick)
    allleaves = frozenset(l.name for l in top.leaves())
    out = {}

    def walk(n):
        s = frozenset([n.name]) if not n.children else frozenset().union(*(walk(c) for c in n.children))
        if n.children and n.parent is not None and n.label not in (None, ""):
            try:
                v = float(n.label)
            except ValueError:
                v = None
            if v is not None:
                out[s] = v
                out[allleaves - s] = v
        return s
    walk(top)
    return out


def _step_decoys(step_decoys_fa):
    """{step: set of ids} from {step: FASTA}; None keeps the union scope."""
    return None if step_decoys_fa is None else {k: set(sp.ids_in_fasta(f)) for k, f in step_decoys_fa.items()}


def _ref_sets(step_refs_fa):
    return {s: set(sp.ids_in_fasta(f)) for s, f in step_refs_fa.items()}


def key_nodes(sides, step_refs, decoys, members, wa, include_variants, step_decoys=None):
    """Named nodes (name -> side or None) whose support is reported.
    target.<S>: the default-rule target of the step; all_references: the smallest side with every reference;
    hepta_refs / hepta_refs_with_wA: the smallest side with the heptaketide references / with them and wA; hepta_inner.<k>: the edges inside the smallest side with the heptaketide references and wA."""
    allrefs = set().union(*step_refs.values())
    targets = sj.define_joint_targets(sides, step_refs, decoys, set(), members, None, step_decoys)
    nodes = {f"target.{s}": t for s, t in sorted(targets.items())}
    nodes["all_references"] = sp.smallest_side(sides, allrefs, ())
    hep = step_refs.get("heptaketide_pks")
    if hep and wa:
        nodes["hepta_refs_with_wA"] = sp.smallest_side(sides, set(hep) | {wa}, ())
        nodes["hepta_refs"] = sp.smallest_side(sides, hep, ())
        for s, v in include_variants.items():
            nodes[f"variant_{s}_target"] = sj.define_joint_targets(sides, step_refs, decoys, set(), members,
                                                                    {"heptaketide_pks": set(v)}, step_decoys)["heptaketide_pks"]
    return nodes


def inner_sides(sides, outer):
    """Sides that are proper subsets of `outer` with at least 2 leaves (the clades inside it), smallest first."""
    if outer is None:
        return []
    return sorted((s for s in sides if s < outer and len(s) >= 2), key=lambda s: (len(s), sorted(s)))


def run_support(best_tree, fbp_tree, tbe_tree, step_refs_fa, decoys_fa, members_fa, wa, include, out_tsv, out_json, step_decoys_fa=None):
    upo = sp.load_upo()
    step_refs = _ref_sets(step_refs_fa)
    decoys = set(sp.ids_in_fasta(decoys_fa))
    members = set(sp.ids_in_fasta(members_fa)) if members_fa else set()
    top, sides = sides_of_tree(upo, Path(best_tree).read_text())
    fbp = edge_support(upo, Path(fbp_tree).read_text())
    tbe = edge_support(upo, Path(tbe_tree).read_text())
    nodes = key_nodes(sides, step_refs, decoys, members, wa, include, _step_decoys(step_decoys_fa))
    labelled = decoys | members | set().union(*step_refs.values())
    rows = []
    for name, side in nodes.items():
        rows.append(_support_row(name, side, fbp, tbe, labelled))
    outer = nodes.get("hepta_refs_with_wA")
    for i, s in enumerate(inner_sides(sides, outer), 1):
        rows.append(_support_row(f"hepta_inner.{i}", s, fbp, tbe, labelled))
    cols = ["node", "n_leaves", "labelled_leaves", "fbp", "tbe"]
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows([r[c] for c in cols] for r in rows)
    summary = {r["node"]: {"n_leaves": r["n_leaves"], "fbp": r["fbp"], "tbe": r["tbe"]} for r in rows}
    Path(out_json).write_text(json.dumps(summary, indent=1) + "\n")
    return rows


def _support_row(name, side, fbp, tbe, labelled):
    if side is None:
        return {"node": name, "n_leaves": "", "labelled_leaves": "none (no such side)", "fbp": "", "tbe": ""}
    f, t = fbp.get(side), tbe.get(side)
    # a side with all leaves but one has a leaf edge, not an internal edge: no support value exists
    return {"node": name, "n_leaves": len(side), "labelled_leaves": ",".join(sorted(side & labelled)),
            "fbp": "" if f is None else f, "tbe": "" if t is None else t}


def stability_of_tree(upo, newick, step_refs, decoys, members, wa, variants, step_decoys=None):
    """Facts about one tree: {step: target side or None}, hepta refs without wA exists, variant targets."""
    _, sides = sides_of_tree(upo, newick)
    targets = sj.define_joint_targets(sides, step_refs, decoys, set(), members, None, step_decoys)
    out = {"targets": targets}
    hep = step_refs.get("heptaketide_pks")
    if hep and wa:
        out["hepta_no_wA"] = wa not in sp.smallest_side(sides, hep, ())
    out["variants"] = {k: sj.define_joint_targets(sides, step_refs, decoys, set(), members, {"heptaketide_pks": set(v)}, step_decoys)["heptaketide_pks"]
                       for k, v in variants.items()}
    return out


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else 1.0


def run_stability(best_tree, bs_trees, step_refs_fa, decoys_fa, members_fa, wa, include, out_tsv, out_json, step_decoys_fa=None):
    upo = sp.load_upo()
    step_refs = _ref_sets(step_refs_fa)
    decoys = set(sp.ids_in_fasta(decoys_fa))
    members = set(sp.ids_in_fasta(members_fa)) if members_fa else set()
    sd = _step_decoys(step_decoys_fa)
    best = stability_of_tree(upo, Path(best_tree).read_text(), step_refs, decoys, members, wa, include, sd)
    trees = [ln.strip() for ln in open(bs_trees) if ln.strip()]
    n = len(trees)
    have, jac, jac_all = Counter(), {}, {}
    hep_no_wa, var_have = 0, Counter()
    for nw in trees:
        r = stability_of_tree(upo, nw, step_refs, decoys, members, wa, include, sd)
        for s, t in r["targets"].items():
            if t is not None:
                have[s] += 1
                if best["targets"][s] is not None:
                    jac.setdefault(s, []).append(jaccard(t, best["targets"][s]))
        hep_no_wa += int(r.get("hepta_no_wA", False))
        for k, t in r["variants"].items():
            var_have[k] += int(t is not None)
    rows, summ = [], {"n_bootstrap_trees": n, "decoy_scope": "union" if sd is None else "per_step", "steps": {}}
    for s in sorted(step_refs):
        j = sorted(jac.get(s, []))
        med = j[len(j) // 2] if j else None
        d = {"target_in_best_tree": best["targets"][s] is not None, "trees_with_target": have[s],
             "fraction_with_target": have[s] / n if n else None, "median_jaccard_with_best_target": med,
             "n_exact_same_target": sum(1 for x in j if x == 1.0)}
        summ["steps"][s] = d
        rows.append([s, "default rule target exists", have[s], n, f"{have[s] / n:.4f}" if n else ""])
        rows.append([s, "target leaf set equals the best-tree target", d["n_exact_same_target"], n, f"{d['n_exact_same_target'] / n:.4f}" if n else ""])
    if wa:
        summ["hepta_refs_clade_without_wA"] = {"trees": hep_no_wa, "of": n, "in_best_tree": best.get("hepta_no_wA")}
        rows.append(["heptaketide_pks", "smallest side with both heptaketide references does not hold wA", hep_no_wa, n, f"{hep_no_wa / n:.4f}" if n else ""])
    for k in sorted(var_have):
        summ.setdefault("variants", {})[k] = {"trees_with_target": var_have[k], "of": n, "in_best_tree": best["variants"][k] is not None}
        rows.append(["heptaketide_pks", f"variant {k} target exists", var_have[k], n, f"{var_have[k] / n:.4f}" if n else ""])
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["step", "measure", "trees", "of", "fraction"])
        w.writerows(rows)
    Path(out_json).write_text(json.dumps(summ, indent=1) + "\n")
    return summ


# ---------- rarefaction ----------

def choose_subset(rep_ids, fraction, seed):
    """Sorted random subset of `rep_ids` (size round(fraction*n), at least 1)."""
    ids = sorted(rep_ids)
    k = max(1, round(len(ids) * fraction))
    return sorted(random.Random(seed).sample(ids, k))


def prune_tree(upo, newick, drop):
    top = upo.parse_newick(newick)
    for leaf in sorted(drop):
        top = upo.prune(top, leaf)
    return upo.to_newick(top)


def run_subset(tree, msa, rep_ids_file, fraction, seed, out_tree, out_msa, out_reps):
    upo = sp.load_upo()
    reps = [ln.strip() for ln in open(rep_ids_file) if ln.strip()]
    top = upo.parse_newick(Path(tree).read_text())
    leaves = {l.name for l in top.leaves()}
    missing = set(reps) - leaves
    if missing:
        raise SystemExit(f"ERROR: representatives not in the tree: {sorted(missing)[:5]}")
    keep = set(choose_subset(reps, fraction, seed))
    drop = set(reps) - keep
    Path(out_tree).write_text(prune_tree(upo, Path(tree).read_text(), drop) + "\n")
    n = 0
    with open(out_msa, "w") as fo:
        for name, seq in upo.read_fasta(msa):
            if name not in drop:
                upo.write_fasta(fo, name, seq)
                n += 1
    Path(out_reps).write_text("".join(f"{r}\n" for r in sorted(keep)))
    return {"representatives_kept": len(keep), "representatives_dropped": len(drop), "msa_sequences": n, "fraction": fraction, "seed": seed}


# ---------- boundary of a target ----------

def target_boundary(upo, newick, target, roles, k=8):
    """Labelled leaves OUTSIDE `target` nearest to it: [(leaf, role, min patristic distance to any target leaf, nearest target leaf)].
    roles: {leaf: role}; a leaf with role 'representative' or 'other' is not labelled."""
    top = upo.parse_newick(newick)
    paths = sp.leaf_depths(top)
    lab = sorted(l for l, r in roles.items() if r not in ("representative", "other") and l not in target and l in paths)
    out = []
    for l in lab:
        d, t = min((sp.tree_distance(paths, l, t), t) for t in sorted(target) if t in paths)
        out.append((l, roles[l], d, t))
    return sorted(out, key=lambda x: (x[2], x[0]))[:k]


def run_boundary(tree, clade_members, out_tsv, k=8):
    upo = sp.load_upo()
    rows = sp.read_tsv(clade_members)
    target = {r["label"] for r in rows if r["our_clade"] == "target"}
    roles = {r["label"]: r["role"] for r in rows}
    res = target_boundary(upo, Path(tree).read_text(), target, roles, k)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["outside_labelled_leaf", "role", "min_distance_to_target", "nearest_target_leaf", "target_size"])
        for l, r, d, t in res:
            w.writerow([l, r, f"{d:.4f}", t, len(target)])
    return res


# ---------- representative coverage and cluster sizes ----------

def cluster_size_stats(clusters):
    """clusters: {rep: [members]} -> sizes summary."""
    sz = sorted(len(m) for m in clusters.values())
    n = sum(sz)
    return {"clusters": len(sz), "domains": n, "singletons": sum(1 for x in sz if x == 1), "largest": sz[-1] if sz else 0,
            "median_size": sz[len(sz) // 2] if sz else 0, "clusters_ge_10": sum(1 for x in sz if x >= 10),
            "domains_in_clusters_ge_10": sum(x for x in sz if x >= 10), "domains_in_largest_10": sum(sz[-10:])}


def rep_coverage(clusters, hits, min_id, min_cov):
    """hits: iterable of (query, target, fident, qcov, tcov) from mmseqs search of all domains against the representatives.
    A hit passes if fident >= min_id and qcov >= min_cov and tcov >= min_cov (the clustering rule, --cov-mode 0).
    Returns counts for the assigned representative and for any representative."""
    assigned = {m: r for r, ms in clusters.items() for m in ms}
    best_any, best_assigned, passing_any, passing_assigned = {}, {}, set(), set()
    for q, t, fi, qc, tc in hits:
        ok = fi >= min_id and qc >= min_cov and tc >= min_cov
        if ok:
            passing_any.add(q)
            if assigned.get(q) == t:
                passing_assigned.add(q)
        if q not in best_any or fi > best_any[q]:
            best_any[q] = fi
        if assigned.get(q) == t:
            best_assigned[q] = fi
    n = len(assigned)
    non_self = [q for q in assigned if assigned[q] != q]
    ident = sorted(best_assigned[q] for q in non_self if q in best_assigned)
    pct = lambda p: ident[min(len(ident) - 1, int(p * len(ident)))] if ident else None
    return {"domains": n, "domains_that_are_representatives": n - len(non_self), "non_representative_domains": len(non_self),
            "with_hit_to_any_rep": len(set(assigned) & set(best_any)), "pass_rule_any_rep": len(passing_any & set(assigned)),
            "pass_rule_assigned_rep": len(passing_assigned & set(assigned)),
            "non_rep_pass_rule_assigned_rep": len(passing_assigned & set(non_self)),
            "non_rep_identity_to_assigned_rep_min": ident[0] if ident else None, "p1": pct(0.01), "p5": pct(0.05), "median": pct(0.5),
            "non_rep_below_min_id_assigned": sum(1 for x in ident if x < min_id),
            "non_rep_without_hit_to_assigned_rep": len(non_self) - len(ident), "min_id": min_id, "min_cov": min_cov}


def read_hits(path):
    with sp.open_text(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 5:
                yield f[0], f[1], float(f[2]), float(f[3]), float(f[4])


# ---------- call comparison ----------

def read_calls(assign_tsv, proteins=None):
    """{query: our_clade} (target / decoy / unassigned) from an upo_place assign table; optionally limited to `proteins`."""
    out = {}
    for r in sp.read_tsv(assign_tsv):
        if proteins is None or r["query"] in proteins:
            out[r["query"]] = r["our_clade"]
    return out


def compare_calls(old, new):
    """Counts of (old call, new call) over the proteins in both. Returns (Counter, summary dict)."""
    common = set(old) & set(new)
    c = Counter((old[p], new[p]) for p in common)
    t_to_not = sum(v for (a, b), v in c.items() if a == "target" and b != "target")
    not_to_t = sum(v for (a, b), v in c.items() if a != "target" and b == "target")
    return c, {"proteins_compared": len(common), "only_in_old": len(set(old) - set(new)), "only_in_new": len(set(new) - set(old)),
               "target_to_not_target": t_to_not, "not_target_to_target": not_to_t,
               "changed_target_status": t_to_not + not_to_t,
               "fraction_changed_target_status": (t_to_not + not_to_t) / len(common) if common else None,
               "any_class_change": sum(v for (a, b), v in c.items() if a != b),
               "old_target": sum(v for (a, _), v in c.items() if a == "target"),
               "new_target": sum(v for (_, b), v in c.items() if b == "target")}


def run_compare(old_assign, new_assign, step_proteins_file, out_tsv, out_json, subset_file=None):
    proteins = {ln.split("\t")[0] for ln in open(step_proteins_file).read().splitlines()[1:] if ln.strip()}
    if subset_file:
        proteins &= {ln.strip() for ln in open(subset_file) if ln.strip()}
    c, summ = compare_calls(read_calls(old_assign, proteins), read_calls(new_assign, proteins))
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["old_call", "new_call", "proteins"])
        for (a, b), v in sorted(c.items()):
            w.writerow([a, b, v])
    Path(out_json).write_text(json.dumps(summ, indent=1) + "\n")
    return summ


def _kv(items):
    return dict(x.split("=", 1) for x in items)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("support", "stability"):
        p = sub.add_parser(name)
        p.add_argument("--best-tree", required=True)
        p.add_argument("--refs", nargs="+", required=True, help="STEP=refs.faa")
        p.add_argument("--decoys", required=True)
        p.add_argument("--members-fa")
        p.add_argument("--decoy-scope", choices=["per-step", "union"], default="per-step")
        p.add_argument("--step-decoys", nargs="*", default=[], help="STEP=decoy FASTA, required for per-step")
        p.add_argument("--wa", default="Q03149", help="id of the decoy wA")
        p.add_argument("--include", nargs="*", default=[], help="NAME=ID,ID heptaketide variants (leaves allowed in the target)")
        p.add_argument("--out-tsv", required=True)
        p.add_argument("--out-json", required=True)
        if name == "support":
            p.add_argument("--fbp-tree", required=True)
            p.add_argument("--tbe-tree", required=True)
        else:
            p.add_argument("--bs-trees", required=True)
    s = sub.add_parser("subset")
    s.add_argument("--tree", required=True)
    s.add_argument("--msa", required=True)
    s.add_argument("--rep-ids", required=True)
    s.add_argument("--fraction", type=float, required=True)
    s.add_argument("--seed", type=int, required=True)
    s.add_argument("--out-tree", required=True)
    s.add_argument("--out-msa", required=True)
    s.add_argument("--out-reps", required=True)
    b = sub.add_parser("boundary")
    b.add_argument("--tree", required=True)
    b.add_argument("--clade-members", required=True, help="clade_members.<STEP>.tsv")
    b.add_argument("--k", type=int, default=8)
    b.add_argument("--out-tsv", required=True)
    k = sub.add_parser("rep-coverage")
    k.add_argument("--clusters", required=True)
    k.add_argument("--hits", required=True, help="mmseqs search table: query target fident qcov tcov")
    k.add_argument("--min-id", type=float, default=0.8)
    k.add_argument("--min-cov", type=float, default=0.8)
    k.add_argument("--out-json", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--old-assign", required=True)
    c.add_argument("--new-assign", required=True)
    c.add_argument("--step-proteins", required=True, help="candidates_step.<S>.tsv (first column: protein id)")
    c.add_argument("--subset", help="file with one protein id per line to restrict to")
    c.add_argument("--out-tsv", required=True)
    c.add_argument("--out-json", required=True)
    a = ap.parse_args(argv)
    sdf = None
    if a.cmd in ("support", "stability") and a.decoy_scope == "per-step":
        sdf = _kv(a.step_decoys)
        if not sdf:
            raise SystemExit("ERROR: --decoy-scope per-step needs --step-decoys STEP=FASTA (or use --decoy-scope union)")
    if a.cmd == "support":
        run_support(a.best_tree, a.fbp_tree, a.tbe_tree, _kv(a.refs), a.decoys, a.members_fa, a.wa,
                    {k: v.split(",") for k, v in _kv(a.include).items()}, a.out_tsv, a.out_json, sdf)
    elif a.cmd == "stability":
        print(json.dumps(run_stability(a.best_tree, a.bs_trees, _kv(a.refs), a.decoys, a.members_fa, a.wa,
                                       {k: v.split(",") for k, v in _kv(a.include).items()}, a.out_tsv, a.out_json, sdf)))
    elif a.cmd == "subset":
        print(json.dumps(run_subset(a.tree, a.msa, a.rep_ids, a.fraction, a.seed, a.out_tree, a.out_msa, a.out_reps)))
    elif a.cmd == "boundary":
        run_boundary(a.tree, a.clade_members, a.out_tsv, a.k)
    elif a.cmd == "rep-coverage":
        cl = sp.parse_clusters(a.clusters)
        out = {"cluster_sizes": cluster_size_stats(cl), "coverage": rep_coverage(cl, read_hits(a.hits), a.min_id, a.min_cov)}
        Path(a.out_json).write_text(json.dumps(out, indent=1) + "\n")
        print(json.dumps(out))
    else:
        print(json.dumps(run_compare(a.old_assign, a.new_assign, a.step_proteins, a.out_tsv, a.out_json, a.subset)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
