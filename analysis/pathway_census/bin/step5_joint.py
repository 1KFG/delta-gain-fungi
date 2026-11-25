#!/usr/bin/env python3.12
"""Step 5 with ONE tree for SEVERAL steps (used for the three PKS steps of melanin_dhn). Run after slurm/step5_prepare_joint.sbatch.

The step5_place.py tools make one tree per step and one target clade per tree. A joint tree holds the references of all
steps as labelled leaves. This script gives each step its own target clade on that tree and keeps the other tools unchanged
(it imports step5_place.py and step5_single_ref_clades.py).

Subcommands:
  clades            target clade of every step. Rule (PROPOSED; the rule of bin/step5_single_ref_clades.py, for one or more
                    references): target(S) = the LARGEST tree side that holds all references of S and no decoy, outgroup,
                    named member or reference of another step. If no such side exists the references of S are not isolated
                    on the tree; that is a result: S gets no clade file, and the report says so. Exit code 3 when NO step has
                    a target. Writes clade_members.<S>.tsv per step, clade_members.tsv (union view: target = in any target),
                    ref.rooted.nwk, clade_report.json, target_overlap.tsv. The targets of two steps cannot share a leaf
                    when a decoy exists and the scope is UNION (a side holding the references of S and none of T, and a side holding
                    those of T and none of S, cannot cover the tree while a decoy lies outside both). With --decoy-scope per-step
                    (the default; the decoys of each step bound only that step) targets may share leaves; the report counts the pairs.
  filter-placements per step S: placements.<S>.tsv, placements_detail.<S>.tsv and candidates_step.<S>.tsv, limited to the proteins
                    that have a row of S in the candidate folder. The files have as many proteins as that folder has for S.
  report            tables for the placed run (calls per step, hit pattern against call, nearest labelled leaves, genomes by class).
Nothing here is an ortholog call. The target is a tree region bounded by the labelled leaves that were put in the tree.
"""
from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import json
import subprocess
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
sr = _load("step5_single_ref_clades")


def define_joint_targets(sides, step_refs, decoys, outgroups, members, include=None, step_decoys=None):
    """{step: target side or None}. step_refs: {step: set of reference ids}.
    include: optional {step: set of labelled leaves (decoys, named members) allowed inside the target of that step}.
    It is a VARIANT, not the default: it widens the target past leaves that the default rule treats as the boundary.
    step_decoys: optional {step: set of decoy ids}. Decoy scope PER STEP (owner decision 2026-11-24): only the decoys of the step itself
    bound its target; a decoy of another step is not a boundary (it may lie inside the target). Without it (scope UNION, the first
    PKS run) every decoy bounds every step. Other steps' references, named leaves and outgroups bound the target in both scopes."""
    allrefs = set().union(*step_refs.values()) if step_refs else set()
    include = include or {}
    out = {}
    for s, refs in step_refs.items():
        d = set(decoys) if step_decoys is None else set(step_decoys.get(s, ()))
        labelled = d | set(outgroups) | set(members)
        out[s] = sr.largest_free_side(sides, refs, (labelled | (allrefs - refs)) - set(include.get(s, ())))
    return out


def _write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def run_clades(tree, step_refs_fa, decoys_fa, members_fa, outgroup_fa, reps_ids, out_dir, joint_id, outgroup_role="tree_outgroup",
               include=None, step_decoys_fa=None):
    """step_decoys_fa: {step: decoy FASTA}; given = decoy scope per step, None = scope union (every decoy bounds every step)."""
    upo = sp.load_upo()
    out_dir = Path(out_dir)
    step_refs = {s: set(sp.ids_in_fasta(f)) for s, f in step_refs_fa.items()}
    refs = set().union(*step_refs.values())
    if sum(len(v) for v in step_refs.values()) != len(refs):
        raise SystemExit("ERROR: a reference id is listed for two steps")
    decoys = set(sp.ids_in_fasta(decoys_fa)) if decoys_fa else set()
    members = set(sp.ids_in_fasta(members_fa)) if members_fa else set()
    outgroups = set(sp.ids_in_fasta(outgroup_fa)) if outgroup_fa else set()
    if members & (refs | decoys | outgroups) or decoys & refs or outgroups & (refs | decoys):
        raise SystemExit("ERROR: a leaf has two kinds (reference, decoy, named member, outgroup)")
    reps = {ln.strip() for ln in open(reps_ids) if ln.strip()}
    top = upo.parse_newick(open(tree).read())
    leaves = {l.name for l in top.leaves()}
    missing = (refs | decoys | members | outgroups) - leaves
    if missing:
        raise SystemExit(f"ERROR: leaves missing from the tree: {sorted(missing)}")
    for i, n in enumerate(sp._nodes(top)[1:]):
        n.edge = i
    bip = upo.edge_bipartitions(top)
    sides = sp.all_sides(bip)
    everything = frozenset(leaves)
    roles = sp.leaf_roles(leaves, refs, decoys, outgroups, reps, members, outgroup_role)
    include = {k: set(v) for k, v in (include or {}).items()}
    bad = set().union(*include.values()) - (decoys | members) if include else set()
    if bad:
        raise SystemExit(f"ERROR: --include leaves are not decoys or named members: {sorted(bad)}")
    step_decoys = None
    if step_decoys_fa is not None:
        step_decoys = {k: set(sp.ids_in_fasta(f)) for k, f in step_decoys_fa.items()}
        stray = set().union(*step_decoys.values()) - decoys if step_decoys else set()
        if stray:
            raise SystemExit(f"ERROR: step decoys not in the --decoys file: {sorted(stray)}")
        for k in step_refs:
            step_decoys.setdefault(k, set())
    targets = define_joint_targets(sides, step_refs, decoys, outgroups, members, include, step_decoys)
    rooting = outgroups if outgroups else decoys if decoys else members
    (out_dir / "ref.rooted.nwk").write_text(sp._reroot(upo, top, rooting) + "\n")
    rep = {"joint_id": joint_id, "n_leaves": len(leaves), "steps": sorted(step_refs),
           "target_rule": "largest side with all references of the step and no decoy, outgroup, named member or reference of another step",
           "rooting_leaves": sorted(rooting),
           "rooting_source": "outgroup file" if outgroups else "decoys" if decoys else "named members",
           "n_references": len(refs), "n_decoys": len(decoys), "n_named_members": len(members), "n_outgroups": len(outgroups),
           "n_representatives_in_tree": len(reps & leaves),
           "references": {s: sorted(v) for s, v in step_refs.items()},
           "decoy_scope": "per_step" if step_decoys is not None else "union",
           "step_decoys": {k: sorted(v) for k, v in step_decoys.items()} if step_decoys is not None else None,
           "variant_include": {k: sorted(v) for k, v in include.items()}, "per_step": {}}
    for s, refs_s in sorted(step_refs.items()):
        smallest = sp.smallest_side(sides, refs_s, ())
        t = targets[s]
        d = {"n_references": len(refs_s), "refs_isolated": t is not None,
             "smallest_side_with_refs_size": len(smallest),
             "smallest_side_labelled_leaves": sorted(l for l in smallest if roles[l] != "representative")}
        if t is not None:
            d.update(target_size=len(t), target_references=len(refs_s & t), target_representatives=len(reps & t),
                     target_named_members=len(members & t), target_decoys=len(decoys & t),
                     decoy_side_size=len(everything - t - outgroups - members),
                     decoy_side_representatives=len(reps & (everything - t - outgroups - members)))
            clades = {"target": t, "decoy": frozenset(everything - t - outgroups - members), "outgroup": frozenset(outgroups)}
            sp.write_clade_members(out_dir / f"clade_members.{s}.tsv", clades, roles)
        rep["per_step"][s] = d
    live = {s: t for s, t in targets.items() if t is not None}
    ov = []
    names = sorted(live)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ov.append([a, b, len(live[a]), len(live[b]), len(live[a] & live[b])])
    _write_tsv(out_dir / "target_overlap.tsv", ["step_a", "step_b", "size_a", "size_b", "shared_leaves"], ov)
    rep["target_pairs_sharing_a_leaf"] = sum(1 for r in ov if r[4] > 0)
    if live:
        union = frozenset().union(*live.values())
        clades = {"target": union, "decoy": frozenset(everything - union - outgroups - members), "outgroup": frozenset(outgroups)}
        sp.write_clade_members(out_dir / "clade_members.tsv", clades, roles)
    (out_dir / "clade_report.json").write_text(json.dumps(rep, indent=1) + "\n")
    return 0 if live else 3


def read_candidate_rows(cand_dir, step):
    """Rows (dicts) of one step from every chunk file of a candidate folder; plus the header list."""
    rows, header = [], None
    for p in sorted(glob.glob(f"{cand_dir}/candidates.chunk*.PROVISIONAL.tsv.zst")):
        txt = subprocess.run(["zstd", "-dc", p], check=True, capture_output=True).stdout.decode().splitlines()
        lines = [l for l in txt if not l.startswith("#") and l.strip()]
        if not lines:
            continue
        header = header or lines[0].split("\t")
        for l in lines[1:]:
            f = dict(zip(header, l.split("\t")))
            if f["step_id"] == step:
                rows.append(f)
    return rows, header


def run_filter_placements(assign_tsv, step, cand_dir, out_plc, out_detail, out_cand):
    rows, header = read_candidate_rows(cand_dir, step)
    hit = {r["protein_id"] for r in rows}
    assign = [r for r in sp.read_tsv(assign_tsv) if r["query"] in hit]
    plc, det = sp.to_placements(assign, step)
    sp._write_dicts(out_plc, plc, ["protein_id", "step_id", "clade", "support"])
    sp._write_dicts(out_detail, det, list(det[0]) if det else ["protein_id"])
    with open(out_cand, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in header) + "\n")
    return {"step": step, "candidate_rows": len(rows), "candidate_proteins": len(hit), "placements": len(plc),
            "queries_missing_from_assign": len(hit - {r["query"] for r in assign})}


def class_of(genomes_tsv):
    out = {}
    for r in sp.read_tsv(genomes_tsv):
        out[r["locustag"]] = (r.get("class_") or "(no label)", r.get("genus") or "(no label)", r.get("phylum") or "(no label)",
                              r.get("has_proteins", "").strip().lower() == "true")
    return out


def run_report(out_dir, scores_tsv, genomes_tsv, min_class_genomes=20, n_neighbours=10):
    d = Path(out_dir)
    rep_dir = d / "joint_report"
    rep_dir.mkdir(exist_ok=True)
    upo = sp.load_upo()
    top = upo.parse_newick((d / "ref.raxml.bestTree").read_text())
    cr = json.loads((d / "clade_report.json").read_text())
    steps = cr["steps"]
    role = {r["label"]: r["role"] for r in sp.read_tsv(d / "clade_members.tsv")}
    labelled = [l for l, r in role.items() if r != "representative"]
    reps = [l for l, r in role.items() if r == "representative"]
    info = class_of(genomes_tsv)
    size = {rep: len(m) for rep, m in sp.parse_clusters(d / "clusters.tsv.zst").items()}
    ref_step = {r: s for s, ids in cr["references"].items() for r in ids}
    gtax = lambda leaf: info.get(leaf.split("_")[0], ("(no label)", "(no label)", "(no label)", False))
    # calls per protein per step clade (all proteins of the union)
    scores = {r["protein_id"]: r for r in sp.read_tsv(scores_tsv)}
    call = {}
    for s in steps:
        p = d / f"assign.{s}.tsv"
        if p.exists():
            call[s] = {r["query"]: r for r in sp.read_tsv(p)}
    live = sorted(call)
    # 1. calls per step (proteins with a row of the step), assign_steps result, genomes
    rows = []
    hit_by_step = {s: {p for p, r in scores.items() if r.get(f"score_{s}")} for s in steps}
    for s in steps:
        if s not in call:
            rows.append([s, len(hit_by_step[s]), "no clade (references not isolated)", "", "", "", "", ""])
            continue
        c = Counter(call[s][p]["our_clade"] for p in hit_by_step[s])
        asg = d / f"assigned.{s}.PROVISIONAL.tsv"
        a = sp.read_tsv(asg) if asg.exists() else []
        ap = {r["protein_id"] for r in a if r["label"] == "assigned"}
        rows.append([s, len(hit_by_step[s]), c.get("target", 0), c.get("decoy", 0) + c.get("outgroup", 0), c.get("unassigned", 0),
                     len(ap), len({p.split("_")[0] for p in ap}), len({p.split("_")[0] for p in hit_by_step[s]})])
    _write_tsv(rep_dir / "calls_per_step.tsv", ["step", "proteins_with_hit", "target", "other_side", "unassigned",
                                                 "assigned_proteins(assign_steps)", "genomes_with_assigned", "genomes_with_hit"], rows)
    # 2. hit pattern against the call of the union of proteins
    pat = Counter()
    for p, r in scores.items():
        hits = "+".join(s for s in steps if r.get(f"score_{s}")) or "none"
        tg = "+".join(s for s in live if call[s].get(p, {}).get("our_clade") == "target") or "no_target"
        pat[(hits, tg)] += 1
    _write_tsv(rep_dir / "hit_pattern_vs_target.tsv", ["hmm_hits", "target_of", "proteins"],
               [[h, t, n] for (h, t), n in sorted(pat.items(), key=lambda x: (-x[1], x[0]))])
    n_multi = sum(1 for p in scores if sum(call[s].get(p, {}).get("our_clade") == "target" for s in live) > 1)
    # 3. nearest labelled leaves of every reference, with distance
    paths = sp.leaf_depths(top)
    nn_rows = []
    refs = [l for l in labelled if role[l] == "reference"]
    for rf in sorted(refs):
        dist = sorted(((sp.tree_distance(paths, rf, o), o) for o in labelled + reps if o != rf), key=lambda x: (x[0], x[1]))
        for rank, (dd, o) in enumerate(dist[:n_neighbours], 1):
            g = gtax(o) if role[o] == "representative" else ("", "", "", False)
            nn_rows.append([rf, ref_step.get(rf, ""), rank, o, role[o], f"{dd:.4f}", g[1] if role[o] == "representative" else "",
                            g[0] if role[o] == "representative" else "", size.get(o, "")])
    _write_tsv(rep_dir / "reference_neighbours.tsv", ["reference", "step", "rank", "neighbour", "neighbour_role", "distance",
                                                       "genus", "class", "n_candidates_in_cluster"], nn_rows)
    ld = []
    for rf in sorted(refs):
        for o in labelled:
            if o != rf:
                ld.append([rf, o, role[o], f"{sp.tree_distance(paths, rf, o):.4f}"])
    _write_tsv(rep_dir / "labelled_distances.tsv", ["reference", "leaf", "role", "distance"], ld)
    # nearest labelled leaf per representative
    near = []
    for q, n, dd in sp.nearest_labelled(top, reps, labelled):
        g = gtax(q)
        near.append([q, n, role[n], f"{dd:.4f}", size.get(q, ""), g[1], g[0]])
    _write_tsv(rep_dir / "nearest_labelled.tsv", ["representative", "nearest_labelled", "role", "distance", "n_candidates_in_cluster",
                                                   "genus", "class"], near)
    # 4. genomes by class: with a hit of the step, with an assigned call
    genomes = [g for g, v in info.items() if v[3]]
    out = []
    for s in live:
        hit_g = {p.split("_")[0] for p in hit_by_step[s]}
        asg = d / f"assigned.{s}.PROVISIONAL.tsv"
        asg_g = {r["protein_id"].split("_")[0] for r in sp.read_tsv(asg) if r["label"] == "assigned"} if asg.exists() else set()
        byc = defaultdict(list)
        for g in genomes:
            byc[info[g][0]].append(g)
        for c, gs in sorted(byc.items(), key=lambda x: -len(x[1])):
            if len(gs) >= min_class_genomes:
                nh, na = sum(g in hit_g for g in gs), sum(g in asg_g for g in gs)
                out.append([s, c, len(gs), nh, na, f"{nh / len(gs):.4f}", f"{na / len(gs):.4f}"])
        out.append([s, "ALL", len(genomes), len(hit_g & set(genomes)), len(asg_g & set(genomes)),
                    f"{len(hit_g & set(genomes)) / len(genomes):.4f}", f"{len(asg_g & set(genomes)) / len(genomes):.4f}"])
    _write_tsv(rep_dir / "genomes_by_class.tsv", ["step", "class", "n_genomes", "with_hit", "with_assigned", "share_hit", "share_assigned"], out)
    summ = {"steps": steps, "steps_with_clade": live, "proteins_target_for_more_than_one_step": n_multi,
            "classes_listed_min_genomes": min_class_genomes, "genomes_with_proteins": len(genomes)}
    (rep_dir / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("clades")
    c.add_argument("--tree", required=True)
    c.add_argument("--refs", required=True, nargs="+", help="STEP=reference FASTA, one per step")
    c.add_argument("--decoys")
    c.add_argument("--members-fa")
    c.add_argument("--outgroup")
    c.add_argument("--outgroup-role", default="tree_outgroup")
    c.add_argument("--include", action="append", default=[], help="VARIANT: STEP=ID,ID decoys or named members allowed inside the target of STEP")
    c.add_argument("--decoy-scope", choices=["per-step", "union"], default="per-step",
                   help="per-step (default): only the decoys of a step bound its target, from --step-decoys; union: every decoy bounds every step (first PKS run)")
    c.add_argument("--step-decoys", nargs="*", default=[], help="STEP=decoy FASTA (refs/<pathway>/decoys/<step>.faa), required for per-step")
    c.add_argument("--reps", required=True)
    c.add_argument("--joint-id", required=True)
    c.add_argument("--out-dir", required=True)
    f = sub.add_parser("filter-placements")
    f.add_argument("--assign", required=True)
    f.add_argument("--step", required=True)
    f.add_argument("--cand-dir", required=True)
    f.add_argument("--out-plc", required=True)
    f.add_argument("--out-detail", required=True)
    f.add_argument("--out-cand", required=True)
    r = sub.add_parser("report")
    r.add_argument("--out-dir", required=True)
    r.add_argument("--scores", required=True)
    r.add_argument("--genomes", required=True)
    r.add_argument("--min-class-genomes", type=int, default=20)
    a = ap.parse_args(argv)
    if a.cmd == "clades":
        step_refs = dict(x.split("=", 1) for x in a.refs)
        inc = {k: v.split(",") for k, v in (x.split("=", 1) for x in a.include)}
        sd = None
        if a.decoy_scope == "per-step":
            sd = dict(x.split("=", 1) for x in a.step_decoys)
            if not sd:
                raise SystemExit("ERROR: --decoy-scope per-step needs --step-decoys STEP=FASTA ... (or use --decoy-scope union)")
        return run_clades(a.tree, step_refs, a.decoys, a.members_fa, a.outgroup, a.reps, a.out_dir, a.joint_id, a.outgroup_role, inc, sd)
    if a.cmd == "filter-placements":
        print(json.dumps(run_filter_placements(a.assign, a.step, a.cand_dir, a.out_plc, a.out_detail, a.out_cand)))
        return 0
    print(json.dumps(run_report(a.out_dir, a.scores, a.genomes, a.min_class_genomes)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
