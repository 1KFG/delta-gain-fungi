#!/usr/bin/env python3.12
"""Clade definition for a step with ONE reference (step 5), run after slurm/step5_prepare.sbatch.

Why: bin/step5_place.py `clades` defines `target` as the SMALLEST tree side that holds all references and no decoy or
outgroup. With one reference that side is the reference leaf itself (size 1). No placement can fall inside a leaf, so every
candidate is called `decoy` or `unassigned`. This script uses another rule (PROPOSED):
  target = the LARGEST tree side that holds the reference and none of the labelled leaves of other kinds
           (decoys, outgroups, named members). It is the clade of the reference bounded by the nearest labelled neighbours.
The sides that hold the reference and no labelled leaf are nested (see `largest_free_side`), so the largest one is unique.
`decoy` = every leaf that is not in `target`, not an outgroup and not a named member (as in step5_place.define_clades).
The result depends on which labelled leaves are in the tree. With no labelled leaf, `target` is the whole tree: the script
stops then (exit 4), because the call carries no information.
Writes clade_members.tsv, clade_report.json, ref.rooted.nwk in --out-dir. The files made by step5_place `clades`
are kept as clade_members.smallest.tsv and clade_report.smallest.json (once; never overwritten).
Exit codes: 0 ok; 4 no labelled leaf other than the reference.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("step5_place", Path(__file__).parent / "step5_place.py")
sp = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sp)


def largest_free_side(sides, must_have, excluded):
    """Largest side that holds all of `must_have` and none of `excluded`; None if absent.
    Such sides are nested: two of them share the reference, are not disjoint, and their union cannot be the whole tree
    (neither holds an excluded leaf), so one contains the other. Ties in size go to the side whose sorted names sort first."""
    must_have, excluded = frozenset(must_have), frozenset(excluded)
    ok = [s for s in sides if must_have <= s and not (s & excluded)]
    return min(ok, key=lambda x: (-len(x), sorted(x))) if ok else None


def define_single_ref_clades(bip, refs, decoys, outgroups, members=()):
    """Return the same dict as step5_place.define_clades, with `target` by the largest-side rule."""
    sides = sp.all_sides(bip)
    everything = frozenset().union(*sides) if sides else frozenset()
    refs, decoys, outgroups, members = map(frozenset, (refs, decoys, outgroups, members))
    target = largest_free_side(sides, refs, decoys | outgroups | members)
    if target is None:
        decoy = decoys - outgroups
    else:
        decoy = everything - target - outgroups - members
    shared = bool(decoys) and sp.smallest_side(sides, decoys, refs) is not None
    return {"target": target, "decoy": frozenset(decoy), "outgroup": outgroups,
            "decoys_share_side_without_refs": shared}


def run(tree, refs_fa, decoys_fa, outgroup_fa, members_fa, reps_ids, out_dir, step_id, outgroup_role="tree_outgroup"):
    upo = sp.load_upo()
    out_dir = Path(out_dir)
    refs = set(sp.ids_in_fasta(refs_fa))
    decoys = set(sp.ids_in_fasta(decoys_fa)) if decoys_fa else set()
    outgroups = set(sp.ids_in_fasta(outgroup_fa)) if outgroup_fa else set()
    members = set(sp.ids_in_fasta(members_fa)) if members_fa else set()
    if len(refs) != 1:
        raise SystemExit(f"ERROR: this script is for one reference, found {len(refs)}")
    if members & (refs | decoys | outgroups):
        raise SystemExit("ERROR: a named member is also a reference, decoy or outgroup")
    if not (decoys | outgroups | members):
        print("ERROR: no labelled leaf other than the reference; the target would be the whole tree", file=sys.stderr)
        return 4
    reps = {ln.strip() for ln in open(reps_ids) if ln.strip()}
    top = upo.parse_newick(open(tree).read())
    leaves = {l.name for l in top.leaves()}
    missing = (refs | decoys | outgroups | members) - leaves
    if missing:
        raise SystemExit(f"ERROR: leaves missing from the tree: {sorted(missing)}")
    for i, n in enumerate(sp._nodes(top)[1:]):
        n.edge = i
    bip = upo.edge_bipartitions(top)
    sides = sp.all_sides(bip)
    for fn in ("clade_members.tsv", "clade_report.json"):
        src = out_dir / fn
        dst = out_dir / fn.replace(".", ".smallest.", 1)
        if src.exists() and not dst.exists():
            shutil.copy(src, dst)
    clades = define_single_ref_clades(bip, refs, decoys, outgroups, members)
    roles = sp.leaf_roles(leaves, refs, decoys, outgroups, reps, members, outgroup_role)
    rooting = outgroups if outgroups else decoys if decoys else members
    (out_dir / "ref.rooted.nwk").write_text(sp._reroot(upo, top, rooting) + "\n")
    t = clades["target"]
    rep = {"step_id": step_id, "n_leaves": len(leaves), "target_rule": "largest side with the reference and no labelled leaf of another kind",
           "rooting_leaves": sorted(rooting), "rooting_source": "outgroup file" if outgroups else "decoys" if decoys else "named members",
           "n_references": len(refs), "n_decoys": len(decoys), "n_outgroups": len(outgroups), "n_named_members": len(members),
           "refs_form_clade_free_of_decoys": t is not None,
           "decoys_share_side_without_refs": clades["decoys_share_side_without_refs"]}
    rep.update(sp.describe_outgroups(sides, outgroups, decoys, roles))
    rep["target_size"] = len(t)
    rep["target_refs"] = len(refs & t)
    rep["target_representatives"] = len(reps & t)
    rep["target_named_members"] = len(members & t)   # 0 by construction
    rep["target_is_reference_only"] = len(t) == 1
    rep["decoy_size"] = len(clades["decoy"])
    rep["decoy_representatives"] = len(reps & clades["decoy"])
    sp.write_clade_members(out_dir / "clade_members.tsv", clades, roles)
    (out_dir / "clade_report.json").write_text(json.dumps(rep, indent=1) + "\n")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--refs", required=True)
    ap.add_argument("--decoys")
    ap.add_argument("--outgroup")
    ap.add_argument("--members-fa")
    ap.add_argument("--outgroup-role", default="tree_outgroup")
    ap.add_argument("--reps", required=True)
    ap.add_argument("--step", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    return run(a.tree, a.refs, a.decoys, a.outgroup, a.members_fa, a.reps, a.out_dir, a.step, a.outgroup_role)


if __name__ == "__main__":
    sys.exit(main())
