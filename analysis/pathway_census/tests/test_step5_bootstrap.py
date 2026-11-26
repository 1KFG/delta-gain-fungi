"""Tests for bin/step5_bootstrap.py. Synthetic trees only."""
import importlib.util
import json
from pathlib import Path

import pytest

BIN = Path(__file__).parents[1] / "bin"
spec = importlib.util.spec_from_file_location("step5_bootstrap", BIN / "step5_bootstrap.py")
sb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sb)
sp = sb.sp
needs_upo = pytest.mark.skipif(not sp.UPO_BIN.exists(), reason="enzymes_UPO upo_place not found")

# S refs a1,a2 ; T ref b1 ; decoy d (=wA) ; member m ; reps r*
NW = "(((a1:1,a2:1):1,r1:1):1,((b1:1,r2:1):1,m:1):1,(d:1,r3:1):1);"
NW_WA = "(((a1:1,d:1):1,a2:1):1,((b1:1,r2:1):1,m:1):1);"   # wA inside the smallest clade of a1,a2


def write_refs(tmp):
    files = {"S": ["a1", "a2"], "T": ["b1"]}
    out = {}
    for k, ids in files.items():
        p = tmp / f"{k}.faa"
        p.write_text("".join(f">{i}\nMKV\n" for i in ids))
        out[k] = str(p)
    (tmp / "dec.faa").write_text(">d\nMKV\n")
    (tmp / "mem.faa").write_text(">m\nMKV\n")
    return out


@needs_upo
def test_edge_support_reads_labels_on_both_sides():
    upo = sp.load_upo()
    sup = sb.edge_support(upo, "(((a1:1,a2:1)0.9:1,r1:1)0.5:1,(b1:1,r2:1)0.5:1);")
    assert sup[frozenset({"a1", "a2"})] == 0.9
    assert sup[frozenset({"r1", "b1", "r2"})] == 0.9


@needs_upo
def test_edge_support_ignores_missing_and_non_numeric_labels():
    upo = sp.load_upo()
    sup = sb.edge_support(upo, "(((a1:1,a2:1):1,r1:1)x:1,(b1:1,r2:1):1);")
    assert sup == {}


@needs_upo
def test_stability_counts_target_and_wa_question(tmp_path):
    refs = write_refs(tmp_path)
    (tmp_path / "best.nwk").write_text(NW + "\n")
    (tmp_path / "bs.nwk").write_text("\n".join([NW, NW, NW_WA]) + "\n")
    s = sb.run_stability(tmp_path / "best.nwk", tmp_path / "bs.nwk", refs, tmp_path / "dec.faa", tmp_path / "mem.faa", "d", {},
                         tmp_path / "o.tsv", tmp_path / "o.json")
    assert s["n_bootstrap_trees"] == 3
    # S target exists when no decoy/member/other ref is inside: true in the first two trees, false in NW_WA (d inside)
    assert s["steps"]["S"]["trees_with_target"] == 2 and s["steps"]["S"]["target_in_best_tree"]
    assert s["steps"]["S"]["n_exact_same_target"] == 2
    assert s["steps"]["T"]["trees_with_target"] == 3
    assert (tmp_path / "o.tsv").read_text().startswith("step\tmeasure")


@needs_upo
def test_wa_measure_counts_trees_with_refs_apart_from_wa(tmp_path):
    refs = write_refs(tmp_path)
    (tmp_path / "best.nwk").write_text(NW_WA + "\n")
    (tmp_path / "bs.nwk").write_text("\n".join([NW, NW_WA]) + "\n")
    # the "heptaketide" step is S here: pass S as the step named heptaketide_pks
    r = {"heptaketide_pks": refs["S"], "T": refs["T"]}
    s = sb.run_stability(tmp_path / "best.nwk", tmp_path / "bs.nwk", r, tmp_path / "dec.faa", tmp_path / "mem.faa", "d",
                         {"wA": ["d"]}, tmp_path / "o.tsv", tmp_path / "o.json")
    # NW: smallest side of a1,a2 is {a1,a2}, no d. NW_WA: it is {a1,d,a2}, holds d.
    assert s["hepta_refs_clade_without_wA"] == {"trees": 1, "of": 2, "in_best_tree": False}
    assert s["variants"]["wA"]["trees_with_target"] == 2     # with d allowed, both trees have a side
    assert s["steps"]["heptaketide_pks"]["trees_with_target"] == 1


@needs_upo
def test_support_report_rows(tmp_path):
    refs = write_refs(tmp_path)
    nw = "(((a1:1,a2:1)0.9:1,r1:1)0.8:1,((b1:1,r2:1)0.7:1,m:1)0.6:1,(d:1,r3:1)0.5:1);"
    (tmp_path / "best.nwk").write_text(NW + "\n")
    (tmp_path / "fbp.nwk").write_text(nw + "\n")
    (tmp_path / "tbe.nwk").write_text(nw.replace("0.", "0.9") + "\n")
    rows = sb.run_support(tmp_path / "best.nwk", tmp_path / "fbp.nwk", tmp_path / "tbe.nwk", refs, tmp_path / "dec.faa",
                          tmp_path / "mem.faa", "d", {}, tmp_path / "s.tsv", tmp_path / "s.json")
    byname = {r["node"]: r for r in rows}
    assert byname["target.S"]["n_leaves"] == 3 and byname["target.S"]["fbp"] == 0.8
    assert byname["target.T"]["fbp"] == 0.7
    assert byname["all_references"]["n_leaves"] >= 4


def test_choose_subset_is_reproducible_and_sized():
    ids = [f"r{i}" for i in range(100)]
    a = sb.choose_subset(ids, 0.5, 1)
    assert a == sb.choose_subset(ids, 0.5, 1) and len(a) == 50
    assert a != sb.choose_subset(ids, 0.5, 2)
    assert len(sb.choose_subset(ids, 0.0, 1)) == 1


@needs_upo
def test_subset_keeps_labelled_leaves_and_prunes_msa(tmp_path):
    (tmp_path / "t.nwk").write_text(NW + "\n")
    (tmp_path / "m.faa").write_text("".join(f">{i}\nMKV\n" for i in ["a1", "a2", "b1", "d", "m", "r1", "r2", "r3"]))
    (tmp_path / "reps.txt").write_text("r1\nr2\nr3\n")
    info = sb.run_subset(tmp_path / "t.nwk", tmp_path / "m.faa", tmp_path / "reps.txt", 0.34, 1, tmp_path / "o.nwk", tmp_path / "o.faa",
                         tmp_path / "o.txt")
    assert info["representatives_kept"] == 1 and info["msa_sequences"] == 6
    upo = sp.load_upo()
    leaves = {l.name for l in upo.parse_newick((tmp_path / "o.nwk").read_text()).leaves()}
    kept = (tmp_path / "o.txt").read_text().split()
    assert leaves == {"a1", "a2", "b1", "d", "m"} | set(kept)


def test_compare_calls_directions():
    old = {"p1": "target", "p2": "target", "p3": "decoy", "p4": "unassigned", "p5": "decoy"}
    new = {"p1": "target", "p2": "decoy", "p3": "target", "p4": "decoy", "p6": "target"}
    c, s = sb.compare_calls(old, new)
    assert s["proteins_compared"] == 4 and s["only_in_old"] == 1 and s["only_in_new"] == 1
    assert s["target_to_not_target"] == 1 and s["not_target_to_target"] == 1 and s["changed_target_status"] == 2
    assert s["any_class_change"] == 3          # p2, p3, p4
    assert c[("unassigned", "decoy")] == 1
    assert s["fraction_changed_target_status"] == 0.5


def test_cluster_size_stats():
    s = sb.cluster_size_stats({"a": ["a", "b", "c"], "d": ["d"], "e": ["e"]})
    assert s["clusters"] == 3 and s["singletons"] == 2 and s["largest"] == 3 and s["domains"] == 5


def test_rep_coverage_applies_identity_and_coverage_rule():
    cl = {"r": ["r", "x", "y", "z"], "s": ["s", "w"]}
    hits = [("r", "r", 1.0, 1.0, 1.0), ("x", "r", 0.9, 0.9, 0.9), ("y", "r", 0.7, 0.9, 0.9), ("z", "r", 0.9, 0.5, 0.9),
            ("z", "s", 0.85, 0.9, 0.9), ("s", "s", 1.0, 1.0, 1.0)]     # w has no hit
    c = sb.rep_coverage(cl, hits, 0.8, 0.8)
    assert c["domains"] == 6 and c["non_representative_domains"] == 4
    assert c["pass_rule_assigned_rep"] == 3            # r, x, s
    assert c["non_rep_pass_rule_assigned_rep"] == 1    # x
    assert c["pass_rule_any_rep"] == 4                 # r, x, s and z (through s)
    assert c["non_rep_below_min_id_assigned"] == 1     # y at 0.7
    assert c["non_rep_without_hit_to_assigned_rep"] == 1   # w


@needs_upo
def test_stability_per_step_scope_changes_target_existence(tmp_path):
    refs = write_refs(tmp_path)
    (tmp_path / "decS.faa").write_text("")          # step S has no decoy of its own
    (tmp_path / "decT.faa").write_text(">d\nMKV\n")
    (tmp_path / "best.nwk").write_text(NW_WA + "\n")
    (tmp_path / "bs.nwk").write_text(NW_WA + "\n")
    sd = {"S": str(tmp_path / "decS.faa"), "T": str(tmp_path / "decT.faa")}
    u = sb.run_stability(tmp_path / "best.nwk", tmp_path / "bs.nwk", refs, tmp_path / "dec.faa", tmp_path / "mem.faa", "d", {},
                         tmp_path / "u.tsv", tmp_path / "u.json")
    p = sb.run_stability(tmp_path / "best.nwk", tmp_path / "bs.nwk", refs, tmp_path / "dec.faa", tmp_path / "mem.faa", "d", {},
                         tmp_path / "p.tsv", tmp_path / "p.json", sd)
    assert u["steps"]["S"]["trees_with_target"] == 0 and p["steps"]["S"]["trees_with_target"] == 1
    assert p["decoy_scope"] == "per_step" and u["decoy_scope"] == "union"


@needs_upo
def test_target_boundary_lists_nearest_outside_labelled_leaves():
    upo = sp.load_upo()
    nw = "(((a1:1,a2:1):1,r1:1):1,(((b1:1,r2:1):1,m:5):1,(d:1,r3:1):1):1);"
    roles = {"a1": "reference", "a2": "reference", "r1": "representative", "b1": "reference", "r2": "representative", "m": "named_member",
             "d": "decoy", "r3": "representative"}
    res = sb.target_boundary(upo, nw, {"a1", "a2", "r1"}, roles, k=3)
    assert [(r[0], r[1], r[2]) for r in res] == [("d", "decoy", 5.0), ("b1", "reference", 6.0), ("m", "named_member", 9.0)]
    assert [r[0] for r in sb.target_boundary(upo, nw, {"a1", "a2", "r1"}, roles, k=1)] == ["d"]
