"""Tests for bin/pks_joint_candidates.py and bin/step5_joint.py. Synthetic inputs only."""
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

from pathway_census import hits

BIN = Path(__file__).parents[1] / "bin"


def load(name):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


pj = load("pks_joint_candidates")
sj = load("step5_joint")
sp = sj.sp
needs_zstd = pytest.mark.skipif(not hits.tool_available("zstd"), reason="zstd not on PATH")
needs_upo = pytest.mark.skipif(not sp.UPO_BIN.exists(), reason="enzymes_UPO upo_place not found")
HDR = "protein_id\thmm_name\tscore\tevalue\tali_from\tali_to\thmm_from\thmm_to\tstep_id"


def row(p, s, step):
    return f"{p}\t{step}\t{s}\t1e-5\t1\t100\t1\t100\t{step}"


def test_union_rows_best_score_wins_and_relabels():
    rows = [row("a", 300, "x"), row("a", 500, "y"), row("b", 100, "x"), row("c", 50, "z")]
    kept, uni, per, best = pj.union_rows(rows, HDR.split("\t"), ["x", "y"], "joint")
    assert len(kept) == 3                       # the row of step z is dropped
    assert {l.split("\t")[0]: l.split("\t")[-1] for l in uni} == {"a": "joint", "b": "joint"}
    assert best == {"a": "y", "b": "x"}
    assert per["a"] == {"x": 300.0, "y": 500.0}


def test_union_rows_tie_goes_to_first_listed_step():
    rows = [row("a", 300, "y"), row("a", 300, "x")]
    _, _, _, best = pj.union_rows(rows, HDR.split("\t"), ["x", "y"], "joint")
    assert best["a"] == "x"


@needs_zstd
def test_candidates_end_to_end(tmp_path):
    ind = tmp_path / "in"
    ind.mkdir()
    pj.write_zst(ind / "candidates.chunk00000.PROVISIONAL.tsv.zst",
                 ["# PROVISIONAL", HDR, row("a_1", 300, "x"), row("a_1", 500, "y"), row("b_1", 100, "x"), row("c_1", 9, "other")])
    pj.write_zst(ind / "candidates.chunk00001.PROVISIONAL.tsv.zst", ["# PROVISIONAL", HDR, row("d_1", 70, "y")])
    rc = pj.main(["--in-dir", str(ind), "--steps", "x,y", "--joint-id", "j", "--steps-dir", str(tmp_path / "s"),
                  "--union-dir", str(tmp_path / "u"), "--scores-out", str(tmp_path / "sc.tsv")])
    assert rc == 0
    summ = json.loads((tmp_path / "u" / "union_summary.json").read_text())
    assert summ["proteins_union"] == 3 and summ["rows_per_step"] == {"x": 2, "y": 2}
    assert summ["proteins_by_n_steps_hit"] == {"1": 2, "2": 1}
    sc = (tmp_path / "sc.tsv").read_text().splitlines()
    assert sc[0].split("\t")[:3] == ["protein_id", "score_x", "score_y"]
    assert sc[1].split("\t")[:6] == ["a_1", "300.0", "500.0", "2", "y", "candidates.chunk00000.PROVISIONAL.tsv.zst"]
    assert (tmp_path / "u" / "candidates.chunk00001.PROVISIONAL.tsv.zst").exists()


def sides_of(newick):
    upo = sp.load_upo()
    top = upo.parse_newick(newick)
    for i, n in enumerate(sp._nodes(top)[1:]):
        n.edge = i
    return sp.all_sides(upo.edge_bipartitions(top))


@needs_upo
def test_joint_targets_are_disjoint_and_bounded():
    # S refs a1,a2 ; T ref b1 ; decoy d ; member m ; representatives r*
    nw = "(((a1:1,a2:1):1,r1:1):1,(((b1:1,r2:1):1,m:1):1,(d:1,r3:1):1):1);"
    sides = sides_of(nw)
    t = sj.define_joint_targets(sides, {"S": {"a1", "a2"}, "T": {"b1"}}, {"d"}, set(), {"m"})
    assert t["S"] == frozenset({"a1", "a2", "r1"})        # largest side with a1,a2 and no excluded leaf
    assert t["T"] == frozenset({"b1", "r2"})
    assert not (t["S"] & t["T"])


@needs_upo
def test_joint_target_is_none_when_references_are_interleaved():
    # a1 and a2 are separated by b1: no side holds a1,a2 without b1
    nw = "((a1:1,b1:1):1,(a2:1,d:1):1);"
    t = sj.define_joint_targets(sides_of(nw), {"S": {"a1", "a2"}, "T": {"b1"}}, {"d"}, set(), set())
    assert t["S"] is None
    assert t["T"] == frozenset({"b1"})


@needs_upo
def test_run_clades_writes_per_step_files_and_report(tmp_path):
    nw = "(((a1:1,a2:1):1,r1:1):1,(((b1:1,r2:1):1,m:1):1,(d:1,r3:1):1):1);"
    (tmp_path / "t.nwk").write_text(nw + "\n")
    for name, ids in {"S.faa": ["a1", "a2"], "T.faa": ["b1"], "dec.faa": ["d"], "mem.faa": ["m"]}.items():
        (tmp_path / name).write_text("".join(f">{i}\nMKV\n" for i in ids))
    (tmp_path / "reps.txt").write_text("r1\nr2\nr3\n")
    rc = sj.run_clades(str(tmp_path / "t.nwk"), {"S": str(tmp_path / "S.faa"), "T": str(tmp_path / "T.faa")},
                       str(tmp_path / "dec.faa"), str(tmp_path / "mem.faa"), None, str(tmp_path / "reps.txt"), tmp_path, "j")
    assert rc == 0
    rep = json.loads((tmp_path / "clade_report.json").read_text())
    assert rep["per_step"]["S"]["target_size"] == 3 and rep["per_step"]["T"]["target_size"] == 2
    assert rep["per_step"]["S"]["target_representatives"] == 1
    assert rep["target_pairs_sharing_a_leaf"] == 0
    rows = sp.read_tsv(tmp_path / "clade_members.S.tsv")
    lab = {r["label"]: r["our_clade"] for r in rows}
    assert lab["r1"] == "target" and lab["r2"] == "decoy" and lab["m"] == ""   # named member gets no clade label
    assert (tmp_path / "clade_members.tsv").exists() and (tmp_path / "ref.rooted.nwk").exists()


@needs_upo
def test_run_clades_exit_3_when_no_step_has_a_target(tmp_path):
    (tmp_path / "t.nwk").write_text("((a1:1,b1:1):1,(a2:1,d:1):1);\n")
    for name, ids in {"S.faa": ["a1", "a2"], "dec.faa": ["d"], "T.faa": ["b1", "b2"]}.items():
        (tmp_path / name).write_text("".join(f">{i}\nMKV\n" for i in ids))
    (tmp_path / "reps.txt").write_text("")
    (tmp_path / "t.nwk").write_text("((a1:1,b1:1):1,((a2:1,b2:1):1,d:1):1);\n")
    rc = sj.run_clades(str(tmp_path / "t.nwk"), {"S": str(tmp_path / "S.faa"), "T": str(tmp_path / "T.faa")},
                       str(tmp_path / "dec.faa"), None, None, str(tmp_path / "reps.txt"), tmp_path, "j")
    assert rc == 3
    rep = json.loads((tmp_path / "clade_report.json").read_text())
    assert rep["per_step"]["S"]["refs_isolated"] is False and rep["per_step"]["T"]["refs_isolated"] is False
    assert not (tmp_path / "clade_members.S.tsv").exists()


@needs_zstd
def test_filter_placements_keeps_only_proteins_with_a_row_of_the_step(tmp_path):
    cd = tmp_path / "c"
    cd.mkdir()
    pj.write_zst(cd / "candidates.chunk00000.PROVISIONAL.tsv.zst",
                 ["# PROVISIONAL", HDR, row("a_1", 300, "x"), row("b_1", 100, "y")])
    assign = tmp_path / "assign.tsv"
    cols = ["query", "our_clade", "our_best", "our_lwr", "our_stem_lwr", "pendant_length", "best_edge_lwr", "flags", "fig1_clade", "fig1_lwr"]
    rows = [["a_1", "target", "target", "0.9", "0.9", "0.1", "0.5", "", "", "0"],
            ["b_1", "decoy", "decoy", "0.8", "0.8", "0.1", "0.5", "", "", "0"]]
    assign.write_text("\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    out = sj.run_filter_placements(str(assign), "x", str(cd), tmp_path / "p.tsv", tmp_path / "d.tsv", tmp_path / "c.tsv")
    assert out["candidate_proteins"] == 1 and out["placements"] == 1 and out["queries_missing_from_assign"] == 0
    p = sp.read_tsv(tmp_path / "p.tsv")
    assert p == [{"protein_id": "a_1", "step_id": "x", "clade": "target", "support": "0.9"}]


@needs_upo
def test_include_variant_widens_target_past_a_decoy():
    nw = "(((a1:1,a2:1):1,d:1):1,((b1:1,r2:1):1,m:1):1);"
    sides = sides_of(nw)
    plain = sj.define_joint_targets(sides, {"S": {"a1", "a2"}, "T": {"b1"}}, {"d"}, set(), {"m"})
    wide = sj.define_joint_targets(sides, {"S": {"a1", "a2"}, "T": {"b1"}}, {"d"}, set(), {"m"}, {"S": {"d"}})
    assert plain["S"] == frozenset({"a1", "a2"})
    assert wide["S"] == frozenset({"a1", "a2", "d"})
    assert wide["T"] == plain["T"]
