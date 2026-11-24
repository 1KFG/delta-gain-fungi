"""Tests for bin/filter_candidate_folder.py and bin/screen_outgroups.py."""
import importlib.util
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).parents[1] / "bin"


def load(name):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


fc = load("filter_candidate_folder")
so = load("screen_outgroups")
HDR = "protein_id\thmm_name\tscore\tevalue\tali_from\tali_to\thmm_from\thmm_to\tstep_id"


def row(p, s, step):
    return f"{p}\t{step}\t{s}\t1e-5\t1\t100\t1\t100\t{step}"


def test_filter_lines_keeps_step_and_floor():
    lines = ["# PROVISIONAL", HDR, row("a", 300, "x"), row("b", 99.9, "x"), row("c", 500, "y"), row("d", 100, "x")]
    header, kept, n_step, n_all = fc.filter_lines(lines, "x", 100.0)
    assert header == HDR
    assert [k.split("\t")[0] for k in kept] == ["a", "d"]      # floor is inclusive
    assert (n_step, n_all) == (3, 4)


def test_filter_lines_stops_on_short_row():
    with pytest.raises(SystemExit):
        fc.filter_lines([HDR, "a\tb"], "x", 1.0)


def test_filter_lines_stops_without_score_column():
    with pytest.raises(SystemExit):
        fc.filter_lines(["protein_id\tstep_id", "a\tx"], "x", 1.0)


def test_run_writes_all_chunk_files_with_header(tmp_path):
    src = tmp_path / "in"
    src.mkdir()
    for i, rows in enumerate([[row("a", 300, "x"), row("b", 50, "x")], [row("c", 10, "y")]]):
        text = "# PROVISIONAL\n" + HDR + "\n" + "".join(r + "\n" for r in rows)
        subprocess.run(["zstd", "-q", "-f", "-o", str(src / f"candidates.chunk0000{i}.PROVISIONAL.tsv.zst")],
                       input=text.encode(), check=True)
    out = tmp_path / "out"
    s = fc.run(src, out, "x", 100.0)
    assert s["rows_kept"] == 1 and s["rows_of_step"] == 2
    files = sorted(p.name for p in out.glob("candidates.chunk*.zst"))
    assert files == ["candidates.chunk00000.PROVISIONAL.tsv.zst", "candidates.chunk00001.PROVISIONAL.tsv.zst"]
    lines = fc.read_zst_lines(out / files[1])
    assert lines[0] == "# PROVISIONAL" and lines[1] == HDR and len(lines) == 2   # empty chunk keeps header
    assert json.loads((out / "filter_summary.json").read_text())["floor"] == 100.0


def test_screen_rows_best_and_union():
    df = pd.DataFrame({"protein_id": ["p", "p", "p", "q"], "score": [40.0, 30.0, 5.0, 8.0],
                       "hmm_from": [1, 51, 91, 1], "hmm_to": [40, 100, 100, 100]})
    recs = [("p", "p sp|P|X_Y d", "M"), ("q", "q", "M"), ("r", "r", "M")]
    rows = {r["protein"]: r for r in so.screen_rows(df, "h", 100, recs)}
    # best domain 40 bits covers 40%; domains >= 10 bits cover 1-40 and 51-100 = 90% of the HMM
    assert rows["p"]["best_cov"] == 0.4 and rows["p"]["pass_best"] is False
    assert rows["p"]["union_cov"] == 0.9 and rows["p"]["pass_union"] is True
    assert rows["p"]["n_domains_ge_min"] == 2
    assert rows["q"]["pass_best"] is False and rows["q"]["pass_union"] is False   # 8 bits < 10
    assert rows["r"]["best_score"] == "" and rows["r"]["pass_best"] is False


# ---- bin/step5_single_ref_clades.py ----
sc = load("step5_single_ref_clades")
s5 = sc.sp
UPO_OK = s5.UPO_BIN.exists()
needs_upo = pytest.mark.skipif(not UPO_OK, reason="enzymes_UPO tools not found")


def bip_of(newick):
    upo = s5.load_upo()
    top = upo.parse_newick(newick)
    for i, n in enumerate(s5._nodes(top)[1:]):
        n.edge = i
    return upo.edge_bipartitions(top)


NW1 = "(((r:1,rep1:1):1,rep2:1):1,(m1:1,rep3:1):1,(og:1,rep4:1):1);"


@needs_upo
def test_smallest_rule_is_degenerate_for_one_reference():
    c = s5.define_clades(bip_of(NW1), {"r"}, set(), {"og"}, {"m1"})
    assert c["target"] == frozenset({"r"})


@needs_upo
def test_largest_rule_stops_at_labelled_neighbours():
    c = sc.define_single_ref_clades(bip_of(NW1), {"r"}, set(), {"og"}, {"m1"})
    assert c["target"] == frozenset({"r", "rep1", "rep2"})
    assert c["decoy"] == frozenset({"rep3", "rep4"})          # outside target, not outgroup, not member
    assert c["outgroup"] == frozenset({"og"})


@needs_upo
def test_largest_rule_with_decoy_counts_as_labelled():
    nw = "((r:1,rep1:1):1,(d:1,rep2:1):1,(og:1,rep3:1):1);"
    c = sc.define_single_ref_clades(bip_of(nw), {"r"}, {"d"}, {"og"})
    assert c["target"] == frozenset({"r", "rep1"})
    assert "d" in c["decoy"] and "rep2" in c["decoy"]


@needs_upo
def test_largest_free_side_is_unique_and_nested():
    sides = s5.all_sides(bip_of(NW1))
    free = [s for s in sides if "r" in s and not (s & {"og", "m1"})]
    big = sc.largest_free_side(sides, {"r"}, {"og", "m1"})
    assert all(s <= big for s in free)


@needs_upo
def test_run_stops_without_labelled_leaves(tmp_path):
    (tmp_path / "refs.faa").write_text(">r\nMK\n")
    (tmp_path / "reps").write_text("rep1\n")
    (tmp_path / "t.nwk").write_text("((r:1,rep1:1):1,rep2:1,rep3:1);")
    assert sc.run(tmp_path / "t.nwk", tmp_path / "refs.faa", None, None, None, tmp_path / "reps", tmp_path, "s") == 4


@needs_upo
def test_run_writes_files_and_keeps_original(tmp_path):
    pytest.importorskip("pandas")
    (tmp_path / "refs.faa").write_text(">r\nMK\n")
    (tmp_path / "og.faa").write_text(">og\nMK\n")
    (tmp_path / "mem.faa").write_text(">m1\nMK\n")
    (tmp_path / "reps").write_text("rep1\nrep2\nrep3\nrep4\n")
    (tmp_path / "t.nwk").write_text(NW1)
    (tmp_path / "clade_report.json").write_text("{}")
    (tmp_path / "clade_members.tsv").write_text("old\n")
    rc = sc.run(tmp_path / "t.nwk", tmp_path / "refs.faa", None, tmp_path / "og.faa", tmp_path / "mem.faa",
                tmp_path / "reps", tmp_path, "s")
    assert rc == 0
    rep = json.loads((tmp_path / "clade_report.json").read_text())
    assert rep["target_size"] == 3 and rep["target_representatives"] == 2 and rep["target_is_reference_only"] is False
    assert (tmp_path / "clade_members.smallest.tsv").read_text() == "old\n"
    rows = (tmp_path / "clade_members.tsv").read_text().splitlines()
    assert "rep1\ttarget\t\trepresentative" in rows and "m1\t\t\tnamed_member" in rows
    assert (tmp_path / "ref.rooted.nwk").exists()


def test_run_rejects_two_references(tmp_path):
    (tmp_path / "refs.faa").write_text(">r1\nMK\n>r2\nMK\n")
    with pytest.raises(SystemExit):
        sc.run(tmp_path / "t.nwk", tmp_path / "refs.faa", None, None, None, tmp_path / "reps", tmp_path, "s")


# ---- bin/select_by_screen.py ----
sel = load("select_by_screen")


def test_select_by_screen_modes_ids_and_exclude(tmp_path):
    scr = tmp_path / "s.tsv"
    scr.write_text("hmm\tprotein\tentry\tbest_score\tbest_cov\tn_domains_ge_min\tunion_cov\tpass_best\tpass_union\n"
                   "h\ta\ta\t50\t0.6\t1\t0.6\tTrue\tTrue\n"
                   "h\tb\tb\t50\t0.2\t3\t0.7\tFalse\tTrue\n"
                   "h\tc\tc\t5\t0.1\t0\t0.0\tFalse\tFalse\n"
                   "g\td\td\t50\t0.6\t1\t0.6\tTrue\tTrue\n")
    assert sel.passing(scr, "h", "best") == {"a"}
    assert sel.passing(scr, "h", "union") == {"a", "b"}
    recs = [("a", "a x", "M"), ("b", "b x", "M"), ("c", "c", "M"), ("d", "d", "M")]
    ok = sel.passing(scr, "h", "union")
    assert [r[0] for r in sel.select(recs, ok)] == ["a", "b"]
    assert [r[0] for r in sel.select(recs, ok, ids={"b", "c"})] == ["b"]
    assert [r[0] for r in sel.select(recs, ok, exclude=["a"])] == ["b"]


# ---- bin/hits_by_class_floors.py ----
hbc = load("hits_by_class_floors")


def test_class_table_merges_small_classes():
    genomes = ["a1", "a2", "a3", "b1"]
    classes = {"a1": "A", "a2": "A", "a3": "A", "b1": "B"}
    rows = hbc.class_table(genomes, classes, {"x": {"a1", "b1"}}, min_n=2)
    d = {r[0]: r for r in rows}
    assert d["A"][1] == 3 and abs(d["A"][2]["x"] - 1 / 3) < 1e-9
    assert d["(other, small classes)"][1] == 1 and d["(other, small classes)"][2]["x"] == 1.0


def test_hits_at_floors(tmp_path):
    text = "# PROVISIONAL\n" + HDR + "\n" + row("G1_1", 300, "x") + "\n" + row("G2_1", 50, "x") + "\n" + row("G3_1", 500, "y") + "\n"
    subprocess.run(["zstd", "-q", "-f", "-o", str(tmp_path / "candidates.chunk00000.PROVISIONAL.tsv.zst")], input=text.encode(), check=True)
    h = hbc.hits_at_floors(tmp_path, "x", [10, 100])
    assert h[10] == {"G1", "G2"} and h[100] == {"G1"}


# ---- bin/single_ref_nearest.py ----
srn = load("single_ref_nearest")


@needs_upo
def test_distances_from_and_nearest():
    upo = s5.load_upo()
    top = upo.parse_newick("((r:1,a:2):1,(m:1,b:3):2,c:5);")
    d = srn.distances_from(top, "r", ["a", "b", "c", "m", "r"])
    assert [x[0] for x in d] == ["a", "m", "b", "c"]          # r itself left out; b and c tie at 7, name order
    assert dict(d) == {"a": 3, "m": 5, "b": 7, "c": 7}
    near = srn.nearest_among(top, ["a", "b"], ["r", "m"])
    assert near[0][:2] == ("a", "r") and near[1][:2] == ("b", "m")


def test_genus_class_missing_label():
    assert srn.genus_class("X", {}) == ("(no label)", "(no label)")
    assert srn.genus_class("X", {"X": {"GENUS": "Mucor", "CLASS": "Mucoromycetes"}}) == ("Mucor", "Mucoromycetes")


# ---- bin/target_overlap.py ----
tov = load("target_overlap")


def test_overlaps_pairs_and_depth():
    sets = {"a": {"G1_1", "G2_1", "G3_1"}, "b": {"G1_1", "G2_1"}, "c": {"G2_1", "G4_1"}}
    rows, depth = tov.overlaps(sets)
    d = {(a, b): (n, g) for a, b, n, g in rows}
    assert d[("a", "b")] == (2, 2) and d[("a", "c")] == (1, 1) and d[("b", "c")] == (1, 1)
    assert depth == {1: 2, 2: 1, 3: 1}      # G3, G4 in one set; G1 in two; G2 in three


def test_target_set_reads_only_target_calls(tmp_path):
    (tmp_path / "placements_detail.tsv").write_text("protein_id\tclade_call\nG1_1\ttarget\nG2_1\tdecoy\nG3_1\tunassigned\n")
    assert tov.target_set(tmp_path) == {"G1_1"}
