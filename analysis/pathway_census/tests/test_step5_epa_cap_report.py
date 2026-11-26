import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
spec = importlib.util.spec_from_file_location("cap", BIN / "step5_epa_cap_report.py")
cap = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cap)


def test_leaf_counts(tmp_path):
    f = tmp_path / "a.afa"
    f.write_text(">a\nAC-\n>b\nac-\n>c\nAAA\n")
    assert cap.leaf_counts(f) == (3, 2)


def test_call_class():
    assert [cap.call_class(x) for x in ("target", "decoy", "unassigned", "")] == ["target", "other", "unassigned", "unassigned"]


def _run(tmp_path, calls, nplc):
    (tmp_path / "ref_set.trim.afa").write_text(">a\nAA\n>b\nAC\n")
    (tmp_path / "placements_detail.tsv").write_text("protein_id\tclade_call\n" + "".join(f"{p}\t{c}\n" for p, c in calls.items()))
    (tmp_path / "assign.tsv").write_text("query\tour_clade\tn_placements\n" + "".join(
        f"{p}\t{c}\t{nplc[p]}\n" for p, c in calls.items()) + "zzz\ttarget\t7\n")
    return tmp_path


def test_survey_one_restricts_to_candidates(tmp_path):
    d = _run(tmp_path, {"p1": "target", "p2": "unassigned", "p3": "decoy"}, {"p1": 1, "p2": 7, "p3": 7})
    r = cap.survey_one("x", d, None)
    assert (r["queries"], r["target"], r["other"], r["unassigned"]) == (3, 1, 1, 1)
    assert (r["cap_hit_queries"], r["cap_hit_unassigned"], r["queries_in_assign"]) == (2, 1, 3)
    assert (r["leaves"], r["distinct_leaf_seqs"]) == (2, 2)


def test_compare_classes():
    r = cap.compare_classes({"a": "target", "b": "unassigned", "c": "other", "d": "target"},
                            {"a": "target", "b": "target", "c": "target", "e": "x"})
    assert r["changed"] == 2 and r["target_other_flips"] == 1 and r["only_old"] == 1 and r["only_new"] == 1
