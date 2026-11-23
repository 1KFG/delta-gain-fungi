"""Tests for bin/filter_calls_busco.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("filter_calls_busco", Path(__file__).parents[1] / "bin" / "filter_calls_busco.py")
fb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fb)


def calls():
    mk = lambda g, cls, s, a, b: dict(locustag=g, CLASS=cls, searched=s, s1=a, s2=b, s1_source="placed", s2_source="hits")
    return [mk("A", "C1", "1", "1", "0"), mk("B", "C1", "1", "0", "0"), mk("C", "", "1", "1", "1"),
            mk("D", "C2", "0", "NA", "NA"), mk("E", "C2", "1", "1", "1")]


def test_filter_sets_all_steps_na_below_threshold_and_keeps_boundary():
    busco = {"A": 95.0, "B": 79.9, "C": 80.0, "D": 10.0, "E": None}
    rows, summ, nb = fb.filter_calls(calls(), busco, 80.0)
    by = {r["locustag"]: r for r in rows}
    assert (by["A"]["s1"], by["A"]["s2"]) == ("1", "0")
    assert (by["B"]["s1"], by["B"]["s2"]) == ("NA", "NA") and by["B"]["busco_filtered"] == "1"
    assert by["C"]["s1"] == "1"  # exactly 80 stays
    assert by["D"]["busco_filtered"] == "0"  # already NA: not counted as removed
    assert by["E"]["s1"] == "NA" and nb == 1  # no busco value is removed and counted
    d = {s["class"]: s for s in summ}
    assert d["C1"]["n_tips_with_data"] == 2 and d["C1"]["n_removed"] == 1
    assert d["(no class)"]["n_removed"] == 0 and d["C2"]["n_tips_with_data"] == 1 and d["C2"]["n_removed"] == 1


def test_parse_busco_handles_empty_and_text():
    assert fb.parse_busco("") is None and fb.parse_busco("x") is None and fb.parse_busco("99.1") == 99.1


def test_relative_rule_median_drop_and_absolute_floor():
    reps = [dict(LOCUSTAG=f"A{i}", CLASS="C1", busco_complete=str(v)) for i, v in enumerate([90, 90, 90, 90, 90, 90, 90, 74, 80, 40])]
    reps += [dict(LOCUSTAG="S1", CLASS="Small", busco_complete="60"), dict(LOCUSTAG="S2", CLASS="Small", busco_complete="45"),
             dict(LOCUSTAG="N1", CLASS="", busco_complete="49"), dict(LOCUSTAG="N2", CLASS="", busco_complete="55"), dict(LOCUSTAG="X", CLASS="C1", busco_complete="")]
    removed, floor_only = fb.relative_removed(reps, drop=15, min_class_n=8, floor=50)
    # C1 median of the 10 values with a value is 90: 74 is 16 below (removed), 80 is 10 below (kept), 40 removed (relative and floor)
    assert "A7" in removed and "A8" not in removed and "A9" in removed
    assert "S1" not in removed and "S2" in removed  # small class: only the floor applies
    assert "N1" in removed and "N2" not in removed  # no class: only the floor applies
    assert "X" in removed  # no busco value
    assert floor_only == {"S2", "N1", "X"}
