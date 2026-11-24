"""Tests for bin/placed_vs_hits_by_class.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("pvh", Path(__file__).parents[1] / "bin" / "placed_vs_hits_by_class.py")
pvh = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pvh)


def test_by_class_shares_and_blank_label_and_order():
    sets = {"a": {"G1", "G2"}, "b": {"G1"}}
    rows = pvh.by_class(["G1", "G2", "G3", "G4"], {"G1": "X", "G2": "X", "G3": "Y"}, sets)
    d = {c: (n, s) for c, n, s in rows}
    assert d["X"] == (2, {"a": 1.0, "b": 0.5}) and d["Y"] == (1, {"a": 0.0, "b": 0.0}) and "(no label)" in d
    assert rows[0][0] == "X"


def test_band_table_counts_calls_and_treats_missing_call_as_unassigned():
    scores = {"p1": 160.0, "p2": 210.0, "p3": 220.0, "p4": 400.0, "p5": 410.0}
    calls = {"p1": "decoy", "p2": "target", "p3": "", "p4": "target"}
    rows = pvh.band_table(scores, calls, [150, 200, 300, 500])
    assert rows[0] == (150, 200, 0, 1, 0)
    assert rows[1] == (200, 300, 1, 0, 1)
    assert rows[2] == (300, 500, 1, 0, 1)
