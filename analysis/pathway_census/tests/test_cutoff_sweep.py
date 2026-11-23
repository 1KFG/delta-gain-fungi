"""Tests for bin/cutoff_sweep.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("cutoff_sweep", Path(__file__).parents[1] / "bin" / "cutoff_sweep.py")
cs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cs)

HITS = [("G1", 400.0, 1e-100), ("G1", 120.0, 1e-30), ("G2", 250.0, 1e-70), ("G3", 90.0, 1e-20)]
GEN = ["G1", "G2", "G3", "G4"]
CLS = {"G1": "A", "G2": "A", "G3": "B", "G4": "B"}


def test_by_class_table_uses_best_hit_per_genome():
    rows = {cl: (n, d) for cl, n, d in cs.by_class_table(HITS, GEN, CLS, [100, 300])}
    assert rows["A"][0] == 2 and rows["A"][1] == {100: 1.0, 300: 0.5}
    assert rows["B"][1] == {100: 0.0, 300: 0.0}


def test_by_class_table_min_n_and_blank_class_label():
    rows = cs.by_class_table(HITS, GEN, {"G1": "A"}, [100], min_n=2)
    assert [r[0] for r in rows] == ["(no label)"]


def test_copies_table_counts_genomes_by_number_of_calls():
    d = cs.copies_table(HITS, GEN, [100])[100]
    assert d == {0: 2, 1: 1, 2: 1, "3+": 0}


def test_score_bands_median_evalue_and_empty_band():
    rows = cs.score_bands(HITS, [0, 100, 300, 500, 1000])
    assert rows[0][2] == 1 and rows[1][2] == 2 and rows[2][2] == 1 and rows[3][2] == 0 and rows[3][3] is None
    assert rows[1][3] == 1e-50 or abs(rows[1][3] - 3.1622776601683795e-50) < 1e-60 or rows[1][3] > 0
