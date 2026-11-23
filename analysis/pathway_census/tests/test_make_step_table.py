"""Tests for reports/carotenoid/bin/make_step_table.py."""
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "make_step_table", Path(__file__).parents[1] / "reports" / "carotenoid" / "bin" / "make_step_table.py")
mt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mt)


def test_count_fasta_missing_file_is_zero(tmp_path):
    assert mt.count_fasta(tmp_path / "none.faa") == 0


def test_count_fasta_counts_records(tmp_path):
    p = tmp_path / "a.faa"
    p.write_text(">a\nMK\n>b\nMK\n")
    assert mt.count_fasta(p) == 2


def test_build_table_marks_lowered_cutoff_and_formats_counts(tmp_path):
    (tmp_path / "decoys").mkdir()
    (tmp_path / "decoys" / "s1.faa").write_text(">d1\nMK\n>d2\nMK\n")
    tab = [dict(step="s1", n_refs=2, rule=1000.0, used=750.0), dict(step="s2", n_refs=1, rule=20.0, used=20.0)]
    t = mt.build_table(tab, {"s1": (12345, 99)}, tmp_path)
    lines = t.splitlines()
    assert lines[2] == "| `s1` | 2 | 2 | 1000 | 750 (lowered) | 12,345 | 99 |"
    assert lines[3] == "| `s2` | 1 | 0 | 20 | 20 | 0 | 0 |"


def test_insert_replaces_between_markers_only():
    text = f"before\n{mt.BEGIN}\nOLD\n{mt.END}\nafter\n"
    out = mt.insert(text, "NEW")
    assert out == f"before\n{mt.BEGIN}\nNEW\n{mt.END}\nafter\n"


def test_insert_without_markers_is_an_error():
    with pytest.raises(ValueError):
        mt.insert("no markers", "x")
