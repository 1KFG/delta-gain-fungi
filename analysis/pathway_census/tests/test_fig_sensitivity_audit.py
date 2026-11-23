"""Tests for reports/carotenoid/figures/bin/fig_sensitivity_audit.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "fig_sensitivity_audit", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_sensitivity_audit.py")
fs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fs)


def test_split_steps_by_entry_count():
    summ = [dict(step="a", n_entries="10"), dict(step="b", n_entries="2"), dict(step="c", n_entries="0")]
    d, f = fs.split_steps(summ, 5)
    assert [s["step"] for s in d] == ["a"] and [s["step"] for s in f] == ["b", "c"]


def test_draw_writes_png_and_pdf_and_reports_which_steps_were_drawn(tmp_path):
    entries = [dict(step="a", accession=f"P{i}", name="n", organism="o", score=str(100 + 50 * i), at_cutoff=str(i > 1), at_rule="True") for i in range(6)]
    entries.append(dict(step="a", accession="Pnone", name="n", organism="o", score="", at_cutoff="False", at_rule="False"))
    summ = [dict(step="a", ec="1.1.1.1", cutoff="200", n_entries="7", n_at_cutoff="4", recall_at_cutoff="0.571"),
            dict(step="b", ec="2.2.2.2", cutoff="50", n_entries="1", n_at_cutoff="1", recall_at_cutoff="1.0")]
    drawn, few = fs.draw(entries, summ, tmp_path, "x")
    assert drawn == ["a"] and few == ["b"]
    assert (tmp_path / "x.pdf").read_bytes().startswith(b"%PDF") and (tmp_path / "x.png").stat().st_size > 1000


def test_num_handles_blank():
    assert fs.num("") is None and fs.num("3.5") == 3.5
