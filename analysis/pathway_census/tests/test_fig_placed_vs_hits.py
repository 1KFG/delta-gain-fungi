"""Tests for reports/carotenoid/figures/bin/fig_placed_vs_hits.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "fig_placed_vs_hits", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_placed_vs_hits.py")
fp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fp)


def write(p):
    p.write_text("class\tn_genomes\tshare_rule_hits\tshare_low_hits\tshare_placed\nA\t500\t0.9\t0.99\t0.98\nB\t100\t0.1\t0.9\t0.8\nALL\t600\t0.8\t0.97\t0.95\n")


def test_read_by_class_filters_by_n_and_drops_all(tmp_path):
    p = tmp_path / "b.tsv"
    write(p)
    rows = fp.read_by_class(p, 150)
    assert [r[0] for r in rows] == ["A"] and rows[0][2:] == (0.9, 0.99, 0.98)


def test_main_writes_png_and_pdf(tmp_path):
    p = tmp_path / "b.tsv"
    write(p)
    out = tmp_path / "o"
    assert fp.main(["--by-class", str(p), "--out-dir", str(out), "--min-n", "50", "--name", "x"]) == 0
    assert (out / "x.pdf").read_bytes().startswith(b"%PDF") and (out / "x.png").stat().st_size > 1000
