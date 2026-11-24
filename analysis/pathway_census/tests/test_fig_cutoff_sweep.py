"""Tests for reports/carotenoid/figures/bin/fig_cutoff_sweep.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "fig_cutoff_sweep", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_cutoff_sweep.py")
fc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fc)


def write(p):
    p.write_text("class\tn_genomes\tfrac_ge_300\tfrac_ge_100\nSaccharomycetes\t10\t0.1\t0.9\nPichiomycetes\t20\t0.0\t0.8\nAgaricomycetes\t5\t0.7\t0.9\nOther\t7\t0.9\t1.0\nALL\t42\t0.4\t0.9\n")


def test_read_sweep_sorts_cutoffs_ascending_and_keeps_alignment(tmp_path):
    p = tmp_path / "s.tsv"
    write(p)
    cut, data, ns = fc.read_sweep(p)
    assert cut == [100.0, 300.0]
    assert data["Saccharomycetes"] == [0.9, 0.1] and ns["ALL"] == 42


def test_main_writes_png_and_pdf(tmp_path):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    write(a)
    write(b)
    out = tmp_path / "o"
    assert fc.main(["--current", str(a), "--wide", str(b), "--out-dir", str(out), "--name", "x"]) == 0
    assert (out / "x.pdf").read_bytes().startswith(b"%PDF") and (out / "x.png").stat().st_size > 1000
