"""Tests for reports/carotenoid/figures/bin/fig_census_search.py (data preparation and figure output)."""
import importlib.util
import subprocess
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "fig_census_search", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_census_search.py")
fc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fc)

needs_zstd = pytest.mark.skipif(subprocess.run(["which", "zstd"], capture_output=True).returncode != 0, reason="zstd not on PATH")


def write_inputs(tmp_path):
    s2 = tmp_path / "step2"
    s2.mkdir()
    (s2 / "step_scores.PROVISIONAL.tsv").write_text(
        "step\tn_refs\tref_min\tref_max\tloo_min\tmax_decoy_score\n"
        "a\t5\t300\t350\t150\t20\n"
        "b\t1\t800\t800\t\t300\n")
    ms = tmp_path / "ms.tsv"
    ms.write_text("step\tmin_score\tbasis\tn_refs\trule_min_score\na\t112.5\tx\t5\t150\nb\t300\tx\t1\t300\n")
    return s2, ms


def write_candidates(d, rows):
    d.mkdir(exist_ok=True)
    txt = "# PROVISIONAL\nprotein_id\thmm_name\tscore\tevalue\tali_from\tali_to\thmm_from\thmm_to\tstep_id\n"
    for p, s, sc in rows:
        txt += f"{p}\th\t{sc}\t1e-5\t1\t10\t1\t10\t{s}\n"
    src = d / "c.tsv"
    src.write_text(txt)
    subprocess.run(["zstd", "-q", "-f", "-o", str(d / "candidates.chunk00000.PROVISIONAL.tsv.zst"), str(src)], check=True)


def test_step_table_joins_rule_and_used_values(tmp_path):
    s2, ms = write_inputs(tmp_path)
    t = {r["step"]: r for r in fc.step_table(s2 / "step_scores.PROVISIONAL.tsv", ms)}
    assert t["a"]["used"] == 112.5 and t["a"]["rule"] == 150 and t["a"]["loo_min"] == 150
    assert t["b"]["loo_min"] is None and t["b"]["max_decoy"] == 300 and t["b"]["n_refs"] == 1


def test_num_handles_blank_and_text():
    assert fc.num("") is None and fc.num(None) is None and fc.num("x") is None and fc.num("2.5") == 2.5


@needs_zstd
def test_count_candidates_counts_rows_and_distinct_genomes(tmp_path):
    d = tmp_path / "s"
    write_candidates(d, [("G1_1", "a", 10), ("G1_2", "a", 11), ("G2_1", "a", 12), ("G1_1", "b", 5)])
    assert fc.count_candidates(d) == {"a": (3, 2), "b": (1, 1)}


@needs_zstd
def test_scores_for_filters_by_step(tmp_path):
    d = tmp_path / "s"
    write_candidates(d, [("G1_1", "a", 10), ("G1_2", "b", 99)])
    assert fc.scores_for(d, "b") == [99.0]


def test_read_tsv_skips_comment_lines(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("# note\na\tb\n1\t2\n")
    assert fc.read_tsv(p) == [{"a": "1", "b": "2"}]


@needs_zstd
def test_main_writes_png_and_pdf_for_each_figure(tmp_path):
    s2, ms = write_inputs(tmp_path)
    write_candidates(tmp_path / "search", [("G1_1", "a", 10), ("G2_1", "b", 12)])
    write_candidates(tmp_path / "side", [("G1_1", "a", 30), ("G1_2", "a", 60), ("G2_1", "b", 400)])
    out = tmp_path / "out"
    rc = fc.main(["--step2-dir", str(s2), "--search-dir", str(tmp_path / "search"), "--min-scores", str(ms),
                  "--side-dir", str(tmp_path / "side"), "--out-dir", str(out)])
    assert rc == 0
    names = sorted(p.name for p in out.iterdir())
    for stem in ("fig1_step2_separation", "fig2_candidates_per_step", "fig3_low_cutoff_side_run"):
        assert f"{stem}.png" in names and f"{stem}.pdf" in names
    assert (out / "fig1_step2_separation.pdf").read_bytes().startswith(b"%PDF")
