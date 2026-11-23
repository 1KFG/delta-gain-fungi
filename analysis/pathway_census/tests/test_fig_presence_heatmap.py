"""Tests for reports/carotenoid/figures/bin/fig_presence_heatmap.py (data preparation and output)."""
import importlib.util
import subprocess
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "fig_presence_heatmap", Path(__file__).parents[1] / "reports" / "carotenoid" / "figures" / "bin" / "fig_presence_heatmap.py")
fh = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fh)
needs_zstd = pytest.mark.skipif(subprocess.run(["which", "zstd"], capture_output=True).returncode != 0, reason="zstd not on PATH")


def rows():
    tax = lambda ph, cl: dict(PHYLUM=ph, SUBPHYLUM="", CLASS=cl, ORDER="o", FAMILY="f")
    mk = lambda g, ph, cl, a, b: dict(locustag=g, **tax(ph, cl), s1=a, s2=b)
    return [mk("G1", "P1", "C1", 1, 0), mk("G2", "P1", "C1", 1, 1), mk("G3", "P1", "C2", 0, 0), mk("G4", "P2", "C3", 1, 1)]


def test_genome_of_takes_locustag_prefix():
    assert fh.genome_of("F8027CAE_004645-T1") == "F8027CAE"


def test_presence_table_marks_calls_and_skips_genomes_without_taxonomy():
    tax = {"G1": dict(PHYLUM="P", SUBPHYLUM="", CLASS="C", ORDER="", FAMILY="")}
    out = fh.presence_table(["G1", "GX"], tax, {"s1": {"G1"}}, ["s1", "s2"])
    assert len(out) == 1 and out[0]["s1"] == 1 and out[0]["s2"] == 0


def test_summarise_fractions_n_filter_and_core_complete():
    s = fh.summarise(rows(), "CLASS", ["s1", "s2"], min_n=1, core=("s1", "s2"))
    by = {d["group"]: d for d in s}
    assert by["C1"]["n"] == 2 and by["C1"]["frac"] == {"s1": 1.0, "s2": 0.5} and by["C1"]["core_complete"] == 0.5
    assert by["C2"]["core_complete"] == 0.0
    s2 = fh.summarise(rows(), "CLASS", ["s1", "s2"], min_n=2)
    assert [d["group"] for d in s2] == ["C1"]


def test_summarise_orders_by_phylum_size_then_lineage():
    s = fh.summarise(rows(), "CLASS", ["s1", "s2"], min_n=1)
    assert [d["group"] for d in s] == ["C1", "C2", "C3"]       # P1 has 3 genomes, P2 has 1


def test_summarise_blank_label_gets_a_visible_name():
    r = rows()
    r[0]["CLASS"] = ""
    assert "(no label)" in {d["group"] for d in fh.summarise(r, "CLASS", ["s1", "s2"], 1)}


def test_floor_steps_reads_basis(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("step\tmin_score\tbasis\tn_refs\ta\t1\thighest decoy score (single reference, floor)\t1\nb\t5\tlowest leave-one-out score\t2\n".replace("step\tmin_score\tbasis\tn_refs\ta", "step\tmin_score\tbasis\tn_refs\na"))
    assert fh.floor_steps(p) == {"a"}


@needs_zstd
def test_placed_genomes_counts_only_assigned(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text("# PROVISIONAL\nprotein_id\tstep_id\tsupport\tlabel\tfusion_candidate\nG1_1\ts\t1.0\tassigned\tFalse\nG2_1\ts\t1.0\tunassigned\tFalse\n")
    assert fh.placed_genomes(p) == {"G1"}


def test_draw_writes_png_and_pdf(tmp_path):
    summ = fh.summarise(rows(), "CLASS", ["s1", "s2"], 1)
    fh.draw(summ, ["s1"], ["s2"], tmp_path, "x", "T", "note")
    assert (tmp_path / "x.pdf").read_bytes().startswith(b"%PDF") and (tmp_path / "x.png").stat().st_size > 1000


def test_core_steps_exclude_ggpp_synthase():
    assert "ggpp_synthase" not in fh.CORE_STEPS and len(fh.CORE_STEPS) == 3
