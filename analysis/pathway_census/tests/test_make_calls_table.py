"""Tests for bin/make_calls_table.py."""
import importlib.util
import subprocess
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("make_calls_table", Path(__file__).parents[1] / "bin" / "make_calls_table.py")
mc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mc)


def test_genome_of():
    assert mc.genome_of("F00043C3_004031-T1") == "F00043C3"


def test_searched_requires_all_rows_searched():
    rows = [dict(locustag="A", status="searched"), dict(locustag="B", status="failed"),
            dict(locustag="C", status="searched"), dict(locustag="C", status="missing_file"),
            dict(locustag="D", status="searched"), dict(locustag="D", status="searched")]
    assert mc.searched_genomes(rows) == {"A", "D"}


def test_placed_only_assigned_label_and_duplicates():
    rows = [dict(protein_id="A_1-T1", label="assigned"), dict(protein_id="A_2-T1", label="assigned"),
            dict(protein_id="B_1-T1", label="unassigned")]
    assert mc.placed_genomes(rows) == {"A"}


def test_hit_genomes_filters_step():
    f1 = [dict(protein_id="A_1", step_id="s1"), dict(protein_id="B_1", step_id="s2")]
    f2 = [dict(protein_id="A_9", step_id="s1")]
    assert mc.hit_genomes([f1, f2], "s1") == {"A"}
    assert mc.hit_genomes([f1, f2], "s3") == set()


def reps():
    mk = lambda g, p="/x.fa": dict(LOCUSTAG=g, pep_path=p, PHYLUM="P", CLASS="C", GENUS=g)
    return [mk("A"), mk("B"), mk("C"), mk("D", ""), mk("A")]  # A duplicated, D has no protein file


def test_build_calls_na_absent_present_and_edge_cases():
    searched = {"A", "B", "D"}  # C missing from the searched files; D has no protein file
    rows, summ = mc.build_calls(reps(), searched, {"s1": "placed", "s2": "hits"},
                                {"s1": {"A", "C", "ZZ"}, "s2": {"B"}})
    by = {r["locustag"]: r for r in rows}
    assert len(rows) == 4  # duplicate representative dropped
    assert (by["A"]["s1"], by["A"]["s2"]) == ("1", "0")
    assert (by["B"]["s1"], by["B"]["s2"]) == ("0", "1")
    assert by["C"]["s1"] == "NA" and by["C"]["s2"] == "NA"  # unsearched genome with a call stays NA
    assert by["D"]["s1"] == "NA"  # no protein file
    s1 = {s["step"]: s for s in summ}["s1"]
    assert (s1["n_present"], s1["n_absent"], s1["n_NA"]) == (1, 1, 2)
    assert s1["n_call_genomes_not_in_reps"] == 1  # ZZ
    assert s1["n_rep_calls_but_NA"] == 1  # C


@pytest.mark.skipif(subprocess.run(["which", "zstd"], capture_output=True).returncode != 0, reason="zstd not on PATH")
def test_read_tsv_zst_skips_comment(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("# PROVISIONAL\nprotein_id\tstep_id\nA_1\ts1\n")
    subprocess.run(["zstd", "-q", str(p)], check=True)
    assert mc.read_tsv(str(p) + ".zst") == [dict(protein_id="A_1", step_id="s1")]
