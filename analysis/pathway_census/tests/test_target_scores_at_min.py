"""Tests for bin/target_scores_at_min.py."""
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("tsam", Path(__file__).parents[1] / "bin" / "target_scores_at_min.py")
tsam = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tsam)

META = {"ref|a|0|R1": ("ref", "a", "R1"), "ref|b|0|R2": ("ref", "b", "R2"),
        "decoy|a|0|D1": ("decoy", "a", "D1"), "decoy|b|0|D2": ("decoy", "b", "D2")}


def test_own_references_are_skipped_and_flag_uses_min_score():
    scores = {"ref|a|0|R1": 500.0, "ref|b|0|R2": 120.0, "decoy|a|0|D1": 100.0, "decoy|b|0|D2": 30.0}
    rows = tsam.flag_rows("a", 100.0, scores, META)
    assert [r["accession"] for r in rows] == ["R2", "D1", "D2"]
    assert [r["at_or_above_min"] for r in rows] == [True, True, False]
    assert all(r["step"] == "a" and r["min_score"] == 100.0 for r in rows)


def test_no_hits_gives_no_rows():
    assert tsam.flag_rows("a", 10.0, {}, META) == []


def test_read_fasta_joins_wrapped_lines(tmp_path):
    f = tmp_path / "x.faa"
    f.write_text(">s1 desc\nAAA\nCC\n>s2\nGG\n")
    assert list(tsam.read_fasta(f)) == [("s1 desc", "AAACC"), ("s2", "GG")]


def test_main_skips_a_step_with_empty_min_score(tmp_path, monkeypatch, capsys):
    """A step with no floor (make_min_scores.py --allow-undefined) is skipped and named, not parsed as a number."""
    import shutil
    if not shutil.which("hmmsearch"):
        import pytest
        pytest.skip("hmmer not available")
    refs = tmp_path / "refs"
    refs.mkdir()
    ms = tmp_path / "min_scores.tsv"
    ms.write_text("step\tmin_score\tbasis\tn_refs\trule_min_score\na\t\tno floor defined: x\t1\t\n")
    (tmp_path / "hmm").mkdir()
    out = tmp_path / "o.tsv"
    rc = tsam.main(["--refs-dir", str(refs), "--hmm-dir", str(tmp_path / "hmm"), "--min-scores", str(ms),
                    "--out", str(out)])
    assert rc == 0
    assert "skipped, min_score empty: ['a']" in capsys.readouterr().err
