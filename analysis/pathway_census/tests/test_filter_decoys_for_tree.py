"""Tests for bin/filter_decoys_for_tree.py."""
import importlib.util
from pathlib import Path

import pytest

from pathway_census import hits

SPEC = importlib.util.spec_from_file_location("filter_decoys_for_tree", Path(__file__).parents[1] / "bin" / "filter_decoys_for_tree.py")
fd = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fd)
DATA = Path(__file__).parent / "data"
needs_hmmer = pytest.mark.skipif(not all(hits.tool_available(t) for t in ["hmmbuild", "hmmsearch", "mafft"]),
                                 reason="hmmer/mafft not loadable")


def test_read_fasta_keeps_header_and_first_token(tmp_path):
    p = tmp_path / "d.faa"
    p.write_text(">A desc one\nMK\nVV\n>B\nMM\n")
    assert list(fd.read_fasta(p)) == [("A", "A desc one", "MKVV"), ("B", "B", "MM")]


def test_select_keeps_at_or_above_score_and_coverage():
    recs = [("a", "a", "M"), ("b", "b", "M"), ("c", "c", "M"), ("d", "d", "M"), ("e", "e", "M")]
    scores = {"a": (30.0, 1, 100), "b": (10.0, 1, 50), "c": (9.9, 1, 100), "d": (50.0, 1, 40)}
    kept, table = fd.select(recs, scores, 10.0, 0.5, 100)
    assert [k[0] for k in kept] == ["a", "b"]
    assert table[0] == ("a", 30.0, 1.0, True)
    assert table[1] == ("b", 10.0, 0.5, True)
    assert table[2][3] is False       # score below the limit
    assert table[3] == ("d", 50.0, 0.4, False)   # high score but short coverage is dropped
    assert table[4] == ("e", None, None, False)  # no hit


def test_select_with_no_scores_keeps_nothing():
    kept, table = fd.select([("a", "a", "M")], {}, 5.0, 0.5, 100)
    assert kept == [] and table == [("a", None, None, False)]


def test_hmm_length_reads_leng_line(tmp_path):
    p = tmp_path / "x.hmm"
    p.write_text("HMMER3/f\nNAME x\nLENG  241\nALPH amino\n")
    assert fd.hmm_length(p) == 241


def test_hmm_length_without_leng_is_an_error(tmp_path):
    p = tmp_path / "x.hmm"
    p.write_text("HMMER3/f\nNAME x\n")
    with pytest.raises(ValueError):
        fd.hmm_length(p)


@needs_hmmer
def test_main_keeps_homologue_and_drops_decoy(tmp_path):
    hmm = tmp_path / "t.hmm"
    hits.build_hmm(DATA / "tiny_refs.faa", hmm)
    dec = tmp_path / "decoys.faa"
    dec.write_text((DATA / "tiny_proteome.faa").read_text())
    out_f, out_t = tmp_path / "kept.faa", tmp_path / "tab.tsv"
    assert fd.main(["--decoys", str(dec), "--hmm", str(hmm), "--min-score", "20", "--min-coverage", "0.3", "--out-fasta", str(out_f), "--out-table", str(out_t)]) == 0
    kept = [n for n, _, _ in fd.read_fasta(out_f)]
    assert "homologue" in kept and "decoy" not in kept
    assert out_t.read_text().splitlines()[0] == "decoy\tbest_domain_score\thmm_coverage\tkept_in_tree"
