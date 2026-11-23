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


def test_select_keeps_at_or_above_threshold_and_drops_missing_scores():
    recs = [("a", "a", "M"), ("b", "b", "M"), ("c", "c", "M"), ("d", "d", "M")]
    kept, table = fd.select(recs, {"a": 30.0, "b": 10.0, "c": 9.9}, 10.0)
    assert [k[0] for k in kept] == ["a", "b"]
    assert table == [("a", 30.0, True), ("b", 10.0, True), ("c", 9.9, False), ("d", None, False)]


def test_select_with_no_scores_keeps_nothing():
    kept, table = fd.select([("a", "a", "M")], {}, 5.0)
    assert kept == [] and table == [("a", None, False)]


@needs_hmmer
def test_main_keeps_homologue_and_drops_decoy(tmp_path):
    hmm = tmp_path / "t.hmm"
    hits.build_hmm(DATA / "tiny_refs.faa", hmm)
    dec = tmp_path / "decoys.faa"
    dec.write_text((DATA / "tiny_proteome.faa").read_text())
    out_f, out_t = tmp_path / "kept.faa", tmp_path / "tab.tsv"
    assert fd.main(["--decoys", str(dec), "--hmm", str(hmm), "--min-score", "20", "--out-fasta", str(out_f), "--out-table", str(out_t)]) == 0
    kept = [n for n, _, _ in fd.read_fasta(out_f)]
    assert "homologue" in kept and "decoy" not in kept
    assert out_t.read_text().splitlines()[0] == "decoy\tbest_domain_score\tkept_in_tree"
