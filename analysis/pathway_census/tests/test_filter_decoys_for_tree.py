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


def test_merge_intervals_merges_overlap_and_touching():
    assert fd.merge_intervals([(10, 20), (1, 5), (15, 30), (31, 40), (50, 60)]) == [(1, 5), (10, 40), (50, 60)]
    assert fd.merge_intervals([]) == []
    assert fd.union_coverage([(1, 50), (40, 100)], 200) == 0.5


def _df(rows):
    import pandas as pd
    return pd.DataFrame(rows, columns=hits.COLUMNS)


def test_union_passes_split_domains_that_fail_best():
    # HMM of 120 positions; domains of 15 and 12 bits cover 1-60 and 61-120; a 3-bit domain is below the limit
    df = _df([("split", "h", 15.0, 1e-3, 1, 50, 1, 60), ("split", "h", 12.0, 1e-3, 60, 110, 61, 120),
              ("split", "h", 3.0, 1e-3, 120, 130, 150, 160)])
    recs = [("split", "split", "M")]
    best = {"split": (15.0, 1, 60)}
    kept_b, _ = fd.select(recs, best, 10.0, 0.6, 120)
    assert kept_b == []                                   # best domain covers 60/120 = 0.5 < 0.6
    un = fd.union_scores(df, 10.0, 120)
    assert un["split"] == (15.0, 1, 120)                  # 60 + 60 positions, the 3-bit domain is ignored
    kept_u, table = fd.select(recs, un, 10.0, 0.6, 120)
    assert [k[0] for k in kept_u] == ["split"] and table[0][2] == 1.0


def test_union_scores_without_domain_above_limit_keeps_best_interval():
    df = _df([("w", "h", 5.0, 1e-3, 1, 50, 1, 30)])
    assert fd.union_scores(df, 10.0, 100) == {"w": (5.0, 1, 30)}
    assert fd.union_scores(_df([]), 10.0, 100) == {}


@needs_hmmer
def test_default_coverage_mode_is_best(tmp_path):
    hmm = tmp_path / "t.hmm"
    hits.build_hmm(DATA / "tiny_refs.faa", hmm)
    dec = tmp_path / "decoys.faa"
    dec.write_text((DATA / "tiny_proteome.faa").read_text())
    outs = {}
    for mode in (None, "best", "union"):
        f, t = tmp_path / f"{mode}.faa", tmp_path / f"{mode}.tsv"
        args = ["--decoys", str(dec), "--hmm", str(hmm), "--min-score", "20", "--min-coverage", "0.3", "--out-fasta", str(f), "--out-table", str(t)]
        if mode:
            args += ["--coverage-mode", mode]
        assert fd.main(args) == 0
        outs[mode] = t.read_text()
    assert outs[None] == outs["best"]
