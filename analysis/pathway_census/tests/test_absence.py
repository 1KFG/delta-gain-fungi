import gzip
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from pathway_census.absence import (blast_available, check_absence, same_genus_relatives, select_absences,
                                    unevaluated_genomes)

DATA = Path(__file__).parent / "data"
ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "bin" / "absence_check.py"
PW_YAML = ROOT / "pathways" / "carotenoid.yaml"
REFS = DATA / "mini_step_refs.faa"
MIN_BITS, MIN_COV = 50.0, 0.7

needs_blast = pytest.mark.skipif(not blast_available(), reason="tblastn not available (ncbi-blast module)")


def _row(lt, step, raw, has=True):
    return dict(locustag=lt, step_id=step, copy_raw=raw, copy_filtered=raw, copy_unknown=0,
                copy_flagged=0, has_proteins=has, genome_flags="")


def _genomes(rows):
    return pd.DataFrame(rows, columns=["locustag", "genus", "genus_source"])


# ---- select_absences -------------------------------------------------------
def test_select_absences_requires_relative_with_copy_and_proteins():
    m = pd.DataFrame([
        _row("A", "s1", 0), _row("A", "s2", 0),
        _row("B", "s1", 2), _row("B", "s2", 0),
        _row("C", "s1", 0, has=False), _row("C", "s2", 0, has=False),
        _row("D", "s1", 0), _row("D", "s2", 0),
    ])
    m.loc[(m.locustag == "D") & (m.step_id == "s1"), ["copy_raw", "copy_flagged"]] = [0, 0]
    rel = {"A": ["B"], "C": ["B"], "D": []}
    out = select_absences(m, None, rel)
    assert out.columns.tolist() == ["locustag", "step_id"]
    # A lacks s1 and B has it. A/s2: B lacks it too. C has no proteins. D has no relatives.
    assert out.values.tolist() == [["A", "s1"]]


def test_select_absences_uses_copy_raw_not_filtered():
    # A relative whose only copy is flagged (copy_filtered 0, copy_raw 1) still counts as a carrier.
    m = pd.DataFrame([_row("A", "s1", 0), _row("B", "s1", 1)])
    m.loc[m.locustag == "B", ["copy_filtered", "copy_flagged"]] = [0, 1]
    assert select_absences(m, None, {"A": ["B"]}).values.tolist() == [["A", "s1"]]
    # A genome whose own copy is flagged-only has copy_raw 1: not an absence.
    m2 = pd.DataFrame([_row("A", "s1", 1), _row("B", "s1", 1)])
    m2.loc[m2.locustag == "A", ["copy_filtered", "copy_flagged"]] = [0, 1]
    assert select_absences(m2, None, {"A": ["B"]}).empty


def test_select_absences_empty_has_columns():
    m = pd.DataFrame([_row("A", "s1", 1)])
    out = select_absences(m, None, {"A": []})
    assert out.empty and out.columns.tolist() == ["locustag", "step_id"]


def test_unevaluated_genomes():
    g = pd.DataFrame({"locustag": ["A", "B", "C"], "has_proteins": [True, False, False]})
    m = pd.DataFrame([_row("A", "s1", 0), _row("B", "s1", 0, has=False), _row("C", "s1", 0, has=False)])
    assert unevaluated_genomes(m, g) == ["B", "C"]


# ---- same_genus_relatives --------------------------------------------------
def test_same_genus_relatives_groups_and_excludes_self():
    g = _genomes([("A", "Mucor", "genus_field"), ("B", "Mucor", "parsed_species"),
                  ("C", "Mucor", "genus_field"), ("D", "Solo", "genus_field")])
    rel = same_genus_relatives(g)
    assert rel == {"A": ["B", "C"], "B": ["A", "C"], "C": ["A", "B"], "D": []}


def test_same_genus_relatives_never_groups_unresolved():
    g = _genomes([("A", "Ascomycota sp.", "unresolved"), ("B", "Ascomycota sp.", "unresolved"),
                  ("C", "Mucor", "genus_field"), ("D", "Mucor", "genus_field")])
    rel = same_genus_relatives(g)
    assert rel["A"] == [] and rel["B"] == []
    assert rel["C"] == ["D"] and rel["D"] == ["C"]


def test_same_genus_relatives_single_genome_genus_empty():
    g = _genomes([("A", "Lone", "genus_field")])
    assert same_genus_relatives(g) == {"A": []}


# ---- check_absence (real tblastn) -----------------------------------------
def test_genome_not_downloaded_is_unconfirmed():
    assert check_absence(None, REFS, MIN_BITS) == "unconfirmed"


@needs_blast
def test_planted_gene_found(tmp_path):
    assert check_absence(DATA / "mini_genome.fna", REFS, MIN_BITS, min_query_cover=MIN_COV,
                         scratch_dir=tmp_path) == "present_unannotated"


@needs_blast
def test_no_homology_returns_absent(tmp_path):
    assert check_absence(DATA / "mini_genome_nohit.fna", REFS, MIN_BITS, min_query_cover=MIN_COV,
                         scratch_dir=tmp_path) == "absent"


@needs_blast
def test_fragmented_when_hit_spans_less_than_fraction_of_query(tmp_path):
    # 45 of 150 aa (cover 0.30) planted. Above min_bitscore, below min_query_cover.
    assert check_absence(DATA / "mini_genome_fragment.fna", REFS, MIN_BITS, min_query_cover=MIN_COV,
                         scratch_dir=tmp_path) == "fragmented"
    # Pin the rule: with a lower cover threshold the same hit is a full call.
    assert check_absence(DATA / "mini_genome_fragment.fna", REFS, MIN_BITS, min_query_cover=0.25,
                         scratch_dir=tmp_path) == "present_unannotated"


@needs_blast
def test_gz_genome_is_streamed(tmp_path):
    gz = tmp_path / "g.scaffolds.fa.gz"
    with gzip.open(gz, "wb") as fh:
        fh.write((DATA / "mini_genome.fna").read_bytes())
    scratch = tmp_path / "scr"
    scratch.mkdir()
    assert check_absence(gz, REFS, MIN_BITS, min_query_cover=MIN_COV,
                         scratch_dir=scratch) == "present_unannotated"
    assert list(scratch.iterdir()) == []  # temp files cleaned up


def test_min_query_cover_validated():
    with pytest.raises(ValueError):
        check_absence(None, REFS, MIN_BITS, min_query_cover=1.5)


# ---- CLI -------------------------------------------------------------------
def test_cli_refuses_unverified_without_provisional(tmp_path):
    r = subprocess.run(
        ["/usr/bin/python3.12", str(CLI), "--pathway", str(PW_YAML),
         "--matrix", str(tmp_path / "m.tsv"), "--genomes", str(tmp_path / "g.tsv"),
         "--seqdb-genomes", str(tmp_path / "s.tsv"), "--refs-dir", str(tmp_path),
         "--out-dir", str(tmp_path / "o"), "--min-bitscore", "50"],
        capture_output=True, text=True)
    assert r.returncode != 0 and "verified: false" in r.stderr
    assert not (tmp_path / "o").exists()
