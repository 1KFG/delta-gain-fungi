import gzip
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from pathway_census.absence import (_covered, blast_available, check_absence, chunk_id, chunk_stem,
                                    same_genus_relatives, select_absences, unevaluated_genomes)

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
    out = select_absences(m, rel)
    assert out.columns.tolist() == ["locustag", "step_id"]
    # A lacks s1 and B has it. A/s2: B lacks it too. C has no proteins. D has no relatives.
    assert out.values.tolist() == [["A", "s1"]]


def test_select_absences_uses_copy_raw_not_filtered():
    # A relative whose only copy is flagged (copy_filtered 0, copy_raw 1) still counts as a carrier.
    m = pd.DataFrame([_row("A", "s1", 0), _row("B", "s1", 1)])
    m.loc[m.locustag == "B", ["copy_filtered", "copy_flagged"]] = [0, 1]
    assert select_absences(m, {"A": ["B"]}).values.tolist() == [["A", "s1"]]
    # A genome whose own copy is flagged-only has copy_raw 1: not an absence.
    m2 = pd.DataFrame([_row("A", "s1", 1), _row("B", "s1", 1)])
    m2.loc[m2.locustag == "A", ["copy_filtered", "copy_flagged"]] = [0, 1]
    assert select_absences(m2, {"A": ["B"]}).empty


def test_select_absences_empty_has_columns():
    m = pd.DataFrame([_row("A", "s1", 1)])
    out = select_absences(m, {"A": []})
    assert out.empty and out.columns.tolist() == ["locustag", "step_id"]


def test_unevaluated_genomes():
    m = pd.DataFrame([_row("A", "s1", 0), _row("B", "s1", 0, has=False), _row("C", "s1", 0, has=False)])
    assert unevaluated_genomes(m) == ["B", "C"]


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


# ---- chunk naming ----------------------------------------------------------
def test_chunk_id_and_stem():
    assert chunk_id("7", None) == "7"
    assert chunk_id("7", "/x/chunk_003.txt") == "7"            # explicit id wins
    assert chunk_id("", "/x/chunk_003.txt") == "chunk_003"      # else locustags basename stem
    assert chunk_id(None, None) == ""
    assert chunk_id("a/b c", None) == "a_b_c"
    assert chunk_stem("carotenoid_absence", "7") == "carotenoid_absence.chunk7"
    assert chunk_stem("carotenoid_absence", "") == "carotenoid_absence"


# ---- bitscore threshold, decoy (real tblastn) ------------------------------
# mini_genome_decoy.fna: refA with 55% of residues substituted. MEASURED 2026-11-20 (ncbi-blast 2.17.0+,
# -evalue 1e-3): one HSP, bitscore 72.8, query 3-132 of 150 (cover 130/150 = 0.867).
@needs_blast
def test_min_bitscore_threshold_decoy(tmp_path):
    g = DATA / "mini_genome_decoy.fna"
    assert check_absence(g, REFS, 72.0, min_query_cover=MIN_COV, scratch_dir=tmp_path) == "present_unannotated"
    assert check_absence(g, REFS, 72.0, min_query_cover=0.9, scratch_dir=tmp_path) == "fragmented"
    assert check_absence(g, REFS, 73.0, min_query_cover=MIN_COV, scratch_dir=tmp_path) == "absent"


# mini_genome_fragment.fna: refB[:45]. MEASURED: one HSP, bitscore 97.1, query 1-45.
@needs_blast
def test_min_bitscore_threshold_fragment(tmp_path):
    g = DATA / "mini_genome_fragment.fna"
    assert check_absence(g, REFS, 97.0, min_query_cover=MIN_COV, scratch_dir=tmp_path) == "fragmented"
    assert check_absence(g, REFS, 97.2, min_query_cover=MIN_COV, scratch_dir=tmp_path) == "absent"


@needs_blast
def test_evalue_is_a_parameter(tmp_path):
    g = DATA / "mini_genome_decoy.fna"
    assert check_absence(g, REFS, 50.0, min_query_cover=MIN_COV, scratch_dir=tmp_path,
                         evalue=1e-30) == "absent"


# mini_genome_stop.fna: refA with a TAA inserted in frame after residue 75. MEASURED: one HSP, query 1-150,
# bitscore 310 -> present_unannotated (stop codon not detected; out of scope, one synthetic case).
@needs_blast
def test_in_frame_stop_not_detected(tmp_path):
    assert check_absence(DATA / "mini_genome_stop.fna", REFS, MIN_BITS, min_query_cover=MIN_COV,
                         scratch_dir=tmp_path) == "present_unannotated"


# ---- multi-HSP union -------------------------------------------------------
def test_covered_union_no_double_count():
    assert _covered([(1, 70), (13, 118)]) == 118
    assert _covered([(1, 10), (11, 20)]) == 20        # adjacent
    assert _covered([(1, 10), (30, 40)]) == 21        # gap
    assert _covered([(5, 50), (10, 20)]) == 46        # nested


# mini_genome_overlap.fna: MEASURED two HSPs for refB, query 1-70 (151 bits) and query 13-118 (110 bits).
# Union = 118/150 = 0.787. Double counting would give 176/150 = 1.17.
@needs_blast
def test_overlapping_hsps_union_exact(tmp_path):
    g = DATA / "mini_genome_overlap.fna"
    assert check_absence(g, REFS, MIN_BITS, min_query_cover=0.8, scratch_dir=tmp_path) == "fragmented"
    assert check_absence(g, REFS, MIN_BITS, min_query_cover=0.78, scratch_dir=tmp_path) == "present_unannotated"


def test_downloaded_but_missing_fasta_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        check_absence(tmp_path / "missing.scaffolds.fa.gz", REFS, MIN_BITS, scratch_dir=tmp_path)


# ---- CLI end to end --------------------------------------------------------
def test_cli_downloaded_false_is_unconfirmed_no_bfd_fasta(tmp_path):
    (tmp_path / "m.tsv").write_text(
        "locustag\tstep_id\tcopy_raw\thas_proteins\n"
        "AAAAAAAA\tggpp_synthase\t0\tTrue\nBBBBBBBB\tggpp_synthase\t2\tTrue\nDDDDDDDD\tggpp_synthase\t0\tFalse\n")
    (tmp_path / "g.tsv").write_text(
        "# stamp line\nlocustag\tgenus\tgenus_source\thas_proteins\tdownloaded\n"
        "AAAAAAAA\tX#1\tgenus_field\tTrue\tFalse\nBBBBBBBB\tX#1\tgenus_field\tTrue\tTrue\n"
        "DDDDDDDD\tX#1\tgenus_field\tFalse\tFalse\n")
    (tmp_path / "s.tsv").write_text(f"BBBBBBBB\tb\t{tmp_path}/b.proteins.fa.gz\tx\n")
    r = subprocess.run(
        ["/usr/bin/python3.12", str(CLI), "--pathway", str(PW_YAML), "--matrix", str(tmp_path / "m.tsv"),
         "--genomes", str(tmp_path / "g.tsv"), "--seqdb-genomes", str(tmp_path / "s.tsv"),
         "--refs-dir", str(DATA), "--out-dir", str(tmp_path / "o"), "--min-bitscore", "50",
         "--chunk-id", "3", "--provisional"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    out = tmp_path / "o" / "carotenoid_absence.chunk3.PROVISIONAL.tsv"
    lines = out.read_text().splitlines()
    assert lines[0].startswith("# PROVISIONAL")
    assert lines[-2:] == ["locustag\tstep_id\toutcome\treason",
                          "AAAAAAAA\tggpp_synthase\tunconfirmed\tno_bfd_fasta"]
    un = (tmp_path / "o" / "carotenoid_unevaluated_genomes.chunk3.PROVISIONAL.tsv").read_text().splitlines()
    assert un[-2:] == ["locustag", "DDDDDDDD"]


# ---- read_tsv compression, chunk id with dots ------------------------------
def _load_cli():
    import importlib.util
    spec = importlib.util.spec_from_file_location("absence_check_cli", CLI)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TSV = b"# stamp line\nlocustag\tgenus\nAAAAAAAA\tX#1\n"


def test_read_tsv_plain_and_gz(tmp_path):
    cli = _load_cli()
    (tmp_path / "a.tsv").write_bytes(TSV)
    with gzip.open(tmp_path / "a.tsv.gz", "wb") as fh:
        fh.write(TSV)
    for name in ("a.tsv", "a.tsv.gz"):
        df = cli.read_tsv(tmp_path / name)
        assert df.columns.tolist() == ["locustag", "genus"]
        assert df.values.tolist() == [["AAAAAAAA", "X#1"]]


@pytest.mark.skipif(__import__("shutil").which("zstd") is None, reason="zstd binary not on PATH")
def test_read_tsv_zst(tmp_path):
    cli = _load_cli()
    (tmp_path / "a.tsv").write_bytes(TSV)
    subprocess.run(["zstd", "-q", str(tmp_path / "a.tsv"), "-o", str(tmp_path / "a.tsv.zst")], check=True)
    df = cli.read_tsv(tmp_path / "a.tsv.zst")
    assert df.columns.tolist() == ["locustag", "genus"]
    assert df.values.tolist() == [["AAAAAAAA", "X#1"]]


def test_chunk_id_dotted_basenames_differ():
    a, b = chunk_id("", "/x/chunk.001.txt"), chunk_id("", "/x/chunk.002.txt")
    assert (a, b) == ("chunk_001", "chunk_002")
    assert chunk_id("", "/x/list.txt") == "list"


# ---- I3: has_proteins given as strings must parse strictly ----
def _probe_matrix(as_str):
    m = pd.DataFrame([
        _row("A", "s1", 0), _row("A", "s2", 0),
        _row("B", "s1", 2), _row("B", "s2", 1),
        _row("C", "s1", 0, has=False), _row("C", "s2", 0, has=False),
    ])
    if as_str:
        m["has_proteins"] = m.has_proteins.map({True: "True", False: "False"})
    return m


def test_select_absences_parses_string_has_proteins_strictly():
    rel = {"A": ["B"], "C": ["B"]}
    real = select_absences(_probe_matrix(False), rel)
    strs = select_absences(_probe_matrix(True), rel)
    assert real.values.tolist() == [["A", "s1"], ["A", "s2"]]
    pd.testing.assert_frame_equal(real, strs)
    bad = _probe_matrix(True)
    bad.loc[0, "has_proteins"] = "nan"
    with pytest.raises(ValueError, match="has_proteins"):
        select_absences(bad, rel)


def test_unevaluated_genomes_parses_string_has_proteins_strictly():
    assert unevaluated_genomes(_probe_matrix(False)) == ["C"]
    assert unevaluated_genomes(_probe_matrix(True)) == ["C"]
    bad = _probe_matrix(True)
    bad.loc[0, "has_proteins"] = "true"
    with pytest.raises(ValueError, match="has_proteins"):
        unevaluated_genomes(bad)
