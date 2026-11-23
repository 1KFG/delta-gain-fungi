"""Tests for bin/make_chunks.py (chunk manifest for the batched BFD search)."""
import csv
import importlib.util
import subprocess
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("make_chunks", Path(__file__).parents[1] / "bin" / "make_chunks.py")
mc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mc)


def genomes(n, bad=()):
    return [dict(locustag=f"L{i:03d}", has_proteins="False" if i in bad else "True") for i in range(n)]


def paths(n, skip=()):
    return {f"L{i:03d}": f"/x/L{i:03d}.proteins.fa.gz" for i in range(n) if i not in skip}


def test_chunk_sizes_and_ids():
    rows, excl = mc.make_chunks(genomes(10), paths(10), size=4, seed=1)
    assert excl == []
    assert [r[0] for r in rows] == [0, 0, 0, 0, 1, 1, 1, 1, 2, 2]


def test_every_genome_listed_once():
    rows, _ = mc.make_chunks(genomes(25), paths(25), size=7, seed=3)
    assert sorted(r[1] for r in rows) == [f"L{i:03d}" for i in range(25)]


def test_same_seed_same_order_other_seed_other_order():
    a, _ = mc.make_chunks(genomes(30), paths(30), size=5, seed=1)
    b, _ = mc.make_chunks(genomes(30), paths(30), size=5, seed=1)
    c, _ = mc.make_chunks(genomes(30), paths(30), size=5, seed=2)
    assert a == b
    assert a != c


def test_input_order_does_not_change_result():
    g = genomes(20)
    a, _ = mc.make_chunks(g, paths(20), size=5, seed=1)
    b, _ = mc.make_chunks(list(reversed(g)), paths(20), size=5, seed=1)
    assert a == b


def test_genome_without_proteins_is_excluded_with_reason():
    rows, excl = mc.make_chunks(genomes(6, bad=(2,)), paths(6), size=3, seed=1)
    assert "L002" not in [r[1] for r in rows]
    assert excl == [("L002", "has_proteins=False")]


def test_genome_without_seqdb_path_is_excluded_with_reason():
    rows, excl = mc.make_chunks(genomes(6), paths(6, skip=(4,)), size=3, seed=1)
    assert "L004" not in [r[1] for r in rows]
    assert excl == [("L004", "has_proteins=True;no_seqdb_path")]


@pytest.mark.parametrize("value", ["true", "TRUE", " True "])
def test_has_proteins_is_parsed_strictly_but_case_insensitively(value):
    rows, excl = mc.make_chunks([dict(locustag="L000", has_proteins=value)], paths(1), size=1, seed=1)
    assert len(rows) == 1 and excl == []


@pytest.mark.parametrize("value", ["", "NA", "<NA>", "1", "yes"])
def test_other_has_proteins_values_are_excluded_not_guessed(value):
    rows, excl = mc.make_chunks([dict(locustag="L000", has_proteins=value)], paths(1), size=1, seed=1)
    assert rows == [] and len(excl) == 1


def test_size_below_one_is_an_error():
    with pytest.raises(ValueError):
        mc.make_chunks(genomes(3), paths(3), size=0, seed=1)


def test_empty_input_gives_no_chunks():
    assert mc.make_chunks([], {}, size=5, seed=1) == ([], [])


@pytest.mark.skipif(subprocess.run(["which", "zstd"], capture_output=True).returncode != 0, reason="zstd not on PATH")
def test_cli_writes_zstd_manifest_and_excluded_list(tmp_path):
    g = tmp_path / "genomes.tsv"
    with open(g, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["locustag", "has_proteins"])
        for r in genomes(5, bad=(1,)):
            w.writerow([r["locustag"], r["has_proteins"]])
    seqdb = tmp_path / "seqdb.tsv"
    seqdb.write_text("".join(f"L{i:03d}\tname{i}\t/x/L{i:03d}.proteins.fa.gz\t/x/cds\n" for i in range(5)))
    out, excl = tmp_path / "chunks.tsv.zst", tmp_path / "excluded.tsv"
    assert mc.main(["--genomes", str(g), "--seqdb", str(seqdb), "--size", "2", "--out", str(out), "--excluded", str(excl)]) == 0
    rows = mc.read_tsv(out)
    assert list(rows[0]) == ["chunk_id", "locustag", "proteins_path"]
    assert len(rows) == 4
    assert {r["chunk_id"] for r in rows} == {"0", "1"}
    assert excl.read_text().splitlines() == ["locustag\treason", "L001\thas_proteins=False"]


def test_read_tsv_plain_file(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text("a\tb\n1\t2\n")
    assert mc.read_tsv(p) == [{"a": "1", "b": "2"}]
