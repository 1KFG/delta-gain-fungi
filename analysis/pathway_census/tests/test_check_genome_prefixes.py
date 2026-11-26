"""Tests for bin/check_genome_prefixes.py."""
import gzip
import importlib.util
import os
import subprocess
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("check_genome_prefixes", Path(__file__).parents[1] / "bin" / "check_genome_prefixes.py")
cg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cg)


def mk(tmp_path, name, text, gz=True):
    p = tmp_path / name
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return str(p)


def test_ok_mismatch_empty_missing(tmp_path):
    ok = mk(tmp_path, "a.fa.gz", ">AAAA0001_000001-T1 desc\nMK\n")
    bad = mk(tmp_path, "b.fa", ">ZZZZ9999_000001-T1\nMK\n", gz=False)
    empty = mk(tmp_path, "c.fa.gz", "")
    r = {x["locustag"]: x for x in cg.check_all(["A", "B", "C", "D", "E"], {"A": ok, "B": bad, "C": empty, "D": str(tmp_path / "none.gz")}, None)}
    # locustags here are short names; use real-looking tags below for the ok case
    assert r["B"]["status"] == "mismatch" and r["B"]["first_prefix"] == "ZZZZ9999"
    assert r["C"]["status"] == "empty"
    assert r["D"]["status"].startswith("missing:")
    assert r["E"]["status"] == "no_seqdb_path"


def test_prefix_needs_underscore_boundary(tmp_path):
    p = mk(tmp_path, "a.fa.gz", ">AAAA0001X_1-T1\nMK\n")
    assert cg.check_one("AAAA0001", p, None)["status"] == "mismatch"
    p2 = mk(tmp_path, "b.fa.gz", ">AAAA0001_1-T1\nMK\n")
    assert cg.check_one("AAAA0001", p2, None)["status"] == "ok"


def test_mtime_flag_and_summary(tmp_path):
    old = mk(tmp_path, "old.fa.gz", ">OLD00001_1-T1\nM\n")
    new = mk(tmp_path, "new.fa.gz", ">NEW00002_1-T1\nM\n")
    os.utime(old, (1000, 1000))
    os.utime(new, (5000, 5000))
    rows = cg.check_all(["OLD00001", "OLD00003"], {"OLD00001": old, "OLD00003": new}, 3000.0)
    assert [r["newer_than_table"] for r in rows] == ["False", "True"]
    s = cg.summarise(rows)
    assert s["n_newer"] == 1 and s["n_newer_mismatch"] == 1 and s["n_older_mismatch"] == 0


def test_threads_are_capped():
    # the cap is a safety rule for the login node
    assert cg.MAX_THREADS == 8


def test_cli_with_zst_table_and_aliases(tmp_path):
    f_ok = mk(tmp_path, "ok.fa.gz", ">AAAA0001_1-T1\nM\n")
    f_bad = mk(tmp_path, "bad.fa.gz", ">BBBB0002_1-T1\nM\n")
    tbl = tmp_path / "genomes.tsv"
    tbl.write_text("locustag\thas_proteins\nAAAA0001\tTrue\nOLDTAG01\tTrue\nNOPROT01\tFalse\n")
    subprocess.run(["zstd", "-q", "-f", str(tbl), "-o", str(tbl) + ".zst"], check=True)
    sdb = tmp_path / "seqdb.tsv"
    sdb.write_text(f"AAAA0001\tx\t{f_ok}\t-\nOLDTAG01\ty\t{f_bad}\t-\n")
    out, al = tmp_path / "r.tsv", tmp_path / "aliases.tsv"
    rc = cg.main(["--genomes", str(tbl) + ".zst", "--seqdb", str(sdb), "--since", "2000-01-01T00:00:00",
                  "--out", str(out), "--aliases-out", str(al), "--threads", "99"])
    assert rc == 0
    assert len(out.read_text().splitlines()) == 3  # header + 2 genomes (has_proteins False skipped)
    lines = al.read_text().splitlines()
    assert len(lines) == 2 and lines[1].startswith("OLDTAG01\tBBBB0002\t")
