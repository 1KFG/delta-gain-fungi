"""bin/merge_absence.py: merge absence_check chunk files (final-review fix I2).

The old README recipe globbed `*_absence.chunk*.tsv` (which also matches `.PROVISIONAL.tsv`) and
stripped every `#` line, so a provisional result could land under the canonical name. The script
must refuse mixed provisional / non-provisional input and carry the tag into the output.
"""
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "bin" / "merge_absence.py"
PROV = "# PROVISIONAL: pathway 'carotenoid' verified: false; not user-reviewed\n"
PARAMS = "# min_bitscore=50.0 min_query_cover=0.7\n"
HEAD = "locustag\tstep_id\toutcome\treason\n"


def _chunk(path, rows, prov=False, params=PARAMS):
    path.write_text((PROV if prov else "") + params + HEAD + "".join("\t".join(r) + "\n" for r in rows))
    return path


def _run(out, *inputs):
    return subprocess.run(["/usr/bin/python3.12", str(CLI), "--out", str(out), *map(str, inputs)],
                          capture_output=True, text=True)


def _read(path):
    return pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)


def test_merges_non_provisional_chunks(tmp_path):
    a = _chunk(tmp_path / "carotenoid_absence.chunk1.tsv", [("AAAA0001", "carx", "absent", "")])
    b = _chunk(tmp_path / "carotenoid_absence.chunk2.tsv",
               [("AAAA0002", "carx", "unconfirmed", "no_bfd_fasta"), ("AAAA0003", "cart", "fragmented", "")])
    out = tmp_path / "res" / "carotenoid_absence_summary.tsv"
    r = _run(out, a, b)
    assert r.returncode == 0, r.stderr
    text = out.read_text()
    assert "PROVISIONAL" not in text
    assert text.startswith(PARAMS)
    assert text.count("locustag\tstep_id") == 1
    df = _read(out)
    assert df.values.tolist() == [["AAAA0001", "carx", "absent", ""],
                                  ["AAAA0002", "carx", "unconfirmed", "no_bfd_fasta"],
                                  ["AAAA0003", "cart", "fragmented", ""]]


def test_provisional_chunks_give_tagged_name_and_first_line(tmp_path):
    a = _chunk(tmp_path / "carotenoid_absence.chunk1.PROVISIONAL.tsv", [("AAAA0001", "carx", "absent", "")], prov=True)
    b = _chunk(tmp_path / "carotenoid_absence.chunk2.PROVISIONAL.tsv", [("AAAA0002", "carx", "absent", "")], prov=True)
    out = tmp_path / "carotenoid_absence_summary.tsv"
    r = _run(out, a, b)
    assert r.returncode == 0, r.stderr
    assert not out.exists()
    tagged = tmp_path / "carotenoid_absence_summary.PROVISIONAL.tsv"
    assert tagged.read_text().startswith("# PROVISIONAL")
    assert len(_read(tagged)) == 2


def test_mixed_input_is_refused(tmp_path):
    # the old glob `*_absence.chunk*.tsv` matched both of these
    a = _chunk(tmp_path / "carotenoid_absence.chunk1.tsv", [("AAAA0001", "carx", "absent", "")])
    b = _chunk(tmp_path / "carotenoid_absence.chunk2.PROVISIONAL.tsv", [("AAAA0002", "carx", "absent", "")], prov=True)
    out = tmp_path / "o" / "carotenoid_absence_summary.tsv"
    r = _run(out, a, b)
    assert r.returncode != 0
    assert "mixed" in r.stderr.lower()
    assert not (tmp_path / "o").exists() or not any((tmp_path / "o").iterdir())


def test_provisional_header_under_canonical_name_counts_as_provisional(tmp_path):
    a = _chunk(tmp_path / "carotenoid_absence.chunk1.tsv", [("AAAA0001", "carx", "absent", "")], prov=True)
    b = _chunk(tmp_path / "carotenoid_absence.chunk2.tsv", [("AAAA0002", "carx", "absent", "")])
    r = _run(tmp_path / "s.tsv", a, b)
    assert r.returncode != 0 and "mixed" in r.stderr.lower()
    r = _run(tmp_path / "s.tsv", a)
    assert r.returncode == 0, r.stderr
    assert not (tmp_path / "s.tsv").exists()
    assert (tmp_path / "s.PROVISIONAL.tsv").read_text().startswith("# PROVISIONAL")


def test_different_params_are_refused(tmp_path):
    a = _chunk(tmp_path / "c1.tsv", [("AAAA0001", "carx", "absent", "")])
    b = _chunk(tmp_path / "c2.tsv", [("AAAA0002", "carx", "absent", "")],
               params="# min_bitscore=40.0 min_query_cover=0.7\n")
    r = _run(tmp_path / "s.tsv", a, b)
    assert r.returncode != 0 and "header" in r.stderr.lower()
    assert not (tmp_path / "s.tsv").exists()


def test_duplicate_pair_across_chunks_is_refused(tmp_path):
    a = _chunk(tmp_path / "c1.tsv", [("AAAA0001", "carx", "absent", "")])
    b = _chunk(tmp_path / "c2.tsv", [("AAAA0001", "carx", "absent", "")])
    r = _run(tmp_path / "s.tsv", a, b)
    assert r.returncode != 0 and "duplicate" in r.stderr.lower()


def test_hash_inside_a_field_is_kept(tmp_path):
    a = _chunk(tmp_path / "c1.tsv", [("AAAA0001", "carx", "absent", "note#1")])
    r = _run(tmp_path / "s.tsv", a)
    assert r.returncode == 0, r.stderr
    assert "AAAA0001\tcarx\tabsent\tnote#1" in (tmp_path / "s.tsv").read_text()
