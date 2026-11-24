"""Tests for bin/search_chunk.py (batched search driver)."""
import csv
import gzip
import importlib.util
import subprocess
from pathlib import Path

import pytest

from pathway_census import hits

ROOT = Path(__file__).parents[1]
DATA = Path(__file__).parent / "data"
SPEC = importlib.util.spec_from_file_location("search_chunk", ROOT / "bin" / "search_chunk.py")
sc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sc)

needs_hmmer = pytest.mark.skipif(not all(hits.tool_available(t) for t in ["hmmbuild", "hmmsearch", "mafft"]),
                                 reason="hmmer/mafft not loadable (module load hmmer/3.4 mafft/7.505)")
needs_zstd = pytest.mark.skipif(not hits.tool_available("zstd"), reason="zstd not on PATH")


def write_genome(path: Path, locustag: str, ids, gz=True, prefix=True):
    text = "".join(f">{locustag + '_' if prefix else ''}{i}\nMKV\n" for i in ids)
    if gz:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def manifest(tmp_path, entries):
    p = tmp_path / "chunks.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chunk_id", "locustag", "proteins_path"])
        w.writerows(entries)
    return p


def test_read_chunk_returns_only_the_requested_chunk(tmp_path):
    m = manifest(tmp_path, [(0, "A", "/a"), (1, "B", "/b"), (1, "C", "/c")])
    assert [r["locustag"] for r in sc.read_chunk(m, 1)] == ["B", "C"]


def test_read_chunk_missing_chunk_stops(tmp_path):
    m = manifest(tmp_path, [(0, "A", "/a")])
    with pytest.raises(SystemExit):
        sc.read_chunk(m, 5)


@needs_zstd
def test_read_chunk_reads_zstd(tmp_path):
    m = manifest(tmp_path, [(0, "A", "/a"), (1, "B", "/b")])
    z = tmp_path / "chunks.tsv.zst"
    subprocess.run(["zstd", "-q", "-o", str(z), str(m)], check=True)
    assert [r["locustag"] for r in sc.read_chunk(z, 0)] == ["A"]


def test_join_proteins_counts_and_joins(tmp_path):
    write_genome(tmp_path / "A.fa.gz", "AAA", [1, 2, 3])
    write_genome(tmp_path / "B.fa", "BBB", [1, 2], gz=False)
    rows = [dict(locustag="AAA", proteins_path=str(tmp_path / "A.fa.gz")),
            dict(locustag="BBB", proteins_path=str(tmp_path / "B.fa"))]
    out = tmp_path / "all.faa"
    man = sc.join_proteins(rows, out)
    assert [(m["locustag"], m["status"], m["n_proteins"], m["n_bad_id_prefix"]) for m in man] == \
        [("AAA", "searched", 3, 0), ("BBB", "searched", 2, 0)]
    assert out.read_text().count(">") == 5


def test_join_proteins_counts_ids_without_locustag_prefix(tmp_path):
    write_genome(tmp_path / "A.fa.gz", "AAA", [1, 2], prefix=False)
    man = sc.join_proteins([dict(locustag="AAA", proteins_path=str(tmp_path / "A.fa.gz"))], tmp_path / "o.faa")
    assert man[0]["n_bad_id_prefix"] == 2


def test_prefix_check_needs_the_underscore(tmp_path):
    # locustag AAA must not accept an id that starts with AAAB_
    (tmp_path / "A.fa").write_text(">AAAB_1\nMKV\n")
    man = sc.join_proteins([dict(locustag="AAA", proteins_path=str(tmp_path / "A.fa"))], tmp_path / "o.faa")
    assert man[0]["n_bad_id_prefix"] == 1


def test_missing_file_is_marked_missing_not_zero_proteins_searched(tmp_path):
    man = sc.join_proteins([dict(locustag="AAA", proteins_path=str(tmp_path / "nope.fa.gz"))], tmp_path / "o.faa")
    assert man[0]["status"].startswith("missing")
    assert man[0]["n_proteins"] == 0


def test_empty_file_is_marked_missing(tmp_path):
    (tmp_path / "e.fa").write_text("")
    man = sc.join_proteins([dict(locustag="AAA", proteins_path=str(tmp_path / "e.fa"))], tmp_path / "o.faa")
    assert man[0]["status"] == "missing"


def test_truncated_gzip_is_marked_missing(tmp_path):
    full = gzip.compress(("".join(f">AAA_{i}\n" + "MKV" * 200 + "\n" for i in range(500))).encode())
    (tmp_path / "t.fa.gz").write_bytes(full[: len(full) // 2])
    man = sc.join_proteins([dict(locustag="AAA", proteins_path=str(tmp_path / "t.fa.gz"))], tmp_path / "o.faa")
    assert man[0]["status"].startswith("missing")


def _build_inputs(tmp_path):
    """Two tiny genomes, two step HMMs built from tests/data/tiny_refs.faa, one chunk."""
    prot = (DATA / "tiny_proteome.faa").read_text().splitlines()
    for lt in ("GA", "GB"):
        out = []
        for line in prot:
            out.append(f">{lt}_{line[1:]}" if line.startswith(">") else line)
        (tmp_path / f"{lt}.faa").write_text("\n".join(out) + "\n")
    m = manifest(tmp_path, [(0, "GA", str(tmp_path / "GA.faa")), (0, "GB", str(tmp_path / "GB.faa")),
                            (0, "GC", str(tmp_path / "missing.faa"))])
    hmm_dir = tmp_path / "hmms"
    hmm_dir.mkdir()
    for step in ("crtYB", "crtY"):
        hits.build_hmm(DATA / "tiny_refs.faa", hmm_dir / f"{step}.PROVISIONAL.hmm")
    ms = tmp_path / "min_scores.tsv"
    ms.write_text("step\tmin_score\tbasis\tn_refs\ncrtYB\t0\tx\t1\ncrtY\t0\tx\t1\n")
    return m, hmm_dir, ms


@needs_hmmer
@needs_zstd
def test_end_to_end_two_genomes_one_missing(tmp_path, monkeypatch):
    (tmp_path / "scr").mkdir(exist_ok=True)
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scr"))
    m, hmm_dir, ms = _build_inputs(tmp_path)
    out = tmp_path / "out"
    rc = sc.main(["--pathway", str(DATA / "pathway_ok.yaml"), "--hmm-dir", str(hmm_dir), "--chunk-manifest", str(m),
                  "--chunk-id", "0", "--min-scores", str(ms), "--out-dir", str(out), "--cpus", "1", "--provisional"])
    assert rc == 0
    names = sorted(p.name for p in out.iterdir())
    assert names == ["candidates.chunk00000.PROVISIONAL.tsv.zst", "searched.chunk00000.tsv", "timing.chunk00000.tsv"]
    raw = subprocess.run(["zstd", "-dc", str(out / names[0])], capture_output=True, check=True).stdout.decode().splitlines()
    assert raw[0] == "# PROVISIONAL"
    cands = list(csv.DictReader(raw[1:], delimiter="\t"))
    assert {c["step_id"] for c in cands} == {"crtYB", "crtY"}
    found = {c["protein_id"] for c in cands}
    assert {"GA_homologue", "GB_homologue"} <= found
    assert not any(p.endswith("_decoy") for p in found)
    srch = [l for l in (out / "searched.chunk00000.tsv").read_text().splitlines()]
    assert srch[0].startswith("# chunk=0 tag=PROVISIONAL")
    rows = list(csv.DictReader(srch[1:], delimiter="\t"))
    assert [(r["locustag"], r["status"].split(":")[0]) for r in rows] == [("GA", "searched"), ("GB", "searched"), ("GC", "missing")]
    timing = (out / "timing.chunk00000.tsv").read_text()
    assert "n_genomes\t3" in timing and "total\t" in timing and "hmmsearch:crtYB" in timing
    assert not (tmp_path / "scr" / "chunk00000" / "proteins.faa").exists()


@needs_hmmer
@needs_zstd
def test_unverified_pathway_without_provisional_is_refused(tmp_path, monkeypatch):
    (tmp_path / "scr").mkdir(exist_ok=True)
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scr"))
    m, hmm_dir, ms = _build_inputs(tmp_path)
    with pytest.raises(hits.UnverifiedPathwayError):
        sc.main(["--pathway", str(DATA / "pathway_ok.yaml"), "--hmm-dir", str(hmm_dir), "--chunk-manifest", str(m),
                 "--chunk-id", "0", "--min-scores", str(ms), "--out-dir", str(tmp_path / "out")])
    assert not (tmp_path / "out").exists()


@needs_hmmer
@needs_zstd
def test_step_without_hmm_stops_before_any_output(tmp_path, monkeypatch):
    (tmp_path / "scr").mkdir(exist_ok=True)
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scr"))
    m, hmm_dir, ms = _build_inputs(tmp_path)
    ms.write_text("step\tmin_score\tbasis\tn_refs\ncrtYB\t0\tx\t1\nghost\t0\tx\t1\n")
    with pytest.raises(SystemExit, match="ghost"):
        sc.main(["--pathway", str(DATA / "pathway_ok.yaml"), "--hmm-dir", str(hmm_dir), "--chunk-manifest", str(m),
                 "--chunk-id", "0", "--min-scores", str(ms), "--out-dir", str(tmp_path / "out"), "--provisional"])
    assert not (tmp_path / "out").exists()


def test_scratch_must_be_set(tmp_path, monkeypatch):
    monkeypatch.delenv("SCRATCH", raising=False)
    m = manifest(tmp_path, [(0, "A", "/a")])
    ms = tmp_path / "ms.tsv"
    ms.write_text("step\tmin_score\tbasis\tn_refs\ncrtYB\t0\tx\t1\n")
    with pytest.raises(SystemExit, match="SCRATCH"):
        sc.main(["--pathway", str(DATA / "pathway_ok.yaml"), "--hmm-dir", str(tmp_path), "--chunk-manifest", str(m),
                 "--chunk-id", "0", "--min-scores", str(ms), "--out-dir", str(tmp_path / "o"), "--provisional"])


@needs_hmmer
@needs_zstd
def test_empty_min_score_stops_the_search_and_names_the_step(tmp_path, monkeypatch):
    (tmp_path / "scr").mkdir(exist_ok=True)
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scr"))
    m, hmm_dir, ms = _build_inputs(tmp_path)
    ms.write_text("step\tmin_score\tbasis\tn_refs\trule_min_score\ncrtYB\t\tno floor defined: x\t1\t\n")
    with pytest.raises(SystemExit, match="crtYB"):
        sc.main(["--pathway", str(DATA / "pathway_ok.yaml"), "--hmm-dir", str(hmm_dir), "--chunk-manifest", str(m),
                 "--chunk-id", "0", "--min-scores", str(ms), "--out-dir", str(tmp_path / "out"), "--provisional"])
    assert not (tmp_path / "out").exists()
