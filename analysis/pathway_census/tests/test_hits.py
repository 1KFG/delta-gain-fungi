import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from pathway_census import hits
from pathway_census.defs import load_pathway

DATA = Path(__file__).parent / "data"
COLS = ["protein_id", "hmm_name", "score", "evalue", "ali_from", "ali_to", "hmm_from", "hmm_to"]

# Tests that run hmmer/mafft/mmseqs fail loudly in an env without the modules: skip with a reason.
def _tools_ok(tools):
    return all(hits.tool_available(t) for t in tools)


needs_hmmer = pytest.mark.skipif(not _tools_ok(["hmmbuild", "hmmsearch", "mafft"]),
                                 reason="hmmer/mafft not loadable (module load hmmer/3.4 mafft/7.505)")
needs_mmseqs = pytest.mark.skipif(not _tools_ok(["mmseqs"]),
                                  reason="mmseqs not loadable (module load mmseqs2/17-b804f)")


def test_parse_domtbl_sample():
    df = hits.parse_domtbl(DATA / "domtbl_sample.txt")
    assert list(df.columns) == COLS
    assert len(df) == 3
    assert df.iloc[0].protein_id == "protA"
    assert df.iloc[0].score == 168.2
    assert df.iloc[0].evalue == 2.0e-50
    assert (df.iloc[0].hmm_from, df.iloc[0].hmm_to, df.iloc[0].ali_from, df.iloc[0].ali_to) == (3, 150, 10, 157)


def test_parse_domtbl_empty(tmp_path):
    p = tmp_path / "e.domtbl"
    p.write_text("# only comments\n#\n")
    df = hits.parse_domtbl(p)
    assert df.empty and list(df.columns) == COLS


def test_candidates_keeps_best_domain_per_protein():
    df = hits.parse_domtbl(DATA / "domtbl_sample.txt")
    c = hits.candidates_for_step(df, min_score=0.0)
    assert len(c) == 2
    assert c.set_index("protein_id").loc["protA", "score"] == 168.2
    c2 = hits.candidates_for_step(df, min_score=100.0)
    assert list(c2.protein_id) == ["protA"]
    assert list(c2.columns) == COLS
    assert hits.candidates_for_step(df, min_score=1e9).empty


@needs_hmmer
def test_build_and_search_tiny(tmp_path):
    hmm = tmp_path / "t.hmm"
    hits.build_hmm(DATA / "tiny_refs.faa", hmm)
    assert hmm.stat().st_size > 0
    out = tmp_path / "t.domtbl"
    hits.search_step(hmm, DATA / "tiny_proteome.faa", out, cpus=1)
    df = hits.candidates_for_step(hits.parse_domtbl(out), min_score=0.0)
    found = set(df.protein_id)
    assert "homologue" in found
    assert "decoy" not in found


@needs_hmmer
def test_zero_hits_returns_empty_frame_with_columns(tmp_path):
    hmm = tmp_path / "t.hmm"
    hits.build_hmm(DATA / "tiny_refs.faa", hmm)
    out = tmp_path / "z.domtbl"
    hits.search_step(hmm, DATA / "unrelated_proteome.faa", out, cpus=1)
    df = hits.parse_domtbl(out)
    assert df.empty
    assert list(df.columns) == COLS


@needs_hmmer
def test_search_gzipped_proteins(tmp_path):
    import gzip
    gz = tmp_path / "p.faa.gz"
    gz.write_bytes(gzip.compress((DATA / "tiny_proteome.faa").read_bytes()))
    hmm = tmp_path / "t.hmm"
    hits.build_hmm(DATA / "tiny_refs.faa", hmm)
    out = tmp_path / "g.domtbl"
    hits.search_step(hmm, gz, out, cpus=1)
    assert "homologue" in set(hits.parse_domtbl(out).protein_id)


@needs_mmseqs
def test_dereplicate_tiny(tmp_path):
    out = tmp_path / "rep.faa"
    hits.dereplicate(DATA / "derep_in.faa", out, identity=0.9)
    ids = [l[1:].strip() for l in out.read_text().splitlines() if l.startswith(">")]
    assert len(ids) == 2
    assert "d" in ids
    assert len({"a", "b", "c"} & set(ids)) == 1


def test_dereplicate_rejects_bad_identity(tmp_path):
    with pytest.raises(ValueError):
        hits.dereplicate(DATA / "derep_in.faa", tmp_path / "x.faa", identity=1.5)


# ---- provisional guard ----
def test_unverified_pathway_refused_without_provisional():
    p = load_pathway(Path(__file__).parents[1] / "pathways" / "carotenoid.yaml")
    assert p.verified is False
    with pytest.raises(hits.UnverifiedPathwayError):
        hits.provisional_tag(p, provisional=False)
    assert hits.provisional_tag(p, provisional=True) == "PROVISIONAL"


def test_verified_pathway_gets_no_tag():
    p = load_pathway(DATA / "pathway_ok.yaml")
    if not p.verified:
        p.verified = True
    assert hits.provisional_tag(p, provisional=False) == ""


def test_tagged_name():
    assert hits.tagged_name("psy", "hmm", "PROVISIONAL") == "psy.PROVISIONAL.hmm"
    assert hits.tagged_name("psy", "hmm", "") == "psy.hmm"


BIN = Path(__file__).parents[1] / "bin"
YAML = Path(__file__).parents[1] / "pathways" / "carotenoid.yaml"


def _run_cli(args):
    return subprocess.run(["/usr/bin/python3.12", *map(str, args)], capture_output=True, text=True,
                          cwd=Path(__file__).parents[1])


def test_build_cli_refuses_unverified(tmp_path):
    r = _run_cli([BIN / "build_step_hmms.py", "--pathway", YAML, "--refs-dir", "refs/carotenoid", "--out-dir", tmp_path])
    assert r.returncode != 0
    assert "verified" in r.stderr.lower() and "--provisional" in r.stderr
    assert not any(tmp_path.iterdir())


def test_search_cli_refuses_unverified(tmp_path):
    r = _run_cli([BIN / "search_step.py", "--pathway", YAML, "--step", "ggpp_synthase", "--hmm-dir", tmp_path,
                  "--proteins", DATA / "tiny_proteome.faa", "--out-dir", tmp_path, "--min-score", "10"])
    assert r.returncode != 0
    assert "--provisional" in r.stderr


@needs_hmmer
def test_cli_provisional_stamps_filenames(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    shutil.copy(DATA / "tiny_refs.faa", refs / "ggpp_synthase.faa")
    hmms = tmp_path / "hmms"
    r = _run_cli([BIN / "build_step_hmms.py", "--pathway", YAML, "--refs-dir", refs, "--out-dir", hmms,
                  "--provisional", "--steps", "ggpp_synthase"])
    assert r.returncode == 0, r.stderr
    assert (hmms / "ggpp_synthase.PROVISIONAL.hmm").exists()
    assert (hmms / "PROVISIONAL.txt").exists()
    res = tmp_path / "res"
    r = _run_cli([BIN / "search_step.py", "--pathway", YAML, "--step", "ggpp_synthase", "--hmm-dir", hmms,
                  "--proteins", DATA / "tiny_proteome.faa", "--out-dir", res, "--min-score", "10", "--provisional"])
    assert r.returncode == 0, r.stderr
    assert (res / "ggpp_synthase.PROVISIONAL.domtbl").exists()
    cand = res / "ggpp_synthase.PROVISIONAL.candidates.tsv"
    assert cand.exists()
    assert cand.read_text().startswith("# PROVISIONAL")
    df = pd.read_csv(cand, sep="\t", comment="#")
    assert "homologue" in set(df.protein_id)


@needs_hmmer
def test_search_cli_candidates_carry_step_id(tmp_path):
    # final-review I4 (part): assign_steps needs step_id; the CLI writes it.
    refs = tmp_path / "refs"
    refs.mkdir()
    shutil.copy(DATA / "tiny_refs.faa", refs / "ggpp_synthase.faa")
    hmms = tmp_path / "hmms"
    r = _run_cli([BIN / "build_step_hmms.py", "--pathway", YAML, "--refs-dir", refs, "--out-dir", hmms,
                  "--provisional", "--steps", "ggpp_synthase"])
    assert r.returncode == 0, r.stderr
    res = tmp_path / "res"
    r = _run_cli([BIN / "search_step.py", "--pathway", YAML, "--step", "ggpp_synthase", "--hmm-dir", hmms,
                  "--proteins", DATA / "tiny_proteome.faa", "--out-dir", res, "--min-score", "10", "--provisional"])
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(res / "ggpp_synthase.PROVISIONAL.candidates.tsv", sep="\t", comment="#")
    assert list(df.columns) == hits.COLUMNS + ["step_id"]
    assert len(df) > 0 and set(df.step_id) == {"ggpp_synthase"}
    from pathway_census.assign import assign_steps
    from pathway_census.defs import load_pathway
    plc = pd.DataFrame(columns=["protein_id", "step_id", "clade", "support"])
    out = assign_steps(df, plc, load_pathway(YAML), 0.5)  # must not raise for a missing step_id
    assert set(out.step_id) == {"ggpp_synthase"}
