"""Tests for the carotenoid pathway definition and its reference provenance (Task 3)."""
import csv
import subprocess
import sys
from pathlib import Path

from pathway_census.defs import load_pathway, validate

ROOT = Path(__file__).resolve().parent.parent
YAML = ROOT / "pathways" / "carotenoid.yaml"
PROV = ROOT / "refs" / "carotenoid" / "provenance.tsv"
PROV_COLS = ["step_id", "accession", "species", "pmid", "source_database",
             "retrieved_on", "verified"]
REQUIRED_STEPS = {"ggpp_synthase", "phytoene_synthase", "phytoene_desaturase",
                  "lycopene_cyclase", "carx", "cart", "astaxanthin_synthase"}


def _steps():
    return {s.id: s for s in load_pathway(YAML).steps}


def _prov_rows():
    with open(PROV, encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        assert rdr.fieldnames == PROV_COLS
        return list(rdr)


def test_carotenoid_yaml_validates():
    assert validate(load_pathway(YAML)) == []


def test_pathway_not_verified():
    assert load_pathway(YAML).verified is False


def test_step_ids_cover_spec_list():
    assert REQUIRED_STEPS <= set(_steps())


def test_fused_psy_cyclase_declared():
    s = _steps()
    assert "lycopene_cyclase" in s["phytoene_synthase"].fused_with
    assert "phytoene_synthase" in s["lycopene_cyclase"].fused_with


def test_every_required_step_has_a_ref():
    s = _steps()
    for sid in REQUIRED_STEPS:
        assert s[sid].refs, sid


def test_provenance_matches_yaml_refs():
    prov = {(r["step_id"], r["accession"]) for r in _prov_rows()}
    for s in load_pathway(YAML).steps:
        for r in s.refs:
            if r.accession:
                assert (s.id, r.accession) in prov, (s.id, r.accession)


def test_provenance_verified_agrees_with_yaml():
    prov = {(r["step_id"], r["accession"]): r for r in _prov_rows()}
    for s in load_pathway(YAML).steps:
        for r in s.refs:
            if r.accession:
                row = prov[(s.id, r.accession)]
                assert row["verified"] == str(r.verified).lower(), (s.id, r.accession)
                assert row["pmid"] == (r.pmid or ""), (s.id, r.accession)


def test_ref_fasta_ids_are_yaml_accessions():
    """Each step FASTA holds only accessions declared for that step in the YAML."""
    s = _steps()
    fdir = PROV.parent
    for sid, step in s.items():
        accs = {r.accession for r in step.refs if r.accession}
        f = fdir / f"{sid}.faa"
        if not accs:
            continue
        assert f.exists(), f
        ids = {ln[1:].split()[0] for ln in f.read_text(encoding="utf-8").splitlines()
               if ln.startswith(">")}
        assert ids <= accs, (sid, ids - accs)


def test_list_unverified_refs_runs():
    out = subprocess.run([sys.executable, str(ROOT / "bin" / "list_unverified_refs.py"),
                          str(YAML)], capture_output=True, text=True, check=True).stdout
    assert "carotenoid" in out
