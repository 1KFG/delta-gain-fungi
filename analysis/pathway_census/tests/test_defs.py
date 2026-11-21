from pathlib import Path

import pytest

from pathway_census.defs import Pathway, Ref, Step, load_pathway, validate

DATA = Path(__file__).parent / "data"


def _step(**kw):
    base = dict(id="s1", name="n", ec=None, refs=[Ref("ACC1", "sp", None, True)],
                domains=[], copy_expect="single", role="catalytic",
                fused_with=[], optional=False)
    base.update(kw)
    return Step(**base)


def test_load_ok_fixture():
    p = load_pathway(DATA / "pathway_ok.yaml")
    assert len(p.steps) == 2
    assert p.steps[0].fused_with == ["crtY"]


def test_validate_ok_returns_empty():
    assert validate(load_pathway(DATA / "pathway_ok.yaml")) == []


def test_ref_without_pmid_or_accession_is_error():
    p = load_pathway(DATA / "pathway_bad_ref.yaml")
    errs = validate(p)
    assert len(errs) == 1
    assert "stepA" in errs[0]


def test_verified_true_requires_all_refs_verified():
    p = Pathway("x", [_step(refs=[Ref("A", "sp", None, False)])], True)
    assert validate(p)


def test_unknown_fused_with_target_is_error():
    p = Pathway("x", [_step(fused_with=["nope"])], False)
    errs = validate(p)
    assert any("nope" in e for e in errs)


def test_unknown_key_rejected(tmp_path):
    f = tmp_path / "p.yaml"
    f.write_text("id: x\nverified: false\nbogus: 1\nsteps: []\n")
    with pytest.raises(ValueError, match="bogus"):
        load_pathway(f)
