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


# --- Task 3: hardening of the loader (Task 1 review minor findings) ---

STEP_OK = """\
id: x
verified: false
steps:
  - id: s1
    name: n
    refs:
      - accession: A1
        species: sp
        verified: false
{extra}
"""


def _write(tmp_path, text):
    f = tmp_path / "p.yaml"
    f.write_text(text, encoding="utf-8")
    return f


@pytest.mark.parametrize("key", ["domains", "fused_with"])
def test_scalar_string_list_field_rejected(tmp_path, key):
    extra = f"    copy_expect: single\n    role: catalytic\n    {key}: PF00494"
    f = _write(tmp_path, STEP_OK.format(extra=extra))
    with pytest.raises(ValueError, match=key):
        load_pathway(f)


@pytest.mark.parametrize("val", ['"false"', "'no'", "1", '"true"'])
def test_pathway_verified_must_be_bool(tmp_path, val):
    text = STEP_OK.format(extra="    copy_expect: single\n    role: catalytic")
    text = text.replace("verified: false\nsteps", f"verified: {val}\nsteps", 1)
    with pytest.raises(ValueError, match="verified"):
        load_pathway(_write(tmp_path, text))


@pytest.mark.parametrize("val", ['"false"', "'no'", "0"])
def test_ref_verified_must_be_bool(tmp_path, val):
    text = STEP_OK.format(extra="    copy_expect: single\n    role: catalytic")
    text = text.replace("        verified: false", f"        verified: {val}")
    with pytest.raises(ValueError, match="verified"):
        load_pathway(_write(tmp_path, text))


@pytest.mark.parametrize("missing", ["copy_expect", "role"])
def test_copy_expect_and_role_required(tmp_path, missing):
    lines = {"copy_expect": "    copy_expect: single", "role": "    role: catalytic"}
    extra = "\n".join(v for k, v in lines.items() if k != missing)
    with pytest.raises(ValueError, match=missing):
        load_pathway(_write(tmp_path, STEP_OK.format(extra=extra)))


def test_verified_pathway_nonoptional_step_without_refs_is_error():
    p = Pathway("x", [_step(id="s1", refs=[]), _step(id="s2", refs=[], optional=True)], True)
    errs = validate(p)
    assert any("s1" in e and "no refs" in e for e in errs)
    assert not any("s2" in e for e in errs)


def test_unverified_pathway_step_without_refs_is_ok():
    assert validate(Pathway("x", [_step(refs=[])], False)) == []
