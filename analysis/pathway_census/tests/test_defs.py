from pathlib import Path

import pytest

from pathway_census.defs import Pathway, Ref, Route, Step, load_pathway, validate

DATA = Path(__file__).parent / "data"


def _step(**kw):
    base = dict(id="s1", name="n", ec=None, refs=[Ref("ACC1", "sp", None, True)],
                domains_any=[], copy_expect="single", role="catalytic",
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


@pytest.mark.parametrize("key", ["domains_any", "fused_with"])
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


# --- Schema v2 (Task S1) ---

V2 = """\
id: x
verified: false
routes:
  - id: r1
    name: route one
    step_ids: [s1, s2]
    note: first
  - id: r2
    name: route two
    step_ids: [s2]
steps:
  - id: s1
    name: n1
    copy_expect: single
    role: catalytic
    route: r1
    ec: "1.1.1.1"
    domains_any: [PF1, PF2]
    reactions: [RHEA:1, RHEA:2]
    lineage_note: only in Pezizomycotina
    refs:
      - accession: A1
        species: sp
        verified: true
        gene_name: crtX
        region: [10, 50]
    decoys:
      - accession: D1
        species: other
        verified: false
  - id: s2
    name: n2
    copy_expect: multi
    role: auxiliary
"""


def test_load_v2_fixture_with_two_routes(tmp_path):
    p = load_pathway(_write(tmp_path, V2))
    assert [r.id for r in p.routes] == ["r1", "r2"]
    assert p.routes[0].step_ids == ["s1", "s2"]
    assert p.routes[0].note == "first" and p.routes[1].note is None
    s1, s2 = p.steps
    assert s1.domains_any == ["PF1", "PF2"]
    assert s1.route == "r1" and s1.ec == "1.1.1.1"
    assert s1.reactions == ["RHEA:1", "RHEA:2"]
    assert s1.lineage_note == "only in Pezizomycotina"
    assert s1.refs[0].gene_name == "crtX" and s1.refs[0].region == (10, 50)
    assert [d.accession for d in s1.decoys] == ["D1"]
    assert s2.route is None and s2.reactions == [] and s2.decoys == []
    assert validate(p) == []


def test_old_domains_key_rejected(tmp_path):
    extra = "    copy_expect: single\n    role: catalytic\n    domains: [PF1]"
    with pytest.raises(ValueError, match="use domains_any"):
        load_pathway(_write(tmp_path, STEP_OK.format(extra=extra)))


def test_route_with_unknown_step_is_error():
    p = Pathway("x", [_step()], False, routes=[Route("r", "R", ["s1", "ghost"], None)])
    errs = validate(p)
    assert any("ghost" in e and "route 'r'" in e for e in errs)


def test_step_route_not_declared_is_error():
    p = Pathway("x", [_step(route="rX")], False, routes=[Route("r", "R", ["s1"], None)])
    assert any("rX" in e for e in validate(p))


def test_duplicate_reaction_in_step_is_error():
    p = Pathway("x", [_step(reactions=["R1", "R2", "R1"])], False)
    assert any("R1" in e and "duplicate" in e for e in validate(p))


def test_decoy_equal_to_reference_is_error():
    p = Pathway("x", [_step(decoys=[Ref("ACC1", "sp", None, False)])], False)
    assert any("ACC1" in e and "decoy" in e for e in validate(p))


def _region_yaml(region):
    return V2.replace("region: [10, 50]", f"region: {region}")


def _with_fasta(tmp_path, region, length):
    (tmp_path / "pathways").mkdir()
    (tmp_path / "refs" / "x").mkdir(parents=True)
    (tmp_path / "refs" / "x" / "s1.faa").write_text(f">A1 sp|A1|X\n{'M' * length}\n")
    f = tmp_path / "pathways" / "x.yaml"
    f.write_text(_region_yaml(region), encoding="utf-8")
    return f


def test_region_inside_sequence_ok(tmp_path):
    assert load_pathway(_with_fasta(tmp_path, "[10, 50]", 50)).steps[0].refs[0].region == (10, 50)


def test_region_outside_sequence_rejected(tmp_path):
    with pytest.raises(ValueError, match="region"):
        load_pathway(_with_fasta(tmp_path, "[10, 51]", 50))


@pytest.mark.parametrize("region", ["[0, 5]", "[20, 10]", "[1]", "10"])
def test_region_malformed_rejected(tmp_path, region):
    with pytest.raises(ValueError, match="region"):
        load_pathway(_write(tmp_path, _region_yaml(region)))


def test_region_length_check_skipped_without_fasta(tmp_path):
    p = load_pathway(_write(tmp_path, _region_yaml("[10, 5000]")))
    assert p.steps[0].refs[0].region == (10, 5000)
