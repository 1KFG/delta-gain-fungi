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
    assert validate(p) == ["step 's1': fused_with target 'nope' is not a step id"]


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
    p = load_pathway(_write(tmp_path, V2), refs_dir=tmp_path / 'norefs')
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
    p = Pathway("x", [_step(route="r")], False, routes=[Route("r", "R", ["s1", "ghost"], None)])
    assert validate(p) == ["route 'r': step 'ghost' is not a step id"]


def test_step_route_not_declared_is_error():
    p = Pathway("x", [_step(route="rX")], False, routes=[Route("r", "R", ["s1"], None)])
    assert validate(p) == ["step 's1': route 'rX' is not declared"]


def test_step_route_must_list_step_in_route_step_ids():
    p = Pathway("x", [_step(route="r")], False, routes=[Route("r", "R", [], None)])
    assert validate(p) == ["step 's1': route 'r' does not list the step in its step_ids"]


def test_shared_step_without_home_route_is_valid():
    p = Pathway("x", [_step(route="r1"), _step(id="s2")], False,
                routes=[Route("r1", "R1", ["s1", "s2"], None), Route("r2", "R2", ["s2"], None)])
    assert validate(p) == []


def test_duplicate_route_id_is_error():
    p = Pathway("x", [_step()], False,
                routes=[Route("r", "A", ["s1"], None), Route("r", "B", ["s1"], None)])
    assert validate(p) == ["duplicate route id 'r'"]


def test_duplicate_reaction_in_step_is_error():
    p = Pathway("x", [_step(reactions=["R1", "R2", "R1"])], False)
    assert validate(p) == ["step 's1': duplicate reaction id 'R1'"]


def test_decoy_equal_to_reference_is_error():
    p = Pathway("x", [_step(decoys=[Ref("ACC1", "sp", None, False)])], False)
    assert validate(p) == ["step 's1': decoy accession 'ACC1' is also a reference of the step"]


def test_duplicate_decoy_accession_is_error():
    d = Ref("D1", "sp", None, False)
    p = Pathway("x", [_step(decoys=[d, Ref("D1", "sp2", None, False)])], False)
    assert validate(p) == ["step 's1': duplicate decoy accession 'D1'"]


def _region_yaml(region):
    return V2.replace("region: [10, 50]", f"region: {region}")


def _root(tmp_path, length=50, faa=True, header="A1 sp|A1|X"):
    (tmp_path / "pathways").mkdir(exist_ok=True)
    refs = tmp_path / "refs" / "x"
    refs.mkdir(parents=True, exist_ok=True)
    if faa:
        (refs / "s1.faa").write_text(f">{header}\n{'M' * length}\n")
    return refs


def _load(tmp_path, region="[10, 50]", refs=None):
    f = tmp_path / "pathways" / "x.yaml"
    f.parent.mkdir(exist_ok=True)
    f.write_text(_region_yaml(region), encoding="utf-8")
    return load_pathway(f, refs_dir=refs)


def test_region_inside_sequence_ok(tmp_path):
    refs = _root(tmp_path, 50)
    p = _load(tmp_path, refs=refs)
    assert p.steps[0].refs[0].region == (10, 50)
    assert validate(p) == []


def test_region_outside_sequence_is_validate_error(tmp_path):
    refs = _root(tmp_path, 50)
    p = _load(tmp_path, "[10, 51]", refs)
    assert validate(p) == ["step 's1' ref 'A1': region [10, 51] lies outside the "
                           "sequence (length 50)"]


def test_region_default_refs_dir_is_root_refs_id(tmp_path):
    _root(tmp_path, 50)
    assert validate(_load(tmp_path, "[10, 51]")) == [
        "step 's1' ref 'A1': region [10, 51] lies outside the sequence (length 50)"]


def test_region_case_a_no_refs_dir_skips(tmp_path):
    p = _load(tmp_path, "[10, 5000]", refs=tmp_path / "norefs")
    assert p.steps[0].refs[0].region == (10, 5000)
    assert validate(p) == []


def test_region_case_b_needs_accession(tmp_path):
    refs = _root(tmp_path, 50)
    p = Pathway("x", [_step(refs=[Ref(None, "sp", "123", False, region=(1, 5))])], False,
                refs_dir=refs)
    assert validate(p) == ["step 's1': region needs an accession"]


def test_region_case_c_missing_fasta(tmp_path):
    refs = _root(tmp_path, faa=False)
    p = Pathway("x", [_step(refs=[Ref("A1", "sp", None, False, region=(1, 5))])], False,
                refs_dir=refs)
    assert validate(p) == ["step 's1' ref 'A1': region cannot be checked, no FASTA file "
                           "s1.faa in the refs directory"]


def test_region_case_c_accession_not_in_fasta(tmp_path):
    refs = _root(tmp_path, 50, header="OTHER sp|OTHER|X")
    p = Pathway("x", [_step(refs=[Ref("A1", "sp", None, False, region=(1, 5))])], False,
                refs_dir=refs)
    assert validate(p) == ["step 's1' ref 'A1': region cannot be checked, accession not "
                           "found in s1.faa"]


def test_validate_on_direct_pathway_runs_region_check(tmp_path):
    refs = _root(tmp_path, 50)
    p = Pathway("x", [_step(refs=[Ref("A1", "sp", None, False, region=(1, 60))])], False,
                refs_dir=refs)
    assert validate(p) == ["step 's1' ref 'A1': region [1, 60] lies outside the sequence "
                           "(length 50)"]


@pytest.mark.parametrize("region", ["[0, 5]", "[20, 10]", "[1]", "10"])
def test_region_malformed_rejected(tmp_path, region):
    with pytest.raises(ValueError, match="region"):
        _load(tmp_path, region, refs=tmp_path / "norefs")


# --- domain normalisation (quality) lives in test_quality.py ---
