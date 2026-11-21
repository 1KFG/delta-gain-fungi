from pathlib import Path

import pandas as pd
import pytest

from pathway_census.assign import assign_steps
from pathway_census.defs import Pathway, Step

D = Path(__file__).parent / "data"


def _step(i, fused=()):
    return Step(id=i, name=i, ec=None, refs=[], domains=[], copy_expect="single",
                role="catalytic", fused_with=list(fused))


def _pw():
    return Pathway(id="t", verified=False, steps=[
        _step("ggpp_synthase"),
        _step("phytoene_synthase", ["lycopene_cyclase"]),
        _step("lycopene_cyclase", ["phytoene_synthase"]),
        _step("phytoene_desaturase"),
        _step("carx"),
    ])


def _load():
    cand = pd.read_csv(D / "assign_hits.tsv", sep="\t")
    plc = pd.read_csv(D / "assign_placements.tsv", sep="\t")
    return cand, plc


def test_fused_protein_gets_two_step_rows():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "fusP"]
    assert sorted(r.step_id) == ["lycopene_cyclase", "phytoene_synthase"]
    assert (r.label == "assigned").all()


def test_overlapping_domains_are_not_fusion():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "overlapQ"]
    assert (r.label == "unassigned").all()


def test_outside_every_clade_is_unassigned_not_dropped():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "outC"]
    assert len(r) == 1 and r.label.iloc[0] == "unassigned"


def test_low_support_is_unassigned():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "lowB"]
    assert len(r) == 1 and r.label.iloc[0] == "unassigned"
    assert out[out.protein_id == "soloA"].label.iloc[0] == "assigned"


def test_assigned_only_to_declared_steps():
    cand, plc = _load()
    bad = pd.concat([plc, pd.DataFrame([dict(protein_id="x", step_id="nope", clade="c", support=0.9)])])
    with pytest.raises(ValueError):
        assign_steps(cand, bad, _pw(), 0.5)


def test_one_row_per_protein_step():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    assert not out.duplicated(["protein_id", "step_id"]).any()
