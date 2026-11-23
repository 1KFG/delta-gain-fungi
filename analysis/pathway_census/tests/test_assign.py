from pathlib import Path

import pandas as pd
import pytest

from pathway_census.assign import assign_steps
from pathway_census.defs import Pathway, Step

D = Path(__file__).parent / "data"


def _step(i, fused=()):
    return Step(id=i, name=i, ec=None, refs=[], domains_any=[], copy_expect="single",
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


def test_fused_protein_assigned_only_where_placed():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "fusP"].set_index("step_id")
    assert sorted(r.index) == ["lycopene_cyclase", "phytoene_synthase"]
    assert r.loc["phytoene_synthase", "label"] == "assigned"
    assert r.loc["lycopene_cyclase", "label"] == "unassigned"
    assert r.fusion_candidate.all()


def test_fused_protein_with_both_placements_counts_once_each():
    cand, plc = _load()
    plc = pd.concat([plc, pd.DataFrame([dict(protein_id="fusP", step_id="lycopene_cyclase",
                                             clade="cladeII", support=0.8)])])
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "fusP"]
    assert len(r) == 2 and (r.label == "assigned").all()


def test_overlapping_domains_are_not_fusion():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "overlapQ"]
    assert len(r) == 2 and (r.label == "unassigned").all() and not r.fusion_candidate.any()


def test_mixed_case_candidate_without_placement_gets_row():
    cand, plc = _load()
    cand = pd.concat([cand, pd.DataFrame([dict(protein_id="soloA", step_id="carx", hmm_name="c",
                      score=100.0, evalue=1e-20, ali_from=1, ali_to=100, hmm_from=1, hmm_to=100)])])
    out = assign_steps(cand, plc, _pw(), 0.5).set_index(["protein_id", "step_id"])
    assert out.loc[("soloA", "ggpp_synthase"), "label"] == "assigned"
    assert out.loc[("soloA", "carx"), "label"] == "unassigned"


def test_placement_without_candidate_is_unassigned():
    cand, plc = _load()
    plc = pd.concat([plc, pd.DataFrame([dict(protein_id="ghost", step_id="carx",
                                             clade="cladeI", support=0.99)])])
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "ghost"]
    assert len(r) == 1 and r.label.iloc[0] == "unassigned"


def test_candidates_without_step_id_raises():
    cand, plc = _load()
    with pytest.raises(ValueError, match="step_id"):
        assign_steps(cand.drop(columns="step_id"), plc, _pw(), 0.5)


def test_min_support_boundary_is_inclusive():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.95)
    assert out[out.protein_id == "soloA"].label.iloc[0] == "assigned"
    out = assign_steps(cand, plc, _pw(), 0.951)
    assert out[out.protein_id == "soloA"].label.iloc[0] == "unassigned"


def test_nan_support_is_unassigned():
    cand, plc = _load()
    plc.loc[plc.protein_id == "soloA", "support"] = float("nan")
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "soloA"]
    assert len(r) == 1 and r.label.iloc[0] == "unassigned"


def test_outside_every_clade_is_unassigned_not_dropped():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "outC"]
    assert len(r) == 1 and r.label.iloc[0] == "unassigned"
    assert r.support.iloc[0] == pytest.approx(0.90)


def test_low_support_is_unassigned():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    r = out[out.protein_id == "lowB"]
    assert len(r) == 1 and r.label.iloc[0] == "unassigned"
    assert r.support.iloc[0] == pytest.approx(0.30)
    a = out[out.protein_id == "soloA"]
    assert len(a) == 1 and a.label.iloc[0] == "assigned" and a.support.iloc[0] == pytest.approx(0.95)


def test_assigned_only_to_declared_steps():
    cand, plc = _load()
    bad = pd.concat([plc, pd.DataFrame([dict(protein_id="x", step_id="nope", clade="c", support=0.9)])])
    with pytest.raises(ValueError):
        assign_steps(cand, bad, _pw(), 0.5)


def test_one_row_per_protein_step():
    cand, plc = _load()
    out = assign_steps(cand, plc, _pw(), 0.5)
    assert not out.duplicated(["protein_id", "step_id"]).any()


# ---- CLI provisional gate (final-review I1) ----
import subprocess  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "bin" / "assign_steps.py"
PW_YAML = ROOT / "pathways" / "carotenoid.yaml"


def _cli(tmp_path, *extra):
    return subprocess.run(
        ["/usr/bin/python3.12", str(CLI), "--pathway", str(PW_YAML), "--candidates", str(D / "assign_hits.tsv"),
         "--placements", str(D / "assign_placements.tsv"), "--min-support", "0.5",
         "--out", str(tmp_path / "o" / "assigned.tsv"), *extra],
        capture_output=True, text=True)


def test_cli_refuses_unverified_without_provisional(tmp_path):
    (tmp_path / "o").mkdir()
    r = _cli(tmp_path)
    assert r.returncode != 0
    assert "verified: false" in r.stderr and "--provisional" in r.stderr
    assert not any((tmp_path / "o").iterdir())


def test_cli_provisional_stamps_name_and_first_line(tmp_path):
    (tmp_path / "o").mkdir()
    r = _cli(tmp_path, "--provisional")
    assert r.returncode == 0, r.stderr
    assert not (tmp_path / "o" / "assigned.tsv").exists()
    out = tmp_path / "o" / "assigned.PROVISIONAL.tsv"
    assert out.exists()
    assert out.read_text().startswith("# PROVISIONAL")
    df = pd.read_csv(out, sep="\t", comment="#")
    assert list(df.columns) == ["protein_id", "step_id", "support", "label", "fusion_candidate"]
    assert ("soloA", "ggpp_synthase", "assigned") in set(zip(df.protein_id, df.step_id, df.label))
