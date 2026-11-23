"""Tests for bin/make_min_scores.py (owner-approved min_score rule, 2026-11-22)."""
import csv
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("make_min_scores", Path(__file__).parents[1] / "bin" / "make_min_scores.py")
mms = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mms)


def row(step, n_refs, loo_min="", max_decoy="", **kw):
    return dict(step=step, n_refs=str(n_refs), loo_min=str(loo_min), max_decoy_score=str(max_decoy), **kw)


def test_multi_reference_step_uses_lowest_leave_one_out_score():
    r = mms.min_scores([row("a", 5, loo_min=174.1, max_decoy=29.3, proposed_min_score=347.5)])
    assert r == [{"step": "a", "min_score": 174.1, "basis": "lowest leave-one-out score", "n_refs": 5,
                  "rule_min_score": 174.1}]


def test_two_references_is_multi_reference():
    assert mms.min_scores([row("a", 2, loo_min=615.8, max_decoy=19.8)])[0]["min_score"] == 615.8


def test_single_reference_step_uses_highest_decoy_score_as_floor():
    r = mms.min_scores([row("a", 1, max_decoy=333.6)])[0]
    assert r["min_score"] == 333.6
    assert "single reference" in r["basis"]


def test_float_formatted_n_refs_is_accepted():
    assert mms.min_scores([row("a", "1.0", max_decoy=8.0)])[0]["min_score"] == 8.0


@pytest.mark.parametrize("bad", [row("a", 1, max_decoy=""), row("a", 1, max_decoy="nan"),
                                 row("a", 3, loo_min=""), row("a", 3, loo_min="nan", max_decoy=5)])
def test_missing_basis_value_stops_with_error_not_a_default(bad):
    with pytest.raises(ValueError, match="step a"):
        mms.min_scores([bad])


def test_multi_reference_does_not_fall_back_to_decoy_when_loo_is_missing():
    with pytest.raises(ValueError):
        mms.min_scores([row("a", 4, loo_min="", max_decoy=10)])


def test_cli_round_trip(tmp_path):
    src = tmp_path / "step_scores.tsv"
    with open(src, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["step", "n_refs", "loo_min", "max_decoy_score"], delimiter="\t")
        w.writeheader()
        w.writerow(dict(step="s1", n_refs=5, loo_min=174.1, max_decoy_score=29.3))
        w.writerow(dict(step="s2", n_refs=1, loo_min="", max_decoy_score=8.0))
    out = tmp_path / "min_scores.tsv"
    assert mms.main(["--step-scores", str(src), "--out", str(out)]) == 0
    got = list(csv.DictReader(open(out), delimiter="\t"))
    assert [(g["step"], float(g["min_score"])) for g in got] == [("s1", 174.1), ("s2", 8.0)]
    assert list(got[0]) == ["step", "min_score", "basis", "n_refs", "rule_min_score"]


def test_scale_lowers_one_step_and_keeps_the_rule_value():
    rows = [row("a", 2, loo_min=1000.0), row("b", 2, loo_min=500.0)]
    r = {x["step"]: x for x in mms.min_scores(rows, {"a": 0.75})}
    assert r["a"]["min_score"] == 750.0 and r["a"]["rule_min_score"] == 1000.0
    assert "scaled x0.75" in r["a"]["basis"]
    assert r["b"]["min_score"] == 500.0 and "scaled" not in r["b"]["basis"]


@pytest.mark.parametrize("factor", [0, -0.5, 1.5])
def test_scale_factor_outside_zero_to_one_is_an_error(factor):
    with pytest.raises(ValueError, match="factor"):
        mms.min_scores([row("a", 2, loo_min=10)], {"a": factor})


def test_scale_naming_an_unknown_step_is_an_error():
    with pytest.raises(ValueError, match="ghost"):
        mms.min_scores([row("a", 2, loo_min=10)], {"ghost": 0.75})


def test_cli_scale_option(tmp_path):
    src = tmp_path / "s.tsv"
    src.write_text("step\tn_refs\tloo_min\tmax_decoy_score\na\t2\t1000\t5\nb\t2\t200\t5\n")
    out = tmp_path / "o.tsv"
    assert mms.main(["--step-scores", str(src), "--out", str(out), "--scale", "a=0.75"]) == 0
    got = {g["step"]: g for g in csv.DictReader(open(out), delimiter="\t")}
    assert float(got["a"]["min_score"]) == 750.0 and float(got["a"]["rule_min_score"]) == 1000.0
    assert float(got["b"]["min_score"]) == 200.0
