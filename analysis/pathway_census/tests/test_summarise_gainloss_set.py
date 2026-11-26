import importlib.util
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location("summarise_gainloss_set", Path(__file__).resolve().parents[1] / "bin" / "summarise_gainloss_set.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_tip_internal_events_splits_by_child_id():
    # 3 tips (ids 1..3); parents 4 and 5 are internal; child 5 is an internal edge
    edge = pd.DataFrame({"parent": [4, 5, 5, 4], "child": [1, 2, 3, 5]})
    edge["s_event"] = ["gain", "loss", "same", "loss"]
    out, ntip = mod.tip_internal_events(edge, ["s"])
    assert ntip == 3
    assert out["s"] == {"gains_tip": 1, "gains_internal": 0, "losses_tip": 1, "losses_internal": 1}


def test_main_writes_tables(tmp_path):
    rows = []
    for model, aic, g, l, u in [("ER", 100.0, 2, 3, 1), ("ARD", 90.0, 1, 4, 2)]:
        rows.append(dict(root="R", step="s", model=model, n_tips=10, n_present=4, n_absent=6, AIC=aic, rate_01=0.1, rate_10=0.2,
                         gains=g, losses=l, uncertain=u, pars_gains_late=1, pars_gains_early=2))
    pd.DataFrame(rows).to_csv(tmp_path / "gainloss_counts.tsv", sep="\t", index=False)
    pd.DataFrame({"parent": [4, 5, 5, 4], "child": [1, 2, 3, 5], "s_event": ["gain", "loss", "same", "loss"]}).to_csv(
        tmp_path / "edge_events_primary.tsv", sep="\t", index=False)
    pd.DataFrame([dict(root="R", step_a="s", step_b="t", event_a="loss", event_b="loss", n_events_a=4, n_events_b=3, observed=2,
                       expected=0.5, p=0.01, q=0.02)]).to_csv(tmp_path / "co_events.tsv", sep="\t", index=False)
    mod.main(["--set-dir", str(tmp_path), "--root", "R"])
    s = pd.read_csv(tmp_path / "summary_steps.tsv", sep="\t")
    assert s.loc[0, "model_lower_AIC"] == "ARD" and s.loc[0, "gains"] == 1 and s.loc[0, "losses"] == 4
    assert s.loc[0, "delta_AIC_ARD_minus_ER"] == -10.0
    p = pd.read_csv(tmp_path / "summary_pairs.tsv", sep="\t")
    assert p.loc[0, "type"] == "loss_loss" and p.loc[0, "observed"] == 2
