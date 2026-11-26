#!/usr/bin/env python3
"""Small summary tables for one gain/loss analysis set (read from the run folder, nothing is recomputed).

Usage: summarise_gainloss_set.py --set-dir results/gainloss/<set> [--null-dir results/gainloss/<set>_null] [--root Blastocladiomycota]

Writes into --set-dir:
  summary_steps.tsv  one row per step: tips present/absent, AIC of ER and ARD, model with the lower AIC, its rates, gains,
                     losses, uncertain branches, and gains/losses on tip and on internal edges (from edge_events_primary.tsv,
                     events of the model with the lower AIC)
  summary_pairs.tsv  one row per step pair and event type (gain/gain, loss/loss, gain/loss, loss/gain): observed co-events,
                     tree-blind permutation p and q, and (when --null-dir has null_summary.tsv) the tree-aware null mean, p and q
                     for both variants (stationary, stationary_match0.1), all-edge and internal scope, plus the Pagel result.
The edge table has one `<step>_event` column per step (the events of the model with the lower AIC, 0.8 rule).
"""
import argparse
import csv
import sys
from pathlib import Path

import pandas as pd


def tip_internal_events(edge_tab, steps):
    """Events on tip and internal edges per step. Tip ids are 1..Ntip and parent ids are above Ntip."""
    ntip = int(edge_tab["parent"].min()) - 1
    out = {}
    for s in steps:
        col = f"{s}_event"
        tip = edge_tab["child"] <= ntip
        ev = edge_tab[col]
        out[s] = {
            "gains_tip": int(((ev == "gain") & tip).sum()), "gains_internal": int(((ev == "gain") & ~tip).sum()),
            "losses_tip": int(((ev == "loss") & tip).sum()), "losses_internal": int(((ev == "loss") & ~tip).sum()),
        }
    return out, ntip


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set-dir", required=True)
    ap.add_argument("--null-dir", default=None)
    ap.add_argument("--root", default="Blastocladiomycota")
    a = ap.parse_args(argv)
    sd = Path(a.set_dir)
    cnt = pd.read_csv(sd / "gainloss_counts.tsv", sep="\t")
    cnt = cnt[cnt["root"] == a.root]
    steps = list(dict.fromkeys(cnt["step"]))
    edges = pd.read_csv(sd / "edge_events_primary.tsv", sep="\t")
    te, ntip = tip_internal_events(edges, steps)
    rows = []
    for s in steps:
        c = cnt[cnt["step"] == s].set_index("model")
        best = "ARD" if c.loc["ARD", "AIC"] < c.loc["ER", "AIC"] else "ER"
        b = c.loc[best]
        rows.append({
            "step": s, "n_tips": int(b["n_tips"]), "present": int(b["n_present"]), "absent": int(b["n_absent"]),
            "prevalence": round(b["n_present"] / b["n_tips"], 4), "AIC_ER": round(c.loc["ER", "AIC"], 1), "AIC_ARD": round(c.loc["ARD", "AIC"], 1),
            "model_lower_AIC": best, "delta_AIC_ARD_minus_ER": round(c.loc["ARD", "AIC"] - c.loc["ER", "AIC"], 1),
            "rate_01": round(b["rate_01"], 4), "rate_10": round(b["rate_10"], 4),
            "gains": int(b["gains"]), "losses": int(b["losses"]), "uncertain": int(b["uncertain"]),
            "gains_tip": te[s]["gains_tip"], "gains_internal": te[s]["gains_internal"],
            "losses_tip": te[s]["losses_tip"], "losses_internal": te[s]["losses_internal"],
            "parsimony_gains_late": int(b["pars_gains_late"]), "parsimony_gains_early": int(b["pars_gains_early"]),
            "gains_ER": int(c.loc["ER", "gains"]), "losses_ER": int(c.loc["ER", "losses"]), "uncertain_ER": int(c.loc["ER", "uncertain"]),
            "gains_ARD": int(c.loc["ARD", "gains"]), "losses_ARD": int(c.loc["ARD", "losses"]), "uncertain_ARD": int(c.loc["ARD", "uncertain"]),
        })
    pd.DataFrame(rows).to_csv(sd / "summary_steps.tsv", sep="\t", index=False)
    ce = pd.read_csv(sd / "co_events.tsv", sep="\t")
    ce = ce[ce["root"] == a.root].copy()
    ce["type"] = ce["event_a"] + "_" + ce["event_b"]
    keep = ["step_a", "step_b", "type", "n_events_a", "n_events_b", "observed", "expected", "p", "q"]
    pr = ce[keep].rename(columns={"expected": "perm_expected", "p": "perm_p", "q": "perm_q"})
    if a.null_dir and (Path(a.null_dir) / "null_summary.tsv").exists():
        ns = pd.read_csv(Path(a.null_dir) / "null_summary.tsv", sep="\t")
        for variant in ns["variant"].unique():
            for scope in ("all_edges", "internal_edges"):
                x = ns[(ns["variant"] == variant) & (ns["scope"] == scope)]
                x = x[["step_a", "step_b", "type", "observed", "null_mean", "null_lo", "null_hi", "p", "q"]]
                tag = f"{variant}.{scope}"
                x.columns = ["step_a", "step_b", "type"] + [f"{c}.{tag}" for c in ["observed", "null_mean", "null_lo", "null_hi", "p", "q"]]
                pr = pr.merge(x, on=["step_a", "step_b", "type"], how="left")
    pg = sd / "pagel"
    if pg.exists():
        parts = [pd.read_csv(f, sep="\t") for f in sorted(pg.glob("pagel.*.tsv"))]
        if parts:
            pdf = pd.concat(parts)
            keepp = [c for c in ["step_a", "step_b", "n_tips", "ok", "lik_ratio", "df", "p", "AIC_indep", "AIC_dep", "secs", "msg"] if c in pdf.columns]
            pdf = pdf[keepp].rename(columns={"p": "pagel_p", "lik_ratio": "pagel_lik_ratio", "n_tips": "pagel_n_tips", "ok": "pagel_ok", "secs": "pagel_secs", "msg": "pagel_msg"})
            pdf.to_csv(sd / "summary_pagel.tsv", sep="\t", index=False)
    pr.to_csv(sd / "summary_pairs.tsv", sep="\t", index=False)
    print(f"{sd}: {len(rows)} steps, {len(pr)} pair-by-type rows, {ntip} tips in the common tree", file=sys.stderr)


if __name__ == "__main__":
    main()
