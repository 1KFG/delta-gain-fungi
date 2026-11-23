#!/usr/bin/env python3
"""Side-by-side tables for the unfiltered run and the BUSCO-filtered runs (70, 80, relative).

Reads gainloss_counts.tsv, edge_events_primary.tsv and co_events.tsv of each run folder (primary root, ARD) and the
genus representatives table (for the per-class removal counts). Writes TSV tables into --outdir.
Usage: compare_busco_versions.py --gl results/gainloss --reps genus_representatives.tsv --outdir <dir>
"""
import argparse
import csv
import importlib.util
import statistics
from pathlib import Path

import pandas as pd

RUNS = [("unfiltered", "carotenoid"), ("busco70", "carotenoid_busco70"), ("busco65", "carotenoid_busco65")]
STEPS = ["phytoene_synthase", "phytoene_desaturase", "lycopene_cyclase", "cart"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gl", required=True)
    ap.add_argument("--reps", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    gl, out = Path(a.gl), Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    cos = []
    for name, d in RUNS:
        c = pd.read_csv(gl / d / "gainloss_counts.tsv", sep="\t")
        c = c[(c.root == "Opisthosporidia") & (c.model == "ARD")].set_index("step")
        e = pd.read_csv(gl / d / "edge_events_primary.tsv", sep="\t")
        tip = e.n_tips == 1
        for s in STEPS:
            ev = e[s + "_event"]
            rows.append(dict(run=name, step=s, n_tips=int(c.loc[s, "n_tips"]), present=int(c.loc[s, "n_present"]), absent=int(c.loc[s, "n_absent"]),
                             gains=int(c.loc[s, "gains"]), losses=int(c.loc[s, "losses"]), uncertain=int(c.loc[s, "uncertain"]),
                             gains_tip=int(((ev == "gain") & tip).sum()), gains_internal=int(((ev == "gain") & ~tip).sum()),
                             losses_tip=int(((ev == "loss") & tip).sum()), losses_internal=int(((ev == "loss") & ~tip).sum())))
        ce = pd.read_csv(gl / d / "co_events.tsv", sep="\t")
        ce = ce[(ce.root == "Opisthosporidia") & (ce.event_a == ce.event_b)]
        for _, r in ce.iterrows():
            # internal-edge co-events from the edge table (same events)
            ea, eb = e[r.step_a + "_event"], e[r.step_b + "_event"]
            both = (ea == r.event_a) & (eb == r.event_b)
            cos.append(dict(run=name, step_a=r.step_a, step_b=r.step_b, event=r.event_a, observed=int(r.observed), observed_internal=int((both & ~tip).sum()),
                            expected=round(r.expected, 3), p=r.p, q=r.q))
    pd.DataFrame(rows).to_csv(out / "counts_by_version.tsv", sep="\t", index=False)
    pd.DataFrame(cos).to_csv(out / "co_events_by_version.tsv", sep="\t", index=False)
    # per-class removal by rule, on all 1,654 genus representatives
    spec = importlib.util.spec_from_file_location("fb", Path(__file__).parent / "filter_calls_busco.py")
    fb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fb)
    reps = fb.read_rows(a.reps)
    rel, _ = fb.relative_removed(reps, 15.0, 8, 50.0)
    cls_n, tab = {}, {}
    for r in reps:
        c = (r.get("CLASS") or "").strip() or "(no class)"
        b = fb.parse_busco(r["busco_complete"])
        t = tab.setdefault(c, dict(n=0, vals=[], b60=0, b65=0, b70=0, b80=0, rel=0))
        t["n"] += 1
        t["vals"].append(b if b is not None else float("nan"))
        t["b60"] += int(b is None or b < 60)
        t["b65"] += int(b is None or b < 65)
        t["b70"] += int(b is None or b < 70)
        t["b80"] += int(b is None or b < 80)
        t["rel"] += int(r["LOCUSTAG"] in rel)
    lines = []
    for c, t in sorted(tab.items(), key=lambda kv: -kv[1]["n"]):
        v = [x for x in t["vals"] if x == x]
        lines.append(dict(CLASS=c, n_representatives=t["n"], median_busco=round(statistics.median(v), 1) if v else "", min_busco=min(v) if v else "",
                          removed_lt60=t["b60"], removed_lt65=t["b65"], removed_lt70=t["b70"], removed_lt80=t["b80"], removed_relative=t["rel"]))
    pd.DataFrame(lines).to_csv(out / "removed_by_class_and_rule.tsv", sep="\t", index=False)
    tot = {k: sum(t[k] for t in tab.values()) for k in ["b60", "b65", "b70", "b80", "rel"]}
    print("total removed: <60", tot["b60"], " <65", tot["b65"], " <70", tot["b70"], " <80", tot["b80"], " relative+floor", tot["rel"])


if __name__ == "__main__":
    main()
