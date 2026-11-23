#!/usr/bin/env python3
"""Side-by-side tables for several gain/loss run folders (primary root, ARD, 0.8 rule).
Usage: compare_runs.py --out <dir> name=<run dir> [name=<run dir> ...]
Reads gainloss_counts.tsv, edge_events_primary.tsv and co_events.tsv of each folder. The primary root is the first root in
gainloss_counts.tsv. Writes counts_by_run.tsv and co_events_by_run.tsv.
"""
import argparse
from pathlib import Path

import pandas as pd

STEPS = ["phytoene_synthase", "phytoene_desaturase", "lycopene_cyclase", "cart"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows, cos = [], []
    for spec in a.runs:
        name, d = spec.split("=", 1)
        d = Path(d)
        c = pd.read_csv(d / "gainloss_counts.tsv", sep="\t")
        root = c.root.iloc[0]
        c = c[(c.root == root) & (c.model == "ARD")].set_index("step")
        e = pd.read_csv(d / "edge_events_primary.tsv", sep="\t")
        tip = e.n_tips == 1
        for s in STEPS:
            ev = e[s + "_event"]
            rows.append(dict(run=name, root=root, step=s, n_tips=int(c.loc[s, "n_tips"]), present=int(c.loc[s, "n_present"]), absent=int(c.loc[s, "n_absent"]),
                             gains=int(c.loc[s, "gains"]), losses=int(c.loc[s, "losses"]), uncertain=int(c.loc[s, "uncertain"]),
                             gains_tip=int(((ev == "gain") & tip).sum()), gains_internal=int(((ev == "gain") & ~tip).sum()),
                             losses_tip=int(((ev == "loss") & tip).sum()), losses_internal=int(((ev == "loss") & ~tip).sum())))
        ce = pd.read_csv(d / "co_events.tsv", sep="\t")
        ce = ce[(ce.root == root) & (ce.event_a == ce.event_b)]
        for _, r in ce.iterrows():
            both = (e[r.step_a + "_event"] == r.event_a) & (e[r.step_b + "_event"] == r.event_b)
            cos.append(dict(run=name, step_a=r.step_a, step_b=r.step_b, event=r.event_a, observed=int(r.observed), observed_internal=int((both & ~tip).sum()),
                            expected=round(r.expected, 3), p=r.p, q=r.q))
    pd.DataFrame(rows).to_csv(out / "counts_by_run.tsv", sep="\t", index=False)
    pd.DataFrame(cos).to_csv(out / "co_events_by_run.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
