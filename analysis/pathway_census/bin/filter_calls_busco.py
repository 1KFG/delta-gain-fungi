#!/usr/bin/env python3
"""Set every step to NA for tips whose genus representative has busco_complete below a threshold.

Reads a calls.tsv from make_calls_table.py and genus_representatives.tsv (LOCUSTAG, busco_complete).
Writes the filtered calls.tsv (adds columns `busco_complete` and `busco_filtered`) and a per-class summary
(`busco_filter_by_class.tsv`: class, n_tips, n_removed, fraction) beside it. A tip with no busco value is removed
and counted in the summary line `no_busco_value`. Tips that were already NA stay NA (not counted as removed).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_rows(path):
    with open(path, newline="") as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def parse_busco(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def step_columns(row):
    return [k[: -len("_source")] for k in row if k.endswith("_source")]


def relative_removed(reps, drop=15.0, min_class_n=8, floor=50.0):
    """Lineage-relative rule. reps: rows with LOCUSTAG, CLASS, busco_complete.
    Remove a tip if its busco_complete is more than `drop` points below the median of its class (classes with at least
    `min_class_n` genus representatives; the median uses every representative of the class) OR below the absolute `floor`
    in any class. A tip without a busco value is removed. Returns (set of removed ids, set removed by the floor only)."""
    import statistics
    by = {}
    for r in reps:
        b = parse_busco(r.get("busco_complete"))
        if b is not None and (r.get("CLASS") or "").strip():
            by.setdefault(r["CLASS"].strip(), []).append(b)
    med = {c: statistics.median(v) for c, v in by.items() if len(v) >= min_class_n}
    removed, floor_only = set(), set()
    for r in reps:
        b = parse_busco(r.get("busco_complete")); c = (r.get("CLASS") or "").strip()
        rel = b is not None and c in med and b < med[c] - drop
        ab = b is None or b < floor
        if rel or ab:
            removed.add(r["LOCUSTAG"])
            if ab and not rel:
                floor_only.add(r["LOCUSTAG"])
    return removed, floor_only


def filter_calls(calls, busco, min_busco, remove_ids=None):
    """calls: list of dict rows. busco: locustag -> float or None. Returns (rows, summary rows, n_no_busco).
    remove_ids: optional set of ids to remove instead of the absolute `min_busco` rule."""
    out, per_class, no_busco = [], {}, 0
    for r in calls:
        r = dict(r)
        steps = step_columns(r)
        b = busco.get(r["locustag"])
        r["busco_complete"] = "" if b is None else f"{b:g}"
        had_data = r.get("searched") == "1"
        remove = (r['locustag'] in remove_ids) if remove_ids is not None else (b is None or b < min_busco)
        r["busco_filtered"] = "1" if (remove and had_data) else "0"
        if b is None and had_data:
            no_busco += 1
        if remove:
            for s in steps:
                r[s] = "NA"
        cls = r.get("CLASS") or "(no class)"
        d = per_class.setdefault(cls, [0, 0])
        if had_data:
            d[0] += 1
            d[1] += int(remove)
        out.append(r)
    summ = [{"class": c, "n_tips_with_data": v[0], "n_removed": v[1], "fraction_removed": f"{(v[1] / v[0]) if v[0] else 0:.3f}"}
            for c, v in sorted(per_class.items(), key=lambda kv: -kv[1][0])]
    return out, summ, no_busco


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--reps", required=True)
    ap.add_argument("--min-busco", type=float, default=80.0)
    ap.add_argument("--relative-drop", type=float, default=None,
                    help="use the lineage-relative rule: drop tips more than this many points below their class median")
    ap.add_argument("--relative-min-class-n", type=int, default=8)
    ap.add_argument("--abs-floor", type=float, default=50.0, help="with --relative-drop: also drop tips below this absolute value")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    calls = read_rows(a.calls)
    reps = read_rows(a.reps)
    busco = {r["LOCUSTAG"]: parse_busco(r.get("busco_complete")) for r in reps}
    remove_ids = None
    if a.relative_drop is not None:
        remove_ids, floor_only = relative_removed(reps, a.relative_drop, a.relative_min_class_n, a.abs_floor)
        print(f"relative rule: removed {len(remove_ids)} of {len(reps)} representatives ({len(floor_only)} only by the absolute floor)", file=sys.stderr)
    rows, summ, nb = filter_calls(calls, busco, a.min_busco, remove_ids)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(out.with_name("busco_filter_by_class.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summ)
    steps = step_columns(rows[0])
    for s in steps:
        v = [r[s] for r in rows]
        print(f"step={s} present={v.count('1')} absent={v.count('0')} NA={v.count('NA')}", file=sys.stderr)
    print(f"min_busco={a.min_busco} removed={sum(r['busco_filtered'] == '1' for r in rows)} no_busco_value={nb}", file=sys.stderr)


if __name__ == "__main__":
    main()
