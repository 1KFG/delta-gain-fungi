#!/usr/bin/env python3
"""Apply the owner-approved min_score rule (2026-11-22) to a step 2 result.

Rule: a step with two or more references uses its lowest leave-one-out score.
A step with one reference uses its highest decoy score as a floor. Tree placement makes the real calls.
Input: step_scores[.TAG].tsv from bin/score_refs_decoys.py. Output: TSV step, min_score, basis, n_refs.
A step without a usable value stops the script with an error. It never gets a silent default.

Optional --floor-override STEP=VALUE (repeatable) replaces the rule value of one step by an explicit value. The
column rule_min_score keeps the rule value (empty if the rule gives none) and the basis column names the override.
--override-basis TEXT sets the basis text (default: "explicit floor override"). The owner decision that
justifies an override belongs in that text.

Optional --allow-undefined: a step whose rule has no value (a single-reference step whose decoys do not score on
the HMM, or a step with no decoy) is written with an EMPTY min_score and the basis "no floor defined: ...". The
script does not invent a value. bin/target_scores_at_min.py skips such a step and says so;
bin/search_chunk.py stops with an error that names the step. Without this option the script stops on such a step.

Optional --scale STEP=FACTOR (repeatable) multiplies the rule value of one step by FACTOR (0 < FACTOR <= 1),
for an exploratory lower cutoff. The column rule_min_score always keeps the unscaled rule value, and the basis
text says that the value was scaled. Because a search keeps every hit at or above min_score, a lower cutoff
gives a superset of the rule hits, so the rule result is recovered by filtering on score.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys


def _num(x):
    if x is None or x == "":
        return None
    v = float(x)
    return None if math.isnan(v) else v


def min_scores(rows, scale=None, floor_override=None, override_basis="explicit floor override", allow_undefined=False):
    scale = scale or {}
    floor_override = floor_override or {}
    unknown = set(scale) - {r["step"] for r in rows}
    if unknown:
        raise ValueError(f"--scale names steps that are not in the input: {sorted(unknown)}")
    unknown = set(floor_override) - {r["step"] for r in rows}
    if unknown:
        raise ValueError(f"--floor-override names steps that are not in the input: {sorted(unknown)}")
    both = set(scale) & set(floor_override)
    if both:
        raise ValueError(f"--scale and --floor-override both name {sorted(both)}")
    for step, v in floor_override.items():
        if not v > 0:
            raise ValueError(f"--floor-override {step}: value must be > 0, got {v}")
    for step, f in scale.items():
        if not 0 < f <= 1:
            raise ValueError(f"--scale {step}: factor must be > 0 and <= 1, got {f}")
    out = []
    for r in rows:
        step, n = r["step"], int(float(r["n_refs"]))
        if n >= 2:
            v, basis = _num(r.get("loo_min")), "lowest leave-one-out score"
        else:
            v, basis = _num(r.get("max_decoy_score")), "highest decoy score (single reference, floor)"
        if v is None and step in floor_override:
            rec = {"step": step, "min_score": floor_override[step], "basis": override_basis, "n_refs": n,
                   "rule_min_score": ""}
            out.append(rec)
            continue
        if v is None:
            if not allow_undefined:
                raise ValueError(f"step {step}: no value for rule basis '{basis}'")
            why = ("the decoys do not score on the HMM or there is no decoy" if n < 2
                   else "no leave-one-out score")
            out.append({"step": step, "min_score": "", "basis": f"no floor defined: {why}", "n_refs": n,
                        "rule_min_score": ""})
            continue
        if step in floor_override:
            out.append({"step": step, "min_score": floor_override[step], "basis": override_basis, "n_refs": n,
                        "rule_min_score": v})
            continue
        rec = {"step": step, "min_score": v, "basis": basis, "n_refs": n, "rule_min_score": v}
        if step in scale:
            rec["min_score"] = round(v * scale[step], 1)
            rec["basis"] = f"{basis}; scaled x{scale[step]} (exploratory)"
        out.append(rec)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step-scores", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--scale", action="append", default=[], metavar="STEP=FACTOR")
    ap.add_argument("--floor-override", action="append", default=[], metavar="STEP=VALUE")
    ap.add_argument("--override-basis", default="explicit floor override")
    ap.add_argument("--allow-undefined", action="store_true")
    a = ap.parse_args(argv)
    override = {}
    for item in a.floor_override:
        step, _, v = item.partition("=")
        override[step] = float(v)
    scale = {}
    for item in a.scale:
        step, _, f = item.partition("=")
        scale[step] = float(f)
    with open(a.step_scores) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    res = min_scores(rows, scale, override, a.override_basis, a.allow_undefined)
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["step", "min_score", "basis", "n_refs", "rule_min_score"], delimiter="\t")
        w.writeheader()
        w.writerows(res)
    print(f"wrote {a.out} ({len(res)} steps)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
