"""Step assignment: tree placement plus fused-gene detection (spec section 5.4 step 3).

Placement itself is run outside this module. See README ("Step assignment") for the tool decision.
"""
from __future__ import annotations

import pandas as pd

from .defs import Pathway

OUT_COLUMNS = ["protein_id", "step_id", "support", "label", "fusion_candidate"]
FUSION_MAX_OVERLAP = 0.5  # PROPOSED: two domains are distinct if they overlap < 50% of the shorter one


def _overlap_frac(a_from: int, a_to: int, b_from: int, b_to: int) -> float:
    ov = min(a_to, b_to) - max(a_from, b_from) + 1
    short = min(a_to - a_from + 1, b_to - b_from + 1)
    return max(ov, 0) / short if short > 0 else 1.0


def assign_steps(candidates: pd.DataFrame, placements: pd.DataFrame, pathway: Pathway,
                 min_support: float) -> pd.DataFrame:
    """Return one row per (protein_id, step_id): protein_id, step_id, support, label, fusion_candidate.

    candidates: Task 4 candidate/domain tables concatenated over steps, Task 4 columns plus `step_id`
        (the step whose HMM gave the row). Build with
        pd.concat([candidates_for_step(...).assign(step_id=sid) for ...]). Missing `step_id` raises
        ValueError. This is the ONLY table read for domain hits.
    placements: columns protein_id, step_id, clade, support. Empty/NaN/"unassigned" clade means
        outside every clade. A step_id not in `pathway` (either table) raises ValueError.
    Rule: label is "assigned" only if (protein, step) is in `candidates` (HMM evidence) AND a placement
        row for it has a clade and support >= min_support (boundary inclusive; NaN support fails).
        Otherwise "unassigned". Rows are never dropped: every candidate pair gets a row, and a
        placement pair absent from `candidates` is kept as "unassigned" (no HMM evidence).
    fusion_candidate: True for pairs whose protein has HMM hits to two steps linked by `fused_with`
        with alignments overlapping < FUSION_MAX_OVERLAP of the shorter one. Fusion decides which
        steps a protein may carry; it never assigns by itself. support is NaN when no placement exists.
    """
    if "step_id" not in candidates.columns:
        raise ValueError("candidates lacks a 'step_id' column; add it with df.assign(step_id=...)")
    known = {s.id for s in pathway.steps}
    for name, df in (("placements", placements), ("candidates", candidates)):
        bad = sorted(set(df.step_id.dropna()) - known)
        if bad:
            raise ValueError(f"{name}: step(s) {bad} not in pathway '{pathway.id}'")

    # best qualifying-or-highest placement per pair
    plc: dict[tuple[str, str], tuple[float, bool]] = {}
    for r in placements.itertuples(index=False):
        no_clade = pd.isna(r.clade) or str(r.clade).strip() in ("", "unassigned")
        ok = (not no_clade) and pd.notna(r.support) and r.support >= min_support
        key = (r.protein_id, r.step_id)
        old = plc.get(key)
        if old is None or (ok and not old[1]) or (ok == old[1] and pd.notna(r.support) and
                                                  (pd.isna(old[0]) or r.support > old[0])):
            plc[key] = (r.support, ok)

    cand_pairs: set[tuple[str, str]] = set()
    fusion: set[tuple[str, str]] = set()
    fused = {s.id: set(s.fused_with) for s in pathway.steps}
    for pid, g in candidates.groupby("protein_id"):
        recs = g.sort_values("score", ascending=False).drop_duplicates("step_id").to_dict("records")
        for rec in recs:
            cand_pairs.add((pid, rec["step_id"]))
        for i, a in enumerate(recs):
            for b in recs[i + 1:]:
                linked = b["step_id"] in fused.get(a["step_id"], ()) or \
                    a["step_id"] in fused.get(b["step_id"], ())
                if linked and _overlap_frac(a["ali_from"], a["ali_to"],
                                            b["ali_from"], b["ali_to"]) < FUSION_MAX_OVERLAP:
                    fusion.update(((pid, a["step_id"]), (pid, b["step_id"])))

    out = []
    for key in cand_pairs | set(plc):
        support, ok = plc.get(key, (float("nan"), False))
        label = "assigned" if (ok and key in cand_pairs) else "unassigned"
        out.append((key[0], key[1], support, label, key in fusion))
    df = pd.DataFrame(out, columns=OUT_COLUMNS)
    return df.sort_values(["protein_id", "step_id"]).reset_index(drop=True)
