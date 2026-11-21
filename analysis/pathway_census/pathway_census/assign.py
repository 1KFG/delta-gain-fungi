"""Step assignment: tree placement plus fused-gene detection (spec section 5.4 step 3).

Placement itself is run outside this module. See README ("Step assignment") for the tool decision.
"""
from __future__ import annotations

import pandas as pd

from .defs import Pathway

OUT_COLUMNS = ["protein_id", "step_id", "support", "label"]
FUSION_MAX_OVERLAP = 0.5  # PROPOSED: two domains are distinct if they overlap < 50% of the shorter one


def _overlap_frac(a_from: int, a_to: int, b_from: int, b_to: int) -> float:
    ov = min(a_to, b_to) - max(a_from, b_from) + 1
    short = min(a_to - a_from + 1, b_to - b_from + 1)
    return max(ov, 0) / short if short > 0 else 1.0


def assign_steps(candidates: pd.DataFrame, placements: pd.DataFrame, pathway: Pathway,
                 min_support: float) -> pd.DataFrame:
    """Return one row per (protein_id, step_id): columns protein_id, step_id, support, label.

    candidates: the Task 4 candidate/domain table, concatenated over steps, with the Task 4 columns
        (protein_id, hmm_name, score, evalue, ali_from, ali_to, hmm_from, hmm_to) plus `step_id`,
        the step whose HMM produced the row. This is the ONLY table read for domain hits.
        Fusion = a protein has candidate rows for step A and step B, B is in A.fused_with (or A in
        B.fused_with), and the two alignments overlap < FUSION_MAX_OVERLAP of the shorter one.
        Both steps then get label "assigned" (domain evidence, not tree placement).
    placements: columns protein_id, step_id, clade, support. An empty/NaN/"unassigned" clade means
        the protein is outside every clade. A step_id not in `pathway` raises ValueError.
    min_support: placement support below this gives "unassigned" (the row is kept).
    Candidate proteins with no placement row get "unassigned" rows for their candidate steps
    (kept, not dropped). For a pair of placement + domain rows, assigned wins.
    """
    known = {s.id for s in pathway.steps}
    for name, df in (("placements", placements), ("candidates", candidates)):
        if "step_id" in df.columns:
            bad = sorted(set(df.step_id.dropna()) - known)
            if bad:
                raise ValueError(f"{name}: step(s) {bad} not in pathway '{pathway.id}'")

    rows: dict[tuple[str, str], tuple[float, str]] = {}

    def put(pid, sid, support, label):
        key = (pid, sid)
        old = rows.get(key)
        if old is None or (label == "assigned" and old[1] != "assigned") or \
           (label == old[1] and pd.notna(support) and (pd.isna(old[0]) or support > old[0])):
            rows[key] = (support, label)

    for r in placements.itertuples(index=False):
        clade = r.clade
        no_clade = pd.isna(clade) or str(clade).strip() in ("", "unassigned")
        ok = (not no_clade) and pd.notna(r.support) and r.support >= min_support
        put(r.protein_id, r.step_id, r.support, "assigned" if ok else "unassigned")

    fused = {s.id: set(s.fused_with) for s in pathway.steps}
    if len(candidates):
        for pid, g in candidates.groupby("protein_id"):
            # best domain per step for this protein
            best = g.sort_values("score", ascending=False).drop_duplicates("step_id")
            recs = best.to_dict("records")
            fused_steps: set[str] = set()
            for i, a in enumerate(recs):
                for b in recs[i + 1:]:
                    linked = b["step_id"] in fused.get(a["step_id"], ()) or \
                        a["step_id"] in fused.get(b["step_id"], ())
                    if linked and _overlap_frac(a["ali_from"], a["ali_to"],
                                                b["ali_from"], b["ali_to"]) < FUSION_MAX_OVERLAP:
                        fused_steps.update((a["step_id"], b["step_id"]))
            for rec in recs:
                sid = rec["step_id"]
                if sid in fused_steps:
                    prev = rows.get((pid, sid))
                    put(pid, sid, prev[0] if prev else float("nan"), "assigned")
                elif (pid, sid) not in rows and not any(k[0] == pid for k in rows):
                    put(pid, sid, float("nan"), "unassigned")

    out = pd.DataFrame([(p, s, v[0], v[1]) for (p, s), v in rows.items()], columns=OUT_COLUMNS)
    return out.sort_values(["protein_id", "step_id"]).reset_index(drop=True)
