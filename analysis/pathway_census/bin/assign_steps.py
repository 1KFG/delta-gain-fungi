#!/usr/bin/env python3.12
"""Assign candidates to steps from a placements TSV (protein_id, step_id, clade, support).

--candidates: TSV of Task 4 candidates concatenated over steps, with an added step_id column.
Tree placement itself is a follow-up (README, "Step assignment"). Input TSVs are read with
pandas (plain or .gz; leading "#" lines skipped with comment="#").
Refuses a pathway with verified: false unless --provisional is given. With --provisional the output
name carries ".PROVISIONAL" (assigned.tsv -> assigned.PROVISIONAL.tsv) and the file starts with a
"# PROVISIONAL" line.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.assign import assign_steps  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--placements", type=Path, required=True)
    ap.add_argument("--min-support", type=float, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")
    cand = pd.read_csv(a.candidates, sep="\t", comment="#")
    plc = pd.read_csv(a.placements, sep="\t", comment="#")
    res = assign_steps(cand, plc, pw, a.min_support)
    out = a.out.with_name(hits.tagged_name(a.out.stem, a.out.suffix.lstrip(".") or "tsv", tag))
    with open(out, "w") as fh:
        if tag:
            fh.write(f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n")
        res.to_csv(fh, sep="\t", index=False)
    print(f"{len(res)} rows -> {out}")


if __name__ == "__main__":
    main()
