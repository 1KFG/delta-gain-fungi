#!/usr/bin/env python3.12
"""Assign candidates to steps from a placements TSV (protein_id, step_id, clade, support).

--candidates: TSV of Task 4 candidates concatenated over steps, with an added step_id column.
Tree placement itself is a follow-up (README, "Step assignment"). Plain or .gz/.zst TSV input.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census.assign import assign_steps  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--placements", type=Path, required=True)
    ap.add_argument("--min-support", type=float, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    cand = pd.read_csv(a.candidates, sep="\t", comment="#")
    plc = pd.read_csv(a.placements, sep="\t", comment="#")
    out = assign_steps(cand, plc, load_pathway(a.pathway), a.min_support)
    out.to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
