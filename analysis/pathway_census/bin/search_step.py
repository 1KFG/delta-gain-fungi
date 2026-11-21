#!/usr/bin/env python3.12
"""Search one protein FASTA (plain or .gz) with one step HMM and write the domtbl and candidate TSV.

Refuses a pathway with verified: false unless --provisional is given. With --provisional the output
names carry ".PROVISIONAL" and the candidate TSV starts with a "# PROVISIONAL" line.
--min-score has no default: the pilot derives it from the reference proteins' own scores (README).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pathway", type=Path, required=True)
    ap.add_argument("--step", required=True)
    ap.add_argument("--hmm-dir", type=Path, required=True)
    ap.add_argument("--proteins", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-score", type=float, required=True)
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")
    if a.step not in {s.id for s in pw.steps}:
        sys.exit(f"ERROR: unknown step '{a.step}'")
    hmm = a.hmm_dir / hits.tagged_name(a.step, "hmm", tag)
    if not hmm.exists():
        sys.exit(f"ERROR: HMM not found: {hmm}")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    dom = a.out_dir / hits.tagged_name(a.step, "domtbl", tag)
    hits.search_step(hmm, a.proteins, dom, a.cpus)
    cand = hits.candidates_for_step(hits.parse_domtbl(dom), a.min_score)
    out = a.out_dir / hits.tagged_name(a.step, "candidates.tsv", tag)
    with open(out, "w") as fh:
        if tag:
            fh.write(f"# {tag}: pathway '{pw.id}' verified: false; not user-reviewed\n")
        fh.write(f"# step={a.step} proteins={a.proteins} min_score={a.min_score}\n")
        cand.to_csv(fh, sep="\t", index=False)
    print(f"{a.step}\t{len(cand)} candidates\t{out}")


if __name__ == "__main__":
    main()
