#!/usr/bin/env python3.12
"""Build one HMM per pathway step from refs/<pathway>/<step>.faa (MAFFT + hmmbuild).

Refuses a pathway with verified: false unless --provisional is given. With --provisional every
output file name carries ".PROVISIONAL" and a PROVISIONAL.txt marker is written to the output dir.
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
    ap.add_argument("--refs-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--steps", nargs="*", help="step ids (default: all)")
    ap.add_argument("--provisional", action="store_true")
    a = ap.parse_args()
    pw = load_pathway(a.pathway)
    try:
        tag = hits.provisional_tag(pw, a.provisional)
    except hits.UnverifiedPathwayError as e:
        sys.exit(f"ERROR: {e}")
    steps = [s.id for s in pw.steps if not a.steps or s.id in a.steps]
    missing = [s for s in (a.steps or []) if s not in {x.id for x in pw.steps}]
    if missing:
        sys.exit(f"ERROR: unknown step(s) {missing}")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    if tag:
        (a.out_dir / "PROVISIONAL.txt").write_text(
            f"PROVISIONAL: built from pathway '{pw.id}' with verified: false. Not user-reviewed (spec 5.8).\n")
    for sid in steps:
        fa = a.refs_dir / f"{sid}.faa"
        out = a.out_dir / hits.tagged_name(sid, "hmm", tag)
        hits.build_hmm(fa, out)
        print(f"{sid}\t{out}")


if __name__ == "__main__":
    main()
