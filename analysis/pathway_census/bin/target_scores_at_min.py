#!/usr/bin/env python3
"""Step 2 add-on: list every decoy and other-step reference that scores at or above the rule min_score.

bin/score_refs_decoys.py counts targets at or above the lowest own-HMM reference score. The owner rule
(bin/make_min_scores.py) uses a lower value for steps with 2 or more references (lowest leave-one-out score).
This script reads min_scores[.TAG].tsv, runs each step HMM against the same target set (every reference and
decoy of the pathway), and writes one row per target that has a hit, with the flag at_or_above_min.
Scores are the best-domain bit scores, as in score_refs_decoys.py. A target with no row had no hit at E <= 10.
"""
from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def flag_rows(step, min_score, scores, meta):
    """One row per non-own target with a score. scores: target id to best score. meta: id to (kind, step, acc)."""
    rows = []
    for tid, v in scores.items():
        kind, tstep, acc = meta[tid]
        if kind == "ref" and tstep == step:
            continue
        rows.append(dict(step=step, target_kind=kind, target_step=tstep, accession=acc, score=v,
                         min_score=min_score, at_or_above_min=v >= min_score))
    return sorted(rows, key=lambda r: -r["score"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs-dir", required=True, type=Path)
    ap.add_argument("--hmm-dir", required=True, type=Path)
    ap.add_argument("--min-scores", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tag", default="")
    ap.add_argument("--cpus", type=int, default=4)
    a = ap.parse_args(argv)
    suffix = f".{a.tag}" if a.tag else ""
    allrows = list(csv.DictReader(open(a.min_scores), delimiter="\t"))
    undefined = sorted(r["step"] for r in allrows if r["min_score"] == "")
    if undefined:  # no floor defined (make_min_scores.py --allow-undefined): nothing to flag against
        print(f"skipped, min_score empty: {undefined}", file=sys.stderr)
    mins = {r["step"]: float(r["min_score"]) for r in allrows if r["min_score"] != ""}
    steps = sorted(p.stem for p in a.refs_dir.glob("*.faa"))
    meta, targets = {}, []
    for step in steps:
        for kind, d in (("ref", a.refs_dir), ("decoy", a.refs_dir / "decoys")):
            f = d / f"{step}.faa"
            if f.exists():
                for i, (hdr, seq) in enumerate(read_fasta(f)):
                    tid = f"{kind}|{step}|{i}|{hdr.split()[0]}"
                    targets.append((tid, seq))
                    meta[tid] = (kind, step, hdr.split()[0])
    out_rows = []
    with tempfile.TemporaryDirectory() as td:
        tfa = Path(td) / "targets.faa"
        tfa.write_text("".join(f">{t}\n{s}\n" for t, s in targets))
        for step, ms in mins.items():
            dt = Path(td) / f"{step}.domtbl"
            hits.search_step(a.hmm_dir / f"{step}{suffix}.hmm", tfa, dt, a.cpus)
            d = hits.parse_domtbl(dt)
            sc = d.groupby("protein_id").score.max().to_dict() if not d.empty else {}
            out_rows += flag_rows(step, ms, sc, meta)
    cols = ["step", "target_kind", "target_step", "accession", "score", "min_score", "at_or_above_min"]
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {a.out} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
