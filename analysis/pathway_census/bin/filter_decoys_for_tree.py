#!/usr/bin/env python3
"""Keep only the decoys that can be aligned with the step domain, for use as tree leaves in step 5.

A decoy with no detectable homology to the step HMM cannot be aligned. In a maximum likelihood tree it lands on an
arbitrary branch, often inside the clade of the references. Such a decoy still serves step 2 (cutoff calibration) and
the absence check; it does not belong in the tree. Rule (PROPOSED): keep a decoy when the best domain score of the step
HMM against it, searched at E <= 1000, is at least --min-score bits AND that domain covers at least --min-coverage of the
HMM length. A short partial match (for example 44 of 241 HMM positions) scores a few bits but gives a sequence that aligns
over a fraction of the domain; in a tree it attaches to the longest branch and can land inside the clade of the references
(MEASURED for Q9KK78 in the lycopene_cyclase tree, 2026-11-22). A decoy without any hit is dropped.
Nothing under refs/ is changed. Output: a FASTA of the kept decoys and a TSV with the score of every decoy.
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pathway_census import hits  # noqa: E402


def read_fasta(path):
    name, hdr, seq = None, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, hdr, "".join(seq)
            hdr = line[1:]
            name, seq = hdr.split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, hdr, "".join(seq)


def hmm_length(hmm: Path) -> int:
    for line in open(hmm):
        if line.startswith("LENG"):
            return int(line.split()[1])
    raise ValueError(f"no LENG line in {hmm}")


def best_scores(hmm: Path, fasta: Path, cpus: int = 2) -> dict[str, tuple[float, int, int]]:
    """name -> (best domain score, hmm_from, hmm_to) of that best domain."""
    with tempfile.TemporaryDirectory() as td:
        dt = Path(td) / "d.domtbl"
        # E <= 1000 and domain E <= 1000, so weak hits of a few bits are still reported
        hits.run_tool(["hmmsearch", "--cpu", str(cpus), "--noali", "-o", "/dev/null", "--domtblout", str(dt),
                       "-E", "1000", "--domE", "1000", str(hmm), str(fasta)])
        df = hits.parse_domtbl(dt)
    if df.empty:
        return {}
    best = df.sort_values("score", ascending=False).drop_duplicates("protein_id")
    return {r.protein_id: (float(r.score), int(r.hmm_from), int(r.hmm_to)) for r in best.itertuples()}


def select(records, scores, min_score, min_coverage=0.0, hmm_len=1):
    """records: list of (name, header, seq). scores: name -> (score, hmm_from, hmm_to).
    Return (kept, table rows (name, score or None, coverage or None, kept))."""
    kept, table = [], []
    for name, hdr, seq in records:
        hit = scores.get(name)
        sc = cov = None
        keep = False
        if hit is not None:
            sc = hit[0]
            cov = (hit[2] - hit[1] + 1) / hmm_len
            keep = sc >= min_score and cov >= min_coverage
        table.append((name, sc, cov, keep))
        if keep:
            kept.append((name, hdr, seq))
    return kept, table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--decoys", required=True, type=Path)
    ap.add_argument("--hmm", required=True, type=Path)
    ap.add_argument("--min-score", required=True, type=float)
    ap.add_argument("--min-coverage", type=float, default=0.5, help="fraction of the HMM the best domain must cover (default 0.5, PROPOSED)")
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-table", required=True, type=Path)
    a = ap.parse_args(argv)
    recs = list(read_fasta(a.decoys))
    tmp = a.out_table.parent
    tmp.mkdir(parents=True, exist_ok=True)
    scores = best_scores(a.hmm, a.decoys)
    kept, table = select(recs, scores, a.min_score, a.min_coverage, hmm_length(a.hmm))
    with open(a.out_fasta, "w") as fh:
        for name, hdr, seq in kept:
            fh.write(f">{hdr}\n{seq}\n")
    with open(a.out_table, "w") as fh:
        fh.write("decoy\tbest_domain_score\thmm_coverage\tkept_in_tree\n")
        for name, sc, cov, keep in table:
            fh.write(f"{name}\t{'' if sc is None else sc}\t{'' if cov is None else round(cov, 3)}\t{keep}\n")
    print(f"kept {len(kept)} of {len(recs)} decoys (min score {a.min_score}, min coverage {a.min_coverage})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
