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


def merge_intervals(intervals):
    """Merge overlapping or touching 1-based inclusive (from, to) intervals; return the merged list."""
    out = []
    for lo, hi in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [(a, b) for a, b in out]


def union_coverage(intervals, hmm_len):
    """Fraction of the HMM positions covered by the union of the intervals."""
    return sum(b - a + 1 for a, b in merge_intervals(intervals)) / hmm_len


def domain_hits(hmm: Path, fasta: Path, cpus: int = 2):
    """Parsed domtbl DataFrame of hmmsearch -E 1000 --domE 1000 (empty when there is no hit)."""
    with tempfile.TemporaryDirectory() as td:
        dt = Path(td) / "d.domtbl"
        hits.run_tool(["hmmsearch", "--cpu", str(cpus), "--noali", "-o", "/dev/null", "--domtblout", str(dt),
                       "-E", "1000", "--domE", "1000", str(hmm), str(fasta)])
        return hits.parse_domtbl(dt)


def union_scores(df, min_score, hmm_len):
    """name -> (best domain score, hmm_from, hmm_to) where the interval is replaced by a surrogate of the same
    length as the union of all domains scoring >= min_score (so select() computes the union coverage).
    The best domain score is kept for the score test. A protein whose domains are all below min_score keeps its best domain."""
    if df.empty:
        return {}
    out = {}
    for pid, g in df.groupby("protein_id"):
        b = g.sort_values("score", ascending=False).iloc[0]
        ok = g[g["score"] >= min_score]
        if ok.empty:
            out[pid] = (float(b.score), int(b.hmm_from), int(b.hmm_to))
        else:
            n = sum(hi - lo + 1 for lo, hi in merge_intervals(zip(ok.hmm_from, ok.hmm_to)))
            out[pid] = (float(b.score), 1, int(n))
    return out


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
    ap.add_argument("--coverage-mode", choices=["best", "union"], default="best",
                    help="best (default): coverage of the best single domain. union: coverage of the union of all domains scoring >= --min-score")
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-table", required=True, type=Path)
    a = ap.parse_args(argv)
    recs = list(read_fasta(a.decoys))
    tmp = a.out_table.parent
    tmp.mkdir(parents=True, exist_ok=True)
    if a.coverage_mode == "union":
        scores = union_scores(domain_hits(a.hmm, a.decoys), a.min_score, hmm_length(a.hmm))
    else:
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
