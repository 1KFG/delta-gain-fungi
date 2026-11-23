#!/usr/bin/env python3
"""Keep only the decoys that can be aligned with the step domain, for use as tree leaves in step 5.

A decoy with no detectable homology to the step HMM cannot be aligned. In a maximum likelihood tree it lands on an
arbitrary branch, often inside the clade of the references. Such a decoy still serves step 2 (cutoff calibration) and
the absence check; it does not belong in the tree. Rule (PROPOSED): keep a decoy when the best domain score of the step
HMM against it, searched at E <= 1000, is at least --min-score bits. A decoy without any hit is dropped.
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


def best_scores(hmm: Path, fasta: Path, cpus: int = 2) -> dict[str, float]:
    with tempfile.TemporaryDirectory() as td:
        dt = Path(td) / "d.domtbl"
        # E <= 1000 and domain E <= 1000, so weak hits of a few bits are still reported
        hits.run_tool(["hmmsearch", "--cpu", str(cpus), "--noali", "-o", "/dev/null", "--domtblout", str(dt),
                       "-E", "1000", "--domE", "1000", str(hmm), str(fasta)])
        df = hits.parse_domtbl(dt)
    if df.empty:
        return {}
    return df.groupby("protein_id").score.max().to_dict()


def select(records, scores, min_score):
    """records: list of (name, header, seq). Return (kept, table rows (name, score or None, kept))."""
    kept, table = [], []
    for name, hdr, seq in records:
        sc = scores.get(name)
        keep = sc is not None and sc >= min_score
        table.append((name, sc, keep))
        if keep:
            kept.append((name, hdr, seq))
    return kept, table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--decoys", required=True, type=Path)
    ap.add_argument("--hmm", required=True, type=Path)
    ap.add_argument("--min-score", required=True, type=float)
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-table", required=True, type=Path)
    a = ap.parse_args(argv)
    recs = list(read_fasta(a.decoys))
    tmp = a.out_table.parent
    tmp.mkdir(parents=True, exist_ok=True)
    scores = best_scores(a.hmm, a.decoys)
    kept, table = select(recs, scores, a.min_score)
    with open(a.out_fasta, "w") as fh:
        for name, hdr, seq in kept:
            fh.write(f">{hdr}\n{seq}\n")
    with open(a.out_table, "w") as fh:
        fh.write("decoy\tbest_domain_score\tkept_in_tree\n")
        for name, sc, keep in table:
            fh.write(f"{name}\t{'' if sc is None else sc}\t{keep}\n")
    print(f"kept {len(kept)} of {len(recs)} decoys (min score {a.min_score})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
