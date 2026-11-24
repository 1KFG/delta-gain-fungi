#!/usr/bin/env python3
"""Alignability screen of candidate tree leaves (outgroups, named leaves) against one or more step HMMs.

For each protein and each HMM: best domain score and best-domain coverage of the HMM (hmmsearch -E 1000 --domE 1000),
the number of domains scoring >= --min-score, and the coverage of the union of those domains. Rule (same as
bin/filter_decoys_for_tree.py): pass_best = best score >= min-score and best-domain coverage >= min-coverage;
pass_union = best score >= min-score and union coverage >= min-coverage. Defaults 10 bits and 0.5 (PROPOSED).
Output: a long TSV (hmm, protein, header_id, best_score, best_cov, n_domains, union_cov, pass_best, pass_union).
Nothing under refs/ is changed. A protein with no hit has empty scores and False.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("filter_decoys_for_tree", Path(__file__).parent / "filter_decoys_for_tree.py")
fd = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fd)

COLS = ["hmm", "protein", "entry", "best_score", "best_cov", "n_domains_ge_min", "union_cov", "pass_best", "pass_union"]


def screen_rows(df, hmm_name, hmm_len, records, min_score=10.0, min_coverage=0.5):
    """df: parsed domtbl (columns protein_id, score, hmm_from, hmm_to); records: [(name, header, seq)].
    Return one row dict per record."""
    rows = []
    for name, hdr, _ in records:
        g = df[df.protein_id == name] if not df.empty else df
        entry = hdr.split()[1].split("|")[-1] if len(hdr.split()) > 1 and "|" in hdr.split()[1] else name
        if g.empty:
            rows.append(dict(hmm=hmm_name, protein=name, entry=entry, best_score="", best_cov="", n_domains_ge_min=0,
                             union_cov="", pass_best=False, pass_union=False))
            continue
        b = g.sort_values("score", ascending=False).iloc[0]
        ok = g[g.score >= min_score]
        bc = (int(b.hmm_to) - int(b.hmm_from) + 1) / hmm_len
        uc = fd.union_coverage(list(zip(ok.hmm_from, ok.hmm_to)), hmm_len) if len(ok) else 0.0
        rows.append(dict(hmm=hmm_name, protein=name, entry=entry, best_score=float(b.score), best_cov=round(bc, 3),
                         n_domains_ge_min=len(ok), union_cov=round(uc, 3),
                         pass_best=bool(b.score >= min_score and bc >= min_coverage),
                         pass_union=bool(b.score >= min_score and uc >= min_coverage)))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--hmm", required=True, nargs="+", type=Path)
    ap.add_argument("--min-score", type=float, default=10.0)
    ap.add_argument("--min-coverage", type=float, default=0.5)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    recs = list(fd.read_fasta(a.fasta))
    allrows = []
    for h in a.hmm:
        df = fd.domain_hits(h, a.fasta)
        allrows += screen_rows(df, h.name.split(".")[0], fd.hmm_length(h), recs, a.min_score, a.min_coverage)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in allrows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")
    print(f"{len(recs)} proteins x {len(a.hmm)} HMMs -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
