#!/usr/bin/env python3
"""Cut tree-leaf proteins (named members, outgroups, decoys) to the step HMM domain, as step 5 cuts the candidates.

Why: a leaf that is longer than the domain (a fusion, or a protein with extra regions) pulls long unaligned regions into
the MSA. Candidates are cut to the best domain (ali_from..ali_to) plus a pad (bin/step5_place.py `extract`, pad 20). The
leaves get the same cut so that all sequences in the tree cover the same region. Rule (PROPOSED): best domain by score of
hmmsearch -E 1000 --domE 1000; start = max(1, ali_from - pad), end = min(length, ali_to + pad). A protein with no hit is
kept whole and listed with an empty domain (it should not occur for leaves that passed the screen).
Output: FASTA with the same ids and headers, and a TSV (id, length, ali_from, ali_to, cut_from, cut_to, cut_length).
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location("filter_decoys_for_tree", Path(__file__).parent / "filter_decoys_for_tree.py")
fd = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fd)


def cut_bounds(length: int, ali_from: int, ali_to: int, pad: int = 20):
    """1-based inclusive (start, end) of the padded domain, clamped to the protein."""
    if ali_from < 1 or ali_to < ali_from:
        raise ValueError(f"bad alignment coordinates {ali_from}-{ali_to}")
    start, end = max(1, ali_from - pad), min(length, ali_to + pad)
    if start > end:
        raise ValueError(f"alignment {ali_from}-{ali_to} is outside a protein of length {length}")
    return start, end


def best_domain_coords(df):
    """DataFrame of domain hits -> {protein_id: (ali_from, ali_to)} of the highest-scoring domain."""
    if df.empty:
        return {}
    b = df.sort_values("score", ascending=False).drop_duplicates("protein_id")
    return {r.protein_id: (int(r.ali_from), int(r.ali_to)) for r in b.itertuples()}


def cut_records(records, coords, pad=20):
    """records [(name, header, seq)] -> (cut records, table rows)."""
    out, table = [], []
    for name, hdr, seq in records:
        c = coords.get(name)
        if c is None:
            out.append((name, hdr, seq))
            table.append((name, len(seq), "", "", 1, len(seq), len(seq)))
            continue
        s, e = cut_bounds(len(seq), c[0], c[1], pad)
        out.append((name, hdr, seq[s - 1:e]))
        table.append((name, len(seq), c[0], c[1], s, e, e - s + 1))
    return out, table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--hmm", required=True, type=Path)
    ap.add_argument("--pad", type=int, default=20)
    ap.add_argument("--out-fasta", required=True, type=Path)
    ap.add_argument("--out-table", required=True, type=Path)
    a = ap.parse_args(argv)
    recs = list(fd.read_fasta(a.fasta))
    if not recs:
        a.out_fasta.write_text("")
        a.out_table.write_text("id\tlength\tali_from\tali_to\tcut_from\tcut_to\tcut_length\n")
        return 0
    coords = best_domain_coords(fd.domain_hits(a.hmm, a.fasta))
    cut, table = cut_records(recs, coords, a.pad)
    with open(a.out_fasta, "w") as fh:
        for n, h, s in cut:
            fh.write(f">{h}\n{s}\n")
    with open(a.out_table, "w") as fh:
        fh.write("id\tlength\tali_from\tali_to\tcut_from\tcut_to\tcut_length\n")
        for r in table:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"cut {sum(1 for r in table if r[2] != '')} of {len(table)} to the domain; {sum(1 for r in table if r[2] == '')} kept whole")
    return 0


if __name__ == "__main__":
    sys.exit(main())
