#!/usr/bin/env python3.12
"""Merge bin/absence_check.py chunk files into one absence table.

Usage: merge_absence.py --out results/carotenoid_absence_summary.tsv CHUNK.tsv [CHUNK.tsv ...]

Rules
- A chunk is provisional if its name contains ".PROVISIONAL." OR one of its leading "#" lines
  contains "PROVISIONAL". Either signal is enough.
- Mixed provisional and non-provisional chunks are refused (exit non-zero, nothing written).
- Provisional chunks: the output name carries ".PROVISIONAL" (X.tsv -> X.PROVISIONAL.tsv) and the
  file starts with the chunks' "# PROVISIONAL" line.
- All chunks must have the same leading "#" lines (for example the same min_bitscore and
  min_query_cover) and the same column header. Otherwise the merge is refused.
- A (locustag, step_id) pair in more than one chunk is refused (overlapping chunks).
- Only LEADING "#" lines are header lines. A "#" inside a data field is kept.
Input: plain or .gz TSV. Output: plain TSV (header lines, column header once, rows in input order).
"""
import argparse
import gzip
import sys
from pathlib import Path

TAG = "PROVISIONAL"


def _read_lines(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return fh.read().splitlines()
    return path.read_text().splitlines()


def _split(path: Path):
    lines = _read_lines(path)
    n = 0
    while n < len(lines) and lines[n].startswith("#"):
        n += 1
    if n >= len(lines):
        raise ValueError(f"{path}: no column header line")
    return lines[:n], lines[n], [ln for ln in lines[n + 1:] if ln.strip()]


def _is_provisional(path: Path, header_lines: list[str]) -> bool:
    return f".{TAG}." in path.name or any(TAG in h for h in header_lines)


def merge(inputs: list[Path], out: Path) -> Path:
    """Merge chunk files. Return the path written. Raise ValueError on a refused merge."""
    if not inputs:
        raise ValueError("no input chunk files")
    parsed = [(p, *_split(p)) for p in inputs]
    prov = {p: _is_provisional(p, h) for p, h, _c, _r in parsed}
    if len(set(prov.values())) > 1:
        yes = [str(p) for p, v in prov.items() if v]
        no = [str(p) for p, v in prov.items() if not v]
        raise ValueError(f"mixed provisional and non-provisional inputs. Provisional: {yes[:3]}; "
                         f"not provisional: {no[:3]}")
    p0, h0, c0, _ = parsed[0]
    for p, h, c, _ in parsed[1:]:
        if h != h0:
            raise ValueError(f"header lines differ between {p0} and {p}: {h0} vs {h}")
        if c != c0:
            raise ValueError(f"column header differs between {p0} and {p}")
    cols = c0.split("\t")
    if cols[:2] != ["locustag", "step_id"]:
        raise ValueError(f"expected columns starting locustag, step_id; got {cols}")
    seen: dict[tuple[str, str], Path] = {}
    rows: list[str] = []
    for p, _h, _c, body in parsed:
        for ln in body:
            key = tuple(ln.split("\t")[:2])
            if key in seen:
                raise ValueError(f"duplicate (locustag, step_id) {key} in {seen[key]} and {p}")
            seen[key] = p
            rows.append(ln)
    is_prov = next(iter(prov.values()))
    header = list(h0)
    if is_prov:
        out = out.with_name(f"{out.stem}.{TAG}{out.suffix or '.tsv'}")
        if not header or not header[0].startswith(f"# {TAG}"):
            # the tag came from the file name only: put a PROVISIONAL line first
            header = [f"# {TAG}: merged from chunk files tagged {TAG}; not user-reviewed"] + \
                [h for h in header if TAG not in h]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(f"{x}\n" for x in header + [c0] + rows))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True,
                    help="canonical output name; '.PROVISIONAL' is inserted for provisional input")
    ap.add_argument("chunks", type=Path, nargs="+")
    a = ap.parse_args()
    try:
        out = merge(a.chunks, a.out)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")
    print(f"merged {len(a.chunks)} chunk files -> {out}")


if __name__ == "__main__":
    main()
