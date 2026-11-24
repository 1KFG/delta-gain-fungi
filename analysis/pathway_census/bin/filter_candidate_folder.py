#!/usr/bin/env python3
"""Write a candidate folder that holds only the rows of one step at or above a score floor.

Input: a candidate folder with candidates.chunkNNNNN[.PROVISIONAL].tsv.zst (bin/search_chunk.py format).
Output: the same file names in --out-dir, first line `# PROVISIONAL`, same header and columns, only rows with
step_id == --step and score >= --floor. No new search is run: the full search already holds every row of each step at
that step's cutoff. Every chunk file is written, also when it has no row, so that readers find the same file set.
Also writes <out-dir>/filter_summary.json with the counts (rows kept per chunk, floor, step).
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from pathlib import Path


def read_zst_lines(path):
    path = str(path)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode().splitlines()
    return Path(path).read_text().splitlines()


def filter_lines(lines, step, floor):
    """lines: decoded file lines. Return (header, kept data lines, n rows of the step, n rows in the file).
    The header is the first non-comment line. A malformed row (too few columns) stops the run."""
    header, kept, n_step, n_all = None, [], 0, 0
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        if header is None:
            header = ln
            cols = header.split("\t")
            for need in ("score", "step_id"):
                if need not in cols:
                    raise SystemExit(f"ERROR: column {need} missing in the header")
            i_score, i_step = cols.index("score"), cols.index("step_id")
            continue
        f = ln.split("\t")
        if len(f) != len(cols):
            raise SystemExit(f"ERROR: row with {len(f)} columns, header has {len(cols)}: {ln[:80]}")
        n_all += 1
        if f[i_step] == step:
            n_step += 1
            if float(f[i_score]) >= floor:
                kept.append(ln)
    return header, kept, n_step, n_all


def run(in_dir, out_dir, step, floor):
    files = sorted(glob.glob(str(Path(in_dir) / "candidates.chunk*.PROVISIONAL.tsv.zst")))
    if not files:
        raise SystemExit(f"ERROR: no candidate files in {in_dir}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    per = {}
    for p in files:
        header, kept, n_step, n_all = filter_lines(read_zst_lines(p), step, floor)
        text = "# PROVISIONAL\n" + (header or "") + "\n" + "".join(k + "\n" for k in kept)
        dest = out_dir / Path(p).name
        subprocess.run(["zstd", "-q", "-T0", "-f", "-o", str(dest)], input=text.encode(), check=True)
        per[Path(p).name] = {"rows_of_step": n_step, "rows_kept": len(kept), "rows_in_file": n_all}
    summ = {"step": step, "floor": floor, "input": str(in_dir),
            "rows_kept": sum(v["rows_kept"] for v in per.values()),
            "rows_of_step": sum(v["rows_of_step"] for v in per.values()), "per_file": per}
    (out_dir / "filter_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--step", required=True)
    ap.add_argument("--floor", required=True, type=float)
    a = ap.parse_args(argv)
    s = run(a.in_dir, a.out_dir, a.step, a.floor)
    print(f"step {s['step']} floor {s['floor']}: kept {s['rows_kept']} of {s['rows_of_step']} rows of the step")
    return 0


if __name__ == "__main__":
    sys.exit(main())
