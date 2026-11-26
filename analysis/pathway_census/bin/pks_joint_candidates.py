#!/usr/bin/env python3
"""Candidate folders for ONE joint tree of several steps (step 5, PKS steps).

Input: a candidate folder (bin/search_chunk.py format: candidates.chunkNNNNN.PROVISIONAL.tsv.zst, first line `# PROVISIONAL`).
It already holds the rows of each step at that step's cutoff. This script writes:
  --steps-dir   the rows of the listed steps only (same file names, same columns). Used to count hits per step.
  --union-dir   one row per PROTEIN: the row with the highest score over the listed steps, with step_id = --joint-id.
                Used as PC_CAND_DIR for step5_prepare.sbatch and step5_place.sbatch (they cut one domain per protein).
  --scores-out  one row per protein: the best score on each step HMM (empty when the protein has no row of that step),
                the number of steps hit, the step of the best row, and the chunk file.
Scores of different HMMs are on different scales. The best row is the highest raw score, as in step5_place `best_rows`.
Ties go to the step listed first. Every chunk file is written, also when empty. Also writes <union-dir>/union_summary.json.
Compressed and plain input is read. Run with /usr/bin/python3.12.
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path


def read_lines(path):
    path = str(path)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode().splitlines()
    return Path(path).read_text().splitlines()


def write_zst(path, lines):
    subprocess.run(["zstd", "-qf", "-T0", "-o", str(path)], input=("\n".join(lines) + "\n").encode(), check=True)


def split_rows(lines):
    """Return (header line, data lines). The header is the first non-comment line."""
    header, rows = None, []
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        if header is None:
            header = ln
            continue
        rows.append(ln)
    if header is None:
        raise SystemExit("ERROR: no header line")
    return header, rows


def union_rows(rows, cols, steps, joint_id):
    """rows: data lines. Return (steps rows kept, union rows, per-protein best score per step)."""
    ip, isc, ist = cols.index("protein_id"), cols.index("score"), cols.index("step_id")
    kept, best, per = [], {}, {}
    for ln in rows:
        f = ln.split("\t")
        if f[ist] not in steps:
            continue
        kept.append(ln)
        s = float(f[isc])
        pid = f[ip]
        d = per.setdefault(pid, {})
        d[f[ist]] = max(s, d.get(f[ist], float("-inf")))
        cur = best.get(pid)
        # strict > keeps the first row on a tie; steps order breaks ties between steps
        key = (s, -steps.index(f[ist]))
        if cur is None or key > cur[0]:
            best[pid] = (key, f)
    out = []
    for pid, (_, f) in best.items():
        g = list(f)
        g[ist] = joint_id
        out.append("\t".join(g))
    return kept, out, per, {pid: f[ist] for pid, (_, f) in best.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True, type=Path)
    ap.add_argument("--steps", required=True, help="comma list of step ids, in tie-break order")
    ap.add_argument("--joint-id", required=True)
    ap.add_argument("--steps-dir", required=True, type=Path)
    ap.add_argument("--union-dir", required=True, type=Path)
    ap.add_argument("--scores-out", required=True, type=Path)
    a = ap.parse_args(argv)
    steps = a.steps.split(",")
    files = sorted(glob.glob(str(a.in_dir / "candidates.chunk*.PROVISIONAL.tsv.zst")))
    if not files:
        raise SystemExit(f"ERROR: no candidate files in {a.in_dir}")
    a.steps_dir.mkdir(parents=True, exist_ok=True)
    a.union_dir.mkdir(parents=True, exist_ok=True)
    allscores, n_rows, per_chunk = {}, Counter(), {}
    seen_chunk = {}
    for p in files:
        name = Path(p).name
        header, rows = split_rows(read_lines(p))
        cols = header.split("\t")
        for need in ("protein_id", "score", "step_id"):
            if need not in cols:
                raise SystemExit(f"ERROR: column {need} missing in {p}")
        kept, uni, per, beststep = union_rows(rows, cols, steps, a.joint_id)
        for ln in kept:
            n_rows[ln.split("\t")[cols.index("step_id")]] += 1
        for pid, d in per.items():
            if pid in allscores:
                raise SystemExit(f"ERROR: protein {pid} in two chunk files ({seen_chunk[pid]} and {name})")
            allscores[pid] = (d, beststep[pid], name)
            seen_chunk[pid] = name
        write_zst(a.steps_dir / name, ["# PROVISIONAL", header, *kept])
        write_zst(a.union_dir / name, ["# PROVISIONAL", header, *uni])
        per_chunk[name] = {"rows_steps": len(kept), "rows_union": len(uni)}
    with open(a.scores_out, "w") as fh:
        fh.write("\t".join(["protein_id", *[f"score_{s}" for s in steps], "n_steps_hit", "best_step", "chunk_file"]) + "\n")
        for pid in sorted(allscores):
            d, b, ch = allscores[pid]
            fh.write("\t".join([pid, *[f"{d[s]}" if s in d else "" for s in steps], str(len(d)), b, ch]) + "\n")
    nh = Counter(len(v[0]) for v in allscores.values())
    summ = {"steps": steps, "joint_id": a.joint_id, "rows_per_step": dict(n_rows), "proteins_union": len(allscores),
            "proteins_by_n_steps_hit": {str(k): nh[k] for k in sorted(nh)}, "per_chunk": per_chunk}
    (a.union_dir / "union_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps({k: summ[k] for k in ("rows_per_step", "proteins_union", "proteins_by_n_steps_hit")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
