#!/usr/bin/env python3
"""Build the genus-representative x step calls table for the gain and loss analysis.

Value per genome and step:
  1   the genome was searched and the step has a call
  0   the genome was searched and the step has no call
  NA  the genome was not searched (not in the searched.* files, or a status other than `searched`,
      or no protein file). Missing data is never absence.

A call is
  - source `hits`: at least one row for the step in candidates.chunk*.tsv.zst (a search hit, NOT placed), or
  - source = path to an assigned.tsv from bin/assign_steps.py: at least one row with label `assigned` (placed).
The config TSV (columns step, source) chooses the source per step. A step switches from hits to placed by editing
the config; no code change is needed. Relative source paths resolve against --root.

Outputs: calls.tsv (one row per genus representative) and calls.summary.tsv (one row per step).
Inputs may be plain, .gz or .zst.
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import io
import subprocess
import sys
from pathlib import Path

TAX_COLS = ["PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS"]


def read_text(path) -> str:
    path = str(path)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return Path(path).read_text()


def read_tsv(path):
    """Rows as dicts; lines starting with `#` are skipped."""
    lines = [l for l in read_text(path).splitlines() if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))


def genome_of(protein_id: str) -> str:
    return protein_id.split("_")[0]


def searched_genomes(rows):
    """Set of locustags whose every searched.* row has status `searched` (conservative for duplicates)."""
    status = {}
    for r in rows:
        status.setdefault(r["locustag"], set()).add(r["status"])
    return {g for g, s in status.items() if s == {"searched"}}


def placed_genomes(rows):
    """Genomes with at least one row labelled `assigned`."""
    return {genome_of(r["protein_id"]) for r in rows if r.get("label") == "assigned"}


def hit_genomes(rows_by_file, step):
    """Genomes with at least one hit row for `step`."""
    return {genome_of(r["protein_id"]) for rows in rows_by_file for r in rows if r["step_id"] == step}


def build_calls(reps, searched, step_sources, step_genomes):
    """reps: list of dicts with LOCUSTAG, pep_path and taxonomy. step_sources: step -> label. step_genomes: step -> set.

    Returns (rows, summary). A representative with an empty pep_path or outside `searched` is NA for every step.
    """
    steps = list(step_sources)
    rows, seen = [], set()
    for r in reps:
        g = r["LOCUSTAG"]
        if g in seen:  # duplicate representative row: keep the first
            continue
        seen.add(g)
        ok = g in searched and bool((r.get("pep_path") or "").strip())
        row = {"locustag": g, **{c: r.get(c, "") for c in TAX_COLS}, "searched": int(ok)}
        for s in steps:
            row[s] = ("1" if g in step_genomes[s] else "0") if ok else "NA"
        rows.append(row)
    summary = []
    rep_ids = {r["locustag"] for r in rows}
    for s in steps:
        vals = [r[s] for r in rows]
        summary.append({"step": s, "source": step_sources[s], "n_reps": len(rows), "n_present": vals.count("1"),
                        "n_absent": vals.count("0"), "n_NA": vals.count("NA"),
                        "n_call_genomes_not_in_reps": len(step_genomes[s] - rep_ids),
                        "n_rep_calls_but_NA": sum(1 for r in rows if r[s] == "NA" and r["locustag"] in step_genomes[s])})
    return rows, summary


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def load_taxonomy(path):
    import pandas as pd
    t = pd.read_parquet(path).drop_duplicates("LOCUSTAG").set_index("LOCUSTAG")
    return t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", required=True, help="genus_representatives.tsv (LOCUSTAG, GENUS, CLASS, busco_complete, pep_path)")
    ap.add_argument("--taxonomy", required=True, help="species.parquet (LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, ORDER, FAMILY, GENUS)")
    ap.add_argument("--search-dir", required=True, help="dir with searched.chunk*.tsv and candidates.chunk*.tsv.zst")
    ap.add_argument("--config", required=True, help="TSV: step, source (`hits` or a path to assigned.tsv)")
    ap.add_argument("--root", default=".", help="base for relative source paths in the config")
    ap.add_argument("--out", required=True, help="calls.tsv path; calls.summary.tsv is written beside it")
    a = ap.parse_args(argv)

    reps = read_tsv(a.reps)
    tax = load_taxonomy(a.taxonomy)
    for r in reps:
        if r["LOCUSTAG"] in tax.index:
            t = tax.loc[r["LOCUSTAG"]]
            for c in TAX_COLS:
                v = t[c]
                r[c] = "" if v is None or v != v else str(v)
        else:
            r["GENUS"] = r.get("GENUS", "")
    sd = Path(a.search_dir)
    srows = [r for f in sorted(glob.glob(str(sd / "searched.chunk*.tsv"))) for r in read_tsv(f)]
    searched = searched_genomes(srows)
    cfg = read_tsv(a.config)
    cand = None
    sources, step_genomes = {}, {}
    for c in cfg:
        step, src = c["step"], c["source"].strip()
        if src == "hits":
            if cand is None:
                cand = [read_tsv(f) for f in sorted(glob.glob(str(sd / "candidates.chunk*.tsv.zst")))]
            step_genomes[step] = hit_genomes(cand, step)
            sources[step] = "hits"
        else:
            p = Path(src) if Path(src).is_absolute() else Path(a.root) / src
            step_genomes[step] = placed_genomes(read_tsv(p))
            sources[step] = "placed"
    rows, summary = build_calls(reps, searched, sources, step_genomes)
    steps = list(sources)
    cols = ["locustag", *TAX_COLS, "searched", *steps]
    for r in rows:
        for s in steps:
            r[f"{s}_source"] = sources[s]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(out, rows, cols + [f"{s}_source" for s in steps])
    write_tsv(out.with_name("calls.summary.tsv"), summary,
              ["step", "source", "n_reps", "n_present", "n_absent", "n_NA", "n_call_genomes_not_in_reps", "n_rep_calls_but_NA"])
    for s in summary:
        print(*(f"{k}={v}" for k, v in s.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
