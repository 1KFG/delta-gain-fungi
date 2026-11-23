#!/usr/bin/env python3
"""Compare, per class, three call sets for one step: search hits at the rule cutoff, search hits at a lower cutoff, and proteins placed in the target clade.

Writes <out-dir>/by_class.tsv (class, n genomes, share with each call set) and <out-dir>/score_bands.tsv (for each score band of the
lower-cutoff candidates: how many fall in the target clade, on the decoy side and unassigned). Genomes are the searched genomes of the
lower-cutoff run. Reads results of bin/search_chunk.py, bin/step5_place.py and bin/assign_steps.py. Nothing here is a call.
"""
from __future__ import annotations

import argparse
import csv
import glob
import io
import subprocess
from pathlib import Path


def read_tsv(path):
    path = str(path)
    if path.endswith(".zst"):
        raw = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    else:
        raw = Path(path).read_text()
    return list(csv.DictReader(io.StringIO("\n".join(l for l in raw.splitlines() if not l.startswith("#"))), delimiter="\t"))


def genomes_with_hit(search_dir, step, cutoff):
    out = set()
    for f in sorted(glob.glob(str(Path(search_dir) / "candidates.chunk*.tsv.zst"))):
        out |= {r["protein_id"].split("_")[0] for r in read_tsv(f) if r["step_id"] == step and float(r["score"]) >= cutoff}
    return out


def genomes_placed(assigned_tsv):
    return {r["protein_id"].split("_")[0] for r in read_tsv(assigned_tsv) if r["label"] == "assigned"}


def by_class(genomes, classes, sets):
    """genomes: list of locustags. classes: locustag -> class. sets: ordered {name: set of genomes}. Return rows (class, n, {name: share})."""
    members = {}
    for g in genomes:
        members.setdefault(classes.get(g) or "(no label)", []).append(g)
    rows = [(c, len(gs), {k: sum(g in s for g in gs) / len(gs) for k, s in sets.items()}) for c, gs in members.items()]
    rows.sort(key=lambda r: -r[1])
    return rows


def band_table(scores, calls, edges):
    """scores: protein -> score. calls: protein -> 'target' | 'decoy' | ''. Rows (low, high, target, decoy, unassigned)."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        c = {"target": 0, "decoy": 0, "": 0}
        for p, s in scores.items():
            if lo <= s < hi:
                c[calls.get(p, "") if calls.get(p, "") in c else ""] += 1
        rows.append((lo, hi, c["target"], c["decoy"], c[""]))
    return rows


def main(argv=None):
    import duckdb
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", required=True)
    ap.add_argument("--rule-search-dir", required=True, type=Path)
    ap.add_argument("--low-search-dir", required=True, type=Path)
    ap.add_argument("--rule-cutoff", required=True, type=float)
    ap.add_argument("--low-cutoff", required=True, type=float)
    ap.add_argument("--assigned", required=True, type=Path)
    ap.add_argument("--placements-detail", required=True, type=Path)
    ap.add_argument("--candidates-step", required=True, type=Path)
    ap.add_argument("--species", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args(argv)
    genomes = []
    for f in sorted(glob.glob(str(a.low_search_dir / "searched.chunk*.tsv"))):
        genomes += [r["locustag"] for r in read_tsv(f) if r["status"] == "searched"]
    classes = {r[0]: r[1] for r in duckdb.connect().execute(f"select LOCUSTAG, \"CLASS\" from read_parquet('{a.species}')").fetchall()}
    sets = {"rule_hits": genomes_with_hit(a.rule_search_dir, a.step, a.rule_cutoff),
            "low_hits": genomes_with_hit(a.low_search_dir, a.step, a.low_cutoff),
            "placed": genomes_placed(a.assigned)}
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rows = by_class(genomes, classes, sets)
    with open(a.out_dir / "by_class.tsv", "w") as fh:
        fh.write("class\tn_genomes\t" + "\t".join(f"share_{k}" for k in sets) + "\n")
        for c, n, d in rows:
            fh.write(f"{c}\t{n}\t" + "\t".join(f"{d[k]:.4f}" for k in sets) + "\n")
        fh.write(f"ALL\t{len(genomes)}\t" + "\t".join(f"{len(sets[k] & set(genomes)) / len(genomes):.4f}" for k in sets) + "\n")
    scores = {r["protein_id"]: float(r["score"]) for r in read_tsv(a.candidates_step)}
    calls = {r["protein_id"]: r["clade_call"] for r in read_tsv(a.placements_detail)}
    with open(a.out_dir / "score_bands.tsv", "w") as fh:
        fh.write("low\thigh\ttarget\tdecoy_side\tunassigned\n")
        for r in band_table(scores, calls, [a.low_cutoff, 200, 250, 300, a.rule_cutoff, 450, 10000]):
            fh.write("\t".join(f"{x:g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    print(f"wrote {a.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
