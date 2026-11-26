#!/usr/bin/env python3
"""Write the per-step summary table of a pathway census report and insert it into REPORT.md.

The table lists, per step: reference count, decoy count, the cutoff of the rule, the cutoff used, candidate rows and
genomes with a candidate. All numbers come from census files (see fig_census_search.py). The text between the markers
<!-- step-table:begin --> and <!-- step-table:end --> in --report is replaced.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fig_census_search", HERE.parent / "figures" / "bin" / "fig_census_search.py")
fc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fc)

BEGIN, END = "<!-- step-table:begin -->", "<!-- step-table:end -->"


def count_fasta(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for line in open(path) if line.startswith(">"))


def build_table(tab, counts, refs_dir: Path) -> str:
    lines = ["| Step | References | Decoys | Cutoff of the rule | Cutoff used | Candidate rows | Genomes with a row |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for r in sorted(tab, key=lambda r: r["step"]):
        n_dec = count_fasta(refs_dir / "decoys" / f"{r['step']}.faa")
        rows, gen = counts.get(r["step"], (0, 0))
        used = f"{r['used']:g}" + ("" if r["used"] == r["rule"] else " (owner override)")
        lines.append(f"| `{r['step']}` | {r['n_refs']} | {n_dec} | {r['rule']:g} | {used} | {rows:,} | {gen:,} |")
    return "\n".join(lines)


def insert(report_text: str, table: str) -> str:
    if BEGIN not in report_text or END not in report_text:
        raise ValueError("markers not found in report")
    head, rest = report_text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    return f"{head}{BEGIN}\n{table}\n{END}{tail}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step2-dir", required=True, type=Path)
    ap.add_argument("--search-dir", required=True, type=Path)
    ap.add_argument("--min-scores", required=True, type=Path)
    ap.add_argument("--refs-dir", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--steps", nargs="*", default=None, help="keep only these steps")
    ap.add_argument("--tag", default="PROVISIONAL")
    a = ap.parse_args(argv)
    sfx = f".{a.tag}" if a.tag else ""
    tab = fc.step_table(a.step2_dir / f"step_scores{sfx}.tsv", a.min_scores, a.steps)
    table = build_table(tab, fc.count_candidates(a.search_dir), a.refs_dir)
    a.report.write_text(insert(a.report.read_text(), table))
    print("table written into", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
