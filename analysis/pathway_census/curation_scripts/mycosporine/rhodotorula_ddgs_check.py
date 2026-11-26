#!/usr/bin/env python3.12
"""Best ddgs HMM hit per genome with no score cutoff (owner question 2026-11-25: Rhodotorula, ddgs).

Subcommands:
  search    hmmsearch the ddgs step HMM against the protein files of chosen genera. Writes per_genome.tsv.
  summarize read per_genome.tsv. Write bins.tsv and sensitivity.tsv.

The report threshold of hmmsearch is -E 10 (pathway_census.hits.search_step). A genome with no reported
hit has best_score empty. This is "no hit at E <= 10", not "score 0".
Pure functions (select_genomes, best_hit, bin_scores, sensitivity) are tested in tests/test_rhodotorula_ddgs.py.
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pathway_census import hits  # noqa: E402

BIN_EDGES = [(None, 50, "<50"), (50, 150, "50-150"), (150, 300, "150-300"), (300, None, "300-cutoff"), (None, None, ">=cutoff")]


def read_zst_tsv(path: Path) -> pd.DataFrame:
    raw = subprocess.run(["zstd", "-dc", str(path)], check=True, capture_output=True).stdout.decode()
    return pd.read_csv(io.StringIO(raw), sep="\t", dtype=str, keep_default_na=False)


def select_genomes(genomes: pd.DataFrame, manifest: pd.DataFrame, groups: dict[str, int | None]) -> pd.DataFrame:
    """Genomes of each genus with a protein file. groups maps genus to a max count (None = all).
    When a cap applies, the first n by sorted locustag are taken. Returns locustag, genus, species, proteins_path, group."""
    paths = dict(zip(manifest.locustag, manifest.proteins_path))
    out = []
    for genus, cap in groups.items():
        g = genomes[(genomes.genus == genus) & (genomes.has_proteins == "True")].sort_values("locustag")
        g = g[g.locustag.isin(paths)]
        if cap is not None:
            g = g.head(cap)
        for r in g.itertuples():
            out.append(dict(group=genus, locustag=r.locustag, genus=r.genus, species=r.species, proteins_path=paths[r.locustag]))
    return pd.DataFrame(out)


def hmm_length(hmm: Path) -> int:
    for line in open(hmm):
        if line.startswith("LENG"):
            return int(line.split()[1])
    raise ValueError(f"no LENG in {hmm}")


def best_hit(domtbl: pd.DataFrame, hmm_len: int) -> dict:
    """Best-scoring domain over all proteins. Coverage = hmm positions of that domain / HMM length."""
    if domtbl.empty:
        return dict(best_protein="", best_score="", best_evalue="", hmm_from="", hmm_to="", hmm_coverage="", n_proteins_hit=0)
    b = domtbl.sort_values(["score", "evalue"], ascending=[False, True]).iloc[0]
    cov = (int(b.hmm_to) - int(b.hmm_from) + 1) / hmm_len
    return dict(best_protein=b.protein_id, best_score=float(b.score), best_evalue=float(b.evalue), hmm_from=int(b.hmm_from),
                hmm_to=int(b.hmm_to), hmm_coverage=round(cov, 3), n_proteins_hit=int(domtbl.protein_id.nunique()))


def bin_scores(scores: pd.Series, cutoff: float) -> dict[str, int]:
    """Counts per bin. NaN = no hit at E<=10. Bins: <50, 50-150, 150-300, 300-cutoff, >=cutoff."""
    s = pd.to_numeric(scores, errors="coerce")
    return {"no_hit_E10": int(s.isna().sum()), "<50": int((s < 50).sum()), "50-150": int(((s >= 50) & (s < 150)).sum()),
            "150-300": int(((s >= 150) & (s < 300)).sum()), f"300-{cutoff}": int(((s >= 300) & (s < cutoff)).sum()),
            f">={cutoff}": int((s >= cutoff).sum())}


def sensitivity(scores: pd.Series, cutoff: float, fracs=(0.0, 0.10, 0.25, 0.50)) -> list[dict]:
    """Genomes with best score >= cutoff*(1-f). Counts domain hits, not census calls."""
    s = pd.to_numeric(scores, errors="coerce")
    return [dict(lowered_by=f"{int(f*100)}%", cutoff=round(cutoff * (1 - f), 2), n_genomes_at_or_above=int((s >= cutoff * (1 - f)).sum())) for f in fracs]


def search_one(hmm: Path, hmm_len: int, row: dict, tmp: Path) -> dict:
    dt = tmp / f"{row['locustag']}.domtbl"
    hits.search_step(hmm, Path(row["proteins_path"]), dt, 1)
    res = best_hit(hits.parse_domtbl(dt), hmm_len)
    dt.unlink()
    return {**row, **res}


def cmd_search(a) -> int:
    sc = os.environ.get("SCRATCH")
    if not sc:
        raise SystemExit("SCRATCH is not set")
    genomes = read_zst_tsv(a.genomes)
    manifest = read_zst_tsv(a.chunk_manifest)
    groups = {}
    for spec in a.group:
        g, _, n = spec.partition(":")
        groups[g] = int(n) if n else None
    sel = select_genomes(genomes, manifest, groups)
    aliases = set(pd.read_csv(a.prefix_aliases, sep="\t", dtype=str).locustag)
    overlap = sorted(set(sel.locustag) & aliases)
    print(f"selected {len(sel)} genomes; in prefix_aliases.tsv: {overlap}")
    L = hmm_length(a.hmm)
    tmp = Path(tempfile.mkdtemp(dir=sc))
    rows = sel.to_dict("records")
    with ThreadPoolExecutor(a.cpus) as ex:
        res = list(ex.map(lambda r: search_one(a.hmm, L, r, tmp), rows))
    a.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(res).to_csv(a.out_dir / "per_genome.tsv", sep="\t", index=False)
    (a.out_dir / "selection_notes.txt").write_text(f"hmm={a.hmm}\nhmm_length={L}\nn_selected={len(sel)}\nprefix_alias_overlap={overlap}\n")
    return 0


def cmd_summarize(a) -> int:
    d = pd.read_csv(a.out_dir / "per_genome.tsv", sep="\t")
    brows, srows = [], []
    for grp, g in d.groupby("group", sort=False):
        b = bin_scores(g.best_score, a.cutoff)
        brows.append(dict(group=grp, n_genomes=len(g), **b,
                          n_ge_max_decoy=int((g.best_score >= a.max_decoy).sum()),
                          median_best_score=g.best_score.median(), max_best_score=g.best_score.max()))
        for r in sensitivity(g.best_score, a.cutoff):
            srows.append(dict(group=grp, n_genomes=len(g), **r))
    pd.DataFrame(brows).to_csv(a.out_dir / "bins.tsv", sep="\t", index=False)
    pd.DataFrame(srows).to_csv(a.out_dir / "sensitivity.tsv", sep="\t", index=False)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search")
    s.add_argument("--hmm", type=Path, required=True)
    s.add_argument("--genomes", type=Path, required=True, help="genomes.corrected.tsv.zst")
    s.add_argument("--chunk-manifest", type=Path, required=True)
    s.add_argument("--prefix-aliases", type=Path, required=True)
    s.add_argument("--group", action="append", required=True, help="GENUS or GENUS:N (first N by locustag)")
    s.add_argument("--out-dir", type=Path, required=True)
    s.add_argument("--cpus", type=int, default=4)
    m = sub.add_parser("summarize")
    m.add_argument("--out-dir", type=Path, required=True)
    m.add_argument("--cutoff", type=float, required=True)
    m.add_argument("--max-decoy", type=float, required=True)
    a = ap.parse_args(argv)
    return cmd_search(a) if a.cmd == "search" else cmd_summarize(a)


if __name__ == "__main__":
    sys.exit(main())
