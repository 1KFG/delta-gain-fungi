#!/usr/bin/env python3
"""Per-genome step calls and pathway presence rules for DHN melanin and L-DOPA melanin.

One row per genome of the BFD species table (species.parquet). Columns:
  searched      1 if the genome is in the searched.chunk*.tsv of the full search with status `searched`, else 0
  prefix_alias  1 for the 5 re-annotated genomes (calls mapped through prefix_aliases.tsv)
  <step>        1 placed call (an `assigned` row of the step's assign_steps output), 0 none, NA genome not searched
  hit_<step>    1 if the step has at least one search-hit row at the min_score (NOT a placed call); NA if not searched
  rule columns  see RULES below (1/0/NA)
A call is placed when at least one protein of the genome has label `assigned` for the step.
Inputs may be plain, .gz or .zst. Protein id prefixes that differ from the table locus tag are mapped with --prefix-aliases.
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

import pandas as pd

# yaml routes of pathways/melanin_dhn.yaml (step_ids of routes pks_thn, pks_hexaketide, pks_heptaketide) and melanin_ldopa.yaml
DHN_ROUTES = {
    "pks_thn": ["thn_synthase", "t4hn_reductase", "scytalone_dehydratase", "t3hn_reductase"],
    "pks_hexaketide": ["hexaketide_pks", "ayg1", "t4hn_reductase", "scytalone_dehydratase", "t3hn_reductase"],
    "pks_heptaketide": ["heptaketide_pks", "ayg1", "t4hn_reductase", "scytalone_dehydratase", "t3hn_reductase"],
}
DHN_PKS = ["thn_synthase", "hexaketide_pks", "heptaketide_pks"]
DHN_STEPS = DHN_PKS + ["ayg1", "t4hn_reductase", "t3hn_reductase", "scytalone_dehydratase", "laccase_abr2", "mco_abr1"]
DHN_OX = ["laccase_abr2", "mco_abr1"]
# pigment-forming L-DOPA-pathway steps (hppd = pyomelanin route and hmga = tyrosine degradation are family present, not counted)
LDOPA_STEPS = ["tyrosinase", "laccase_diphenol", "aspulvinone_e_synthase", "asp_melanin_tyrosinase"]
FAMILY_PRESENT = ["hppd", "hmga"]
TAX = ["PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"]

RULES = {
    # DHN
    "D_any": "any placed call in any of the 9 DHN steps",
    "D_pks": "any placed call in a DHN PKS step (thn_synthase, hexaketide_pks or heptaketide_pks) (headline)",
    "D_route": "all placed steps of at least one yaml PKS route (strict)",
    "D_route_ox": "D_route and a placed laccase_abr2 or mco_abr1",
    "D_hit_any": "search-hit level: at least one hit row in any of the 9 DHN steps (not placed)",
    # L-DOPA
    "L_any": "any placed call in tyrosinase, laccase_diphenol, aspulvinone_e_synthase or asp_melanin_tyrosinase",
    "L_dopa": "tyrosinase or laccase_diphenol placed (the two routes that oxidise L-tyrosine or L-DOPA or diphenols)",
    "L_tyr": "tyrosinase placed",
    "L_nontyr": "laccase_diphenol placed, or both Asp-melanin steps placed (excludes the tyrosinase family)",
    "L_route": "complete route: tyrosinase, or laccase_diphenol, or both Asp-melanin steps (headline)",
    "L_hit_any": "search-hit level: at least one hit row in the 4 pigment-forming L-DOPA steps (not placed)",
}


def read_text(path) -> str:
    path = str(path)
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return Path(path).read_text()


def read_tsv(path):
    lines = [l for l in read_text(path).splitlines() if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))


def load_aliases(path):
    out = {}
    for r in read_tsv(path):
        a, lt = r["alias_prefix"], r["locustag"]
        if out.get(a, lt) != lt:
            raise ValueError(f"alias prefix {a} maps to two locustags")
        out[a] = lt
    return out


def genome_of(protein_id, aliases):
    pre = protein_id.split("_")[0]
    return aliases.get(pre, pre)


def apply_rules(df):
    """Add the rule columns to a frame whose step columns are 0/1 (searched genomes only)."""
    s = df
    df["D_any"] = (s[DHN_STEPS].sum(axis=1) > 0).astype(int)
    df["D_pks"] = (s[DHN_PKS].sum(axis=1) > 0).astype(int)
    route_ok = None
    for steps in DHN_ROUTES.values():
        ok = (s[steps].sum(axis=1) == len(steps))
        route_ok = ok if route_ok is None else (route_ok | ok)
    df["D_route"] = route_ok.astype(int)
    df["D_route_ox"] = (route_ok & (s[DHN_OX].sum(axis=1) > 0)).astype(int)
    for k, steps_ in DHN_ROUTES.items():
        df["route_" + k] = (s[steps_].sum(axis=1) == len(steps_)).astype(int)
    df["D_hit_any"] = (df[[f"hit_{x}" for x in DHN_STEPS]].sum(axis=1) > 0).astype(int)
    df["L_any"] = (s[LDOPA_STEPS].sum(axis=1) > 0).astype(int)
    df["L_dopa"] = (s[["tyrosinase", "laccase_diphenol"]].sum(axis=1) > 0).astype(int)
    df["L_tyr"] = s["tyrosinase"].astype(int)
    df["L_route"] = ((s["tyrosinase"] == 1) | (s["laccase_diphenol"] == 1) |
                     ((s["aspulvinone_e_synthase"] == 1) & (s["asp_melanin_tyrosinase"] == 1))).astype(int)
    df["L_nontyr"] = ((s["laccase_diphenol"] == 1) |
                      ((s["aspulvinone_e_synthase"] == 1) & (s["asp_melanin_tyrosinase"] == 1))).astype(int)
    df["L_hit_any"] = (df[[f"hit_{x}" for x in LDOPA_STEPS]].sum(axis=1) > 0).astype(int)
    return df


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-parquet", required=True)
    ap.add_argument("--steps", required=True, help="steps.tsv (step, pathway, root, assigned_file, role)")
    ap.add_argument("--root", action="append", required=True, help="KEY=results dir, repeatable (keys mel, s5b)")
    ap.add_argument("--search-dir", required=True, help="dir with searched.chunk*.tsv and candidates.chunk*.tsv.zst")
    ap.add_argument("--prefix-aliases", required=True)
    ap.add_argument("--out", required=True, help="per-genome TSV (.zst by extension is not written; plain or .gz)")
    a = ap.parse_args(argv)

    roots = dict(x.split("=", 1) for x in a.root)
    aliases = load_aliases(a.prefix_aliases)
    alias_genomes = set(aliases.values())
    import duckdb
    tax = duckdb.connect().execute(
        f"select LOCUSTAG, PHYLUM, SUBPHYLUM, \"CLASS\", \"ORDER\", FAMILY, GENUS, SPECIES from read_parquet('{a.species_parquet}')").df()
    tax.columns = ["locustag"] + TAX
    assert tax.locustag.is_unique

    sd = Path(a.search_dir)
    srows = [r for f in sorted(glob.glob(str(sd / "searched.chunk*.tsv"))) for r in read_tsv(f)]
    status = {}
    for r in srows:
        status.setdefault(r["locustag"], set()).add(r["status"])
    searched = {g for g, st in status.items() if st == {"searched"}}
    print(f"searched genomes {len(searched)}; in species table {len(searched & set(tax.locustag))}", file=sys.stderr)

    steps = read_tsv(a.steps)
    placed = {}
    for r in steps:
        p = Path(roots[r["root"]]) / r["assigned_file"]
        rows = read_tsv(p)
        gs = {genome_of(x["protein_id"], aliases) for x in rows if x.get("label") == "assigned"}
        placed[r["step"]] = gs
        print(f"{r['step']:24s} assigned rows {sum(1 for x in rows if x.get('label') == 'assigned'):6d} genomes {len(gs):6d} "
              f"not in searched {len(gs - searched)}", file=sys.stderr)
    hits = {}
    for f in sorted(glob.glob(str(sd / "candidates.chunk*.tsv.zst"))):
        for r in read_tsv(f):
            hits.setdefault(r["step_id"], set()).add(genome_of(r["protein_id"], aliases))

    df = tax.copy()
    df["searched"] = df.locustag.isin(searched).astype(int)
    df["prefix_alias"] = df.locustag.isin(alias_genomes).astype(int)
    for r in steps:
        s = r["step"]
        df[s] = df.locustag.isin(placed[s]).astype(int)
        df[f"hit_{s}"] = df.locustag.isin(hits.get(s, set())).astype(int)
    for s in DHN_STEPS + LDOPA_STEPS:
        assert s in df.columns, s
    ok = df.searched == 1
    sub = apply_rules(df.loc[ok].copy())
    rule_cols = list(RULES) + ["route_" + k for k in DHN_ROUTES]
    for c in rule_cols:
        df[c] = pd.NA
        df.loc[ok, c] = sub[c].values
    # not-searched genomes: NA for all step columns, never absence
    stepcols = [r["step"] for r in steps] + [f"hit_{r['step']}" for r in steps]
    for c in stepcols:
        df[c] = df[c].astype("Int64")
        df.loc[~ok, c] = pd.NA
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False, na_rep="NA")
    print(f"wrote {out}: {len(df)} genomes, {int(ok.sum())} searched", file=sys.stderr)
    stray = {s: len(placed[s] - set(df.locustag)) for s in placed}
    print("placed-call genomes not in species table:", stray, file=sys.stderr)


if __name__ == "__main__":
    main()
