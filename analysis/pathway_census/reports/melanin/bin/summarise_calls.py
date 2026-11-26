#!/usr/bin/env python3
"""Summary tables for the melanin report from the per-genome table of build_genome_calls.py.

Writes TSV files into --out-dir:
  rule_counts.tsv                 genomes positive per rule (denominator: searched genomes)
  species_collapse.tsv            genomes per species, species counts
  sensitivity_four_class.tsv      DHN-only / L-DOPA-only / both / neither for several rule pairs, genome and species level
  four_class_by_<rank>.tsv        primary rule pair by phylum, class, order (genome and species counts)
  step_presence_by_<rank>.tsv     genomes with a placed call per step, by phylum and class
Species level: a species is positive for a pathway when at least one of its searched genomes is positive
(`any`, primary) or when all its searched genomes are positive (`all`, sensitivity). Genomes without proteins
(not searched) are not assessed and are left out of every denominator.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import importlib.util
_spec = importlib.util.spec_from_file_location("bgc", Path(__file__).with_name("build_genome_calls.py"))
bgc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bgc)

PRIMARY = ("D_pks", "L_route")
PAIRS = [
    ("D_pks", "L_route", "HEADLINE (PROPOSED): any DHN PKS placed; complete L-DOPA-pathway route"),
    ("D_any", "L_any", "wide: any placed step of either pathway"),
    ("D_route", "L_route", "strict: all placed steps of a yaml PKS route; complete L-DOPA-pathway route"),
    ("D_pks", "L_nontyr", "headline DHN rule; L-DOPA without the tyrosinase family"),
    ("D_route_ox", "L_route", "DHN route plus an oxidase (laccase_abr2 or mco_abr1); complete L-DOPA route"),
    ("D_hit_any", "L_hit_any", "search-hit level (not placed): any hit row"),
]
NA_LABEL = "(no label)"


def four_class(d, p, l):
    d_, l_ = d[p].astype(int) == 1, d[l].astype(int) == 1
    return pd.Series({"DHN_only": int((d_ & ~l_).sum()), "LDOPA_only": int((~d_ & l_).sum()),
                      "both": int((d_ & l_).sum()), "neither": int((~d_ & ~l_).sum()), "n": len(d)})


def species_table(s, cols):
    """Per species: any / all over searched genomes for each column in cols; n_genomes."""
    g = s.groupby("SPECIES")
    any_ = g[cols].max()
    all_ = g[cols].min()
    n = g.size().rename("n_genomes")
    lin = g[["PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS"]].first()
    return any_, all_, n, lin


def by_rank(s, sp_any, lin, rank, p, l):
    s = s.copy()
    s[rank] = s[rank].fillna(NA_LABEL)
    lin = lin.copy()
    lin[rank] = lin[rank].fillna(NA_LABEL)
    rows = []
    for key, gg in s.groupby(rank):
        ss = sp_any.loc[lin.index[lin[rank] == key]]
        gc, sc = four_class(gg, p, l), four_class(ss, p, l)
        rows.append({rank: key, "PHYLUM": "(mixed)" if key == NA_LABEL and rank != "PHYLUM" else gg.PHYLUM.iloc[0], "n_genomes": len(gg), "n_species": len(ss),
                     **{f"genome_{k}": v for k, v in gc.items() if k != "n"},
                     **{f"species_{k}": v for k, v in sc.items() if k != "n"}})
    return pd.DataFrame(rows).sort_values("n_species", ascending=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    d = pd.read_csv(a.calls, sep="\t", na_values="NA", keep_default_na=False)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    s = d[d.searched == 1].copy()
    rule_cols = list(bgc.RULES) + ["route_" + k for k in bgc.DHN_ROUTES]
    step_cols = [x for x in bgc.DHN_STEPS + bgc.LDOPA_STEPS + bgc.FAMILY_PRESENT]
    for c in rule_cols + step_cols + [f"hit_{x}" for x in step_cols]:
        s[c] = s[c].astype(int)

    # rule counts
    rc = [{"rule": r, "definition": bgc.RULES.get(r, "all placed steps of yaml route " + r[6:]), "genomes_positive": int(s[r].sum()),
           "genomes_searched": len(s)} for r in rule_cols]
    pd.DataFrame(rc).to_csv(out / "rule_counts.tsv", sep="\t", index=False)
    sc = [{"step": x, "role": "family_present" if x in bgc.FAMILY_PRESENT else "placed",
           "genomes_placed_call": int(s[x].sum()), "genomes_hit_row": int(s[f"hit_{x}"].sum()), "genomes_searched": len(s)}
          for x in step_cols]
    pd.DataFrame(sc).to_csv(out / "step_genome_counts.tsv", sep="\t", index=False)

    # species collapse
    any_, all_, n, lin = species_table(s, rule_cols + step_cols)
    coll = {"genomes_in_species_table": len(d), "genomes_searched": len(s), "genomes_not_searched": int((d.searched == 0).sum()),
            "species_in_table": d.SPECIES.nunique(), "species_with_searched_genome": s.SPECIES.nunique(),
            "species_without_searched_genome": d.SPECIES.nunique() - s.SPECIES.nunique(),
            "species_one_genome": int((n == 1).sum()), "species_2_to_9": int(((n >= 2) & (n <= 9)).sum()),
            "species_10_or_more": int((n >= 10).sum()), "max_genomes_per_species": int(n.max()),
            "median_genomes_per_species": float(n.median()), "genomes_in_species_with_10_or_more": int(n[n >= 10].sum()),
            "searched_genomes_no_phylum": int(s.PHYLUM.isna().sum()), "searched_genomes_no_class": int(s.CLASS.isna().sum()),
            "searched_genomes_no_order": int(s.ORDER.isna().sum()),
            "species_with_more_than_one_class_label": int((s.groupby("SPECIES")["CLASS"].nunique(dropna=False) > 1).sum()),
            "species_with_more_than_one_order_label": int((s.groupby("SPECIES")["ORDER"].nunique(dropna=False) > 1).sum()),
            "re_annotated_genomes_searched": int(s.prefix_alias.sum())}
    pd.Series(coll).to_csv(out / "species_collapse.tsv", sep="\t", header=["value"], index_label="quantity")
    top = n.sort_values(ascending=False).head(10).rename_axis("species").reset_index()
    top.to_csv(out / "species_most_genomes.tsv", sep="\t", index=False)

    # sensitivity
    rows = []
    for p, l, note in PAIRS:
        g = four_class(s, p, l)
        sa = four_class(any_, p, l)
        sl = four_class(all_, p, l)
        sx = four_class(any_[n == 1], p, l)
        for lvl, v in (("genome", g), ("species_any", sa), ("species_all", sl), ("species_1genome_only", sx)):
            rows.append({"DHN_rule": p, "LDOPA_rule": l, "note": note, "level": lvl, **v.to_dict()})
    pd.DataFrame(rows).to_csv(out / "sensitivity_four_class.tsv", sep="\t", index=False)

    # excluding the 5 re-annotated genomes (genome level and species level)
    s5 = s[s.prefix_alias == 0]
    a5, _, n5, lin5 = species_table(s5, rule_cols)
    ex = [{"level": "genome_without_5", **four_class(s5, *PRIMARY).to_dict()},
          {"level": "species_any_without_5_genomes", **four_class(a5, *PRIMARY).to_dict()}]
    pd.DataFrame(ex).to_csv(out / "primary_without_5_reannotated.tsv", sep="\t", index=False)
    s.loc[s.prefix_alias == 1, ["locustag", "SPECIES", "CLASS"] + step_cols + rule_cols].to_csv(
        out / "reannotated_5_genomes.tsv", sep="\t", index=False)

    # why the strict route is small: genomes with all non-t3hn steps of a route placed, and how many also have t3hn
    gap = []
    for k, steps_ in bgc.DHN_ROUTES.items():
        rest = [x for x in steps_ if x != "t3hn_reductase"]
        m = s[rest].sum(axis=1) == len(rest)
        gap.append({"route": k, "steps_without_t3hn": ",".join(rest), "genomes_with_those_steps": int(m.sum()),
                    "of_which_with_t3hn_reductase": int((m & (s.t3hn_reductase == 1)).sum())})
    pd.DataFrame(gap).to_csv(out / "route_gap_t3hn.tsv", sep="\t", index=False)
    pks_sum = int(s[bgc.DHN_PKS].sum().sum())
    pd.Series({"sum_of_three_pks_genome_counts": pks_sum, "genomes_with_any_pks": int(s.D_pks.sum())}).to_csv(
        out / "pks_union.tsv", sep="\t", header=["value"], index_label="quantity")
    hit5 = s[s.prefix_alias == 1][["locustag", "SPECIES", "D_hit_any", "L_hit_any"] + [f"hit_{x}" for x in step_cols]]
    hit5.to_csv(out / "reannotated_5_hit_rows.tsv", sep="\t", index=False)

    # by rank, primary pair
    for rank in ("PHYLUM", "CLASS", "ORDER"):
        by_rank(s, any_, lin.join(n), rank, *PRIMARY).to_csv(out / f"four_class_by_{rank.lower()}.tsv", sep="\t", index=False)
    # DHN and L pair at other rule pairs by class, for the sensitivity figure/table
    for p, l, _ in PAIRS[1:3]:
        by_rank(s, any_, lin.join(n), "CLASS", p, l).to_csv(out / f"four_class_by_class.{p}.{l}.tsv", sep="\t", index=False)

    # step presence by phylum and class (genome level)
    for rank in ("PHYLUM", "CLASS"):
        ss = s.copy()
        ss[rank] = ss[rank].fillna(NA_LABEL)
        agg = ss.groupby(rank)[step_cols + ["D_any", "D_pks", "D_route", "L_any", "L_route"]].sum()
        agg.insert(0, "n_genomes", ss.groupby(rank).size())
        agg.insert(1, "PHYLUM", ss.groupby(rank)["PHYLUM"].first())
        agg.sort_values("n_genomes", ascending=False).to_csv(out / f"step_presence_by_{rank.lower()}.tsv", sep="\t")


if __name__ == "__main__":
    main()
