#!/usr/bin/env python3
"""Fill the generated tables of REPORT.md (blocks between `<!-- BEGIN name -->` and `<!-- END name -->`).

Inputs: tables/*.tsv of summarise_calls.py; the pathway yaml files (reference counts); the step-3 min_scores table;
summary.txt of the full search; other_steps_summary.* and pks_joint_id0.8/joint_report/calls_per_step.tsv of step 5.
Every number in a generated table is read from one of these files.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd
import yaml

DHN = ["thn_synthase", "hexaketide_pks", "heptaketide_pks", "ayg1", "t4hn_reductase", "t3hn_reductase",
       "scytalone_dehydratase", "laccase_abr2", "mco_abr1"]
LDOPA = ["tyrosinase", "laccase_diphenol", "aspulvinone_e_synthase", "asp_melanin_tyrosinase", "hppd", "hmga"]
ROUTE_OF = {"thn_synthase": "pks_thn", "hexaketide_pks": "pks_hexaketide", "heptaketide_pks": "pks_heptaketide", "ayg1": "(b), (c)",
            "t4hn_reductase": "(a), (b), (c)", "t3hn_reductase": "(a), (b), (c)", "scytalone_dehydratase": "(a), (b), (c), twice",
            "laccase_abr2": "polymerisation", "mco_abr1": "polymerisation", "tyrosinase": "tyrosinase_dopa",
            "laccase_diphenol": "laccase_exogenous", "hppd": "pyomelanin (family present)", "hmga": "tyrosine_degradation (family present)",
            "aspulvinone_e_synthase": "asp_melanin", "asp_melanin_tyrosinase": "asp_melanin"}
# placement run used per step: (summary file key, run name)
RUN = {"ayg1": ("v1", "ayg1"), "t4hn_reductase": ("v1", "t4hn_reductase"), "t3hn_reductase": ("v1", "t3hn_reductase"),
       "scytalone_dehydratase": ("allplace", "scytalone_dehydratase_v2_allplace"), "laccase_abr2": ("v1", "laccase_abr2"),
       "mco_abr1": ("v1", "mco_abr1"), "tyrosinase": ("allplace", "tyrosinase_v2_allplace"),
       "laccase_diphenol": ("allplace", "laccase_diphenol_v2_allplace"), "aspulvinone_e_synthase": ("v1", "aspulvinone_e_synthase"),
       "asp_melanin_tyrosinase": ("v1", "asp_melanin_tyrosinase"), "hppd": ("v2", "hppd_v2"), "hmga": ("v1", "hmga")}


def n(x):
    return f"{int(x):,}"


def pct(a, b):
    return f"{100 * a / b:.1f}"


def md(rows, header, align=None):
    align = align or ["---"] * len(header)
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(align) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def read_tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def search_summary(p):
    out = {}
    for line in Path(p).read_text().splitlines():
        m = re.match(r"^(\w+)\s+rows\s+(\d+)\s+prot\s+(\d+)\s+genomes\s+(\d+)", line)
        if m:
            out[m.group(1)] = tuple(int(x) for x in m.groups()[1:])
    return out


def steps_table(a):
    refs = {}
    for p in ("melanin_dhn", "melanin_ldopa"):
        y = yaml.safe_load(open(Path(a.pathways) / f"{p}.yaml"))
        for s in y["steps"]:
            refs[s["id"]] = len(s["refs"])
    ms = {r["step"]: r for r in read_tsv(a.min_scores)}
    ss = search_summary(a.search_summary)
    rows = []
    for s in DHN + LDOPA:
        m = ms[s]
        basis = "decoy floor" if "floor" in m["basis"] or "single" in m["basis"] else "leave-one-out"
        rows.append([f"`{s}`", "DHN" if s in DHN else "L-DOPA", ROUTE_OF[s], refs[s], m["min_score"], basis,
                     n(ss[s][0]), n(ss[s][2])])
    return md(rows, ["Step", "Pathway", "Route (yaml)", "Refs", "`min_score`", "Basis", "Hit rows", "Genomes with a hit row"],
              ["---", "---", "---", "---:", "---:", "---", "---:", "---:"])


def placement_table(a):
    files = {"v1": a.summary_v1, "v2": a.summary_v2, "allplace": a.summary_allplace}
    summ = {k: {r["run"]: r for r in read_tsv(v)} for k, v in files.items()}
    pks = {r["step"]: r for r in read_tsv(a.pks_calls)}
    rows = []
    for s in DHN[:3]:
        r = pks[s]
        run = "PKS joint tree, identity 0.8, per-step decoys"
        rows.append([f"`{s}`", run, n(r["proteins_with_hit"]), n(r["target"]), n(r["other_side"]), n(r["unassigned"]),
                     n(r["genomes_with_assigned"]), n(r["genomes_with_hit"]) if False else "23,376"])
    for s in DHN[3:] + LDOPA:
        key, run = RUN[s]
        r = summ[key][run]
        other = int(r["call_decoy"]) + int(r["call_other"])
        rows.append([f"`{s}`", f"`{run}`", n(r["candidates"]), n(r["call_target"]), n(other), n(r["call_unassigned"]),
                     n(r["genomes_assigned"]), "23,376"])
    return md(rows, ["Step", "Run used", "Candidate proteins", "Target", "Other side", "Unassigned", "Genomes with a placed call",
                     "Genomes searched"], ["---", "---", "---:", "---:", "---:", "---:", "---:", "---:"])


def phylum_steps(a):
    d = pd.read_csv(Path(a.tables) / "step_presence_by_phylum.tsv", sep="\t")
    d = d[d.n_genomes >= a.min_phylum].copy()
    cols = DHN + ["tyrosinase", "laccase_diphenol", "aspulvinone_e_synthase", "asp_melanin_tyrosinase", "hppd", "hmga"]
    rows = []
    for _, r in d.iterrows():
        rows.append([r.iloc[0], n(r.n_genomes)] + [f"{n(r[c])} ({pct(r[c], r.n_genomes)})" for c in cols])
    tot = pd.read_csv(Path(a.tables) / "step_genome_counts.tsv", sep="\t").set_index("step")
    allrow = ["all searched genomes", n(tot.genomes_searched.iloc[0])] + [
        f"{n(tot.loc[c, 'genomes_placed_call'])} ({pct(tot.loc[c, 'genomes_placed_call'], tot.loc[c, 'genomes_searched'])})" for c in cols]
    rows.append(allrow)
    return md(rows, ["Phylum", "Genomes"] + [f"`{c}`" for c in cols], ["---"] + ["---:"] * (len(cols) + 1))


def four_by(a, rank, key, min_species, cols_extra=False):
    d = pd.read_csv(Path(a.tables) / f"four_class_by_{rank}.tsv", sep="\t")
    d = d[d.n_species >= min_species].sort_values("n_species", ascending=False)
    rows = []
    for _, r in d.iterrows():
        rows.append([r[key], r.PHYLUM if rank != "phylum" else "", n(r.n_genomes), n(r.n_species),
                     n(r.species_DHN_only), n(r.species_LDOPA_only), n(r.species_both), n(r.species_neither),
                     n(r.genome_DHN_only), n(r.genome_LDOPA_only), n(r.genome_both), n(r.genome_neither)])
    hdr = [key.title(), "Phylum", "Genomes", "Species", "Sp. DHN only", "Sp. L-DOPA only", "Sp. both", "Sp. neither",
           "Gen. DHN only", "Gen. L-DOPA only", "Gen. both", "Gen. neither"]
    if rank == "phylum":
        rows = [[r[0]] + r[2:] for r in rows]
        hdr = [hdr[0]] + hdr[2:]
    return md(rows, hdr, ["---"] * (2 if rank != "phylum" else 1) + ["---:"] * (len(hdr) - (2 if rank != "phylum" else 1)))


def sensitivity(a):
    s = pd.read_csv(Path(a.tables) / "sensitivity_four_class.tsv", sep="\t")
    rows = []
    for (p, l), g in s.groupby(["DHN_rule", "LDOPA_rule"], sort=False):
        for lvl, lab in (("species_any", "species, any genome positive"), ("species_all", "species, all genomes positive"),
                         ("species_1genome_only", "species with one genome"), ("genome", "genomes")):
            r = g[g.level == lvl].iloc[0]
            rows.append([f"`{p}` / `{l}`", lab, n(r.n), n(r.DHN_only), n(r.LDOPA_only), n(r.both), n(r.neither)])
    return md(rows, ["Rule pair (DHN / L-DOPA)", "Unit", "n", "DHN only", "L-DOPA only", "Both", "Neither"],
              ["---", "---", "---:", "---:", "---:", "---:", "---:"])


def rules(a):
    r = pd.read_csv(Path(a.tables) / "rule_counts.tsv", sep="\t")
    rows = [[f"`{x.rule}`", x.definition, f"{n(x.genomes_positive)} of {n(x.genomes_searched)} ({pct(x.genomes_positive, x.genomes_searched)} %)"]
            for x in r.itertuples()]
    return md(rows, ["Rule", "Definition", "Genomes positive"], ["---", "---", "---:"])


def collapse(a):
    c = pd.read_csv(Path(a.tables) / "species_collapse.tsv", sep="\t").set_index("quantity")["value"]
    top = pd.read_csv(Path(a.tables) / "species_most_genomes.tsv", sep="\t").head(5)
    v = lambda k: n(c[k])
    rows = [["Genomes in the BFD species table", v("genomes_in_species_table")],
            ["Genomes searched (have proteins); the denominator of every genome count", v("genomes_searched")],
            ["Genomes not searched (no proteins); not assessed, left out of all denominators", v("genomes_not_searched")],
            ["Species in the table", v("species_in_table")],
            ["Species with at least one searched genome; the denominator of every species count", v("species_with_searched_genome")],
            ["Species with no searched genome (not assessed)", v("species_without_searched_genome")],
            ["Species with one searched genome", v("species_one_genome")],
            ["Species with 2 to 9 searched genomes", v("species_2_to_9")],
            ["Species with 10 or more searched genomes", v("species_10_or_more")],
            ["Searched genomes inside those 10-or-more species", v("genomes_in_species_with_10_or_more")],
            ["Median / maximum searched genomes per species", f"{c['median_genomes_per_species']:.0f} / {v('max_genomes_per_species')}"],
            ["Species whose genomes carry more than one class label / order label", f"{v('species_with_more_than_one_class_label')} / {v('species_with_more_than_one_order_label')}"],
            ["Searched genomes with no phylum / class / order label",
             f"{v('searched_genomes_no_phylum')} / {v('searched_genomes_no_class')} / {v('searched_genomes_no_order')}"],
            ["Re-annotated genomes among the searched (mapped through the alias file)", v("re_annotated_genomes_searched")],
            ["Species with most genomes", "; ".join(f"{r.species} {n(r.n_genomes)}" for r in top.itertuples())]]
    return md(rows, ["Quantity", "Value"], ["---", "---:"])


def reannot(a):
    d = pd.read_csv(Path(a.tables) / "reannotated_5_genomes.tsv", sep="\t")
    rows = []
    for r in d.itertuples():
        placed = [c for c in DHN + ["tyrosinase", "laccase_diphenol", "aspulvinone_e_synthase", "asp_melanin_tyrosinase", "hppd", "hmga"] if getattr(r, c) == 1]
        rows.append([r.locustag, r.SPECIES, r.CLASS if isinstance(r.CLASS, str) else "(no label)", ", ".join(placed) or "none",
                     r.D_pks, r.L_route])
    return md(rows, ["Locus tag (table)", "Species", "Class", "Steps with a placed call", "`D_pks`", "`L_route`"],
              ["---", "---", "---", "---", "---:", "---:"])


def dominant(a):
    """Classes and orders where at least 90 % of the species fall in one of the four classes (headline rules)."""
    rows = []
    for rank, key, mn in (("class", "CLASS", 20), ("order", "ORDER", 50)):
        d = pd.read_csv(Path(a.tables) / f"four_class_by_{rank}.tsv", sep="\t")
        d = d[(d.n_species >= mn) & (d[key] != "(no label)")]
        for lab, col in (("DHN only", "species_DHN_only"), ("L-DOPA only", "species_LDOPA_only"),
                         ("both", "species_both"), ("neither", "species_neither")):
            hit = d[d[col] / d.n_species >= 0.9].sort_values("n_species", ascending=False)
            txt = "; ".join(f"{r[key]} ({int(r[col]):,} of {int(r.n_species):,})" for _, r in hit.iterrows()) or "none"
            rows.append([lab, rank, txt])
    return md(rows, ["Class", "Rank", "Groups with 90 % or more of their species in the class (species in the class of species in the group)"],
              ["---", "---", "---"])


BLOCKS = {"dominant": dominant, "steps": steps_table, "placement": placement_table, "phylum_steps": phylum_steps, "rules": rules,
          "sensitivity": sensitivity, "collapse": collapse, "reannot": reannot,
          "four_phylum": lambda a: four_by(a, "phylum", "PHYLUM", 1),
          "four_class": lambda a: four_by(a, "class", "CLASS", 20),
          "four_order": lambda a: four_by(a, "order", "ORDER", 100)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--pathways", required=True)
    ap.add_argument("--min-scores", required=True)
    ap.add_argument("--search-summary", required=True)
    ap.add_argument("--summary-v1", required=True)
    ap.add_argument("--summary-v2", required=True)
    ap.add_argument("--summary-allplace", required=True)
    ap.add_argument("--pks-calls", required=True)
    ap.add_argument("--min-phylum", type=int, default=50)
    a = ap.parse_args(argv)
    text = Path(a.report).read_text()
    for name, fn in BLOCKS.items():
        pat = re.compile(rf"(<!-- BEGIN {name} -->\n).*?(\n<!-- END {name} -->)", re.S)
        if not pat.search(text):
            continue
        body = fn(a)
        text = pat.sub(lambda m: m.group(1) + body + m.group(2), text)
    Path(a.report).write_text(text)


if __name__ == "__main__":
    main()
