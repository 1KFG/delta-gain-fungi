#!/usr/bin/env python3
"""Pathway diagram of one census pathway, drawn from its knowledge-base tables (PNG, PDF and SVG through Graphviz `dot`).

Reads kb/<pathway>/{compound,reaction,reaction_participant,step_reaction}.tsv and compound names from the ChEBI cache
(.cache/chebi/<id>.json, written by the knowledge-base fetchers). Nodes are compounds. An edge is a curated reaction; its
label is the step (gene family) that does it, and its style is the BEST evidence tier among the links of that reaction.
A reaction whose substrate or product is not in the compound table cannot be drawn as an edge. It is listed in a note box.

Reading trap, repeated in the report caption: the tier of an edge is the best link of that single reaction. The tier of a
whole route is the weakest step on it (`path_possible`), so an edge can be T1a while its route is T3.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shutil
import subprocess
from pathlib import Path

TIER_ORDER = ["T1a", "T1b", "T1c", "T2", "T3", "T4"]
# tier -> (colour, pen width, style): one hue, strongest evidence darkest and thickest; weakest dashed or dotted
TIER_STYLE = {
    "T1a": ("#08306b", 2.6, "solid"), "T1b": ("#2171b5", 2.2, "solid"), "T1c": ("#6baed6", 2.0, "solid"),
    "T2": ("#6baed6", 1.6, "dashed"), "T3": ("#6f6e6a", 1.5, "dashed"), "T4": ("#8a8985", 1.2, "dotted"),
}
# compounds that feed the chain but are not part of it; the first reaction gets one grey node instead
FEEDERS = {"CHEBI:16584", "CHEBI:17407"}   # isopentenyl diphosphate, farnesyl diphosphate
SIDE_PRODUCTS = {"CHEBI:15825"}            # 3-methylbut-2-enal, released by a cleavage


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def best_tier(tiers):
    ranked = [t for t in tiers if t in TIER_ORDER]
    return min(ranked, key=TIER_ORDER.index) if ranked else None


def short_name(name):
    """Drop stereo prefixes that make a node very wide, and give GGPP its usual short name."""
    name = re.sub(r"^\d+-trans,\d+-trans,\d+-trans-", "", name)
    return "geranylgeranyl diphosphate (GGPP)" if name == "geranylgeranyl diphosphate" else name


def compound_name(chebi_id, cache_dir):
    f = Path(cache_dir) / "chebi" / f"{chebi_id.split(':')[1]}.json"
    if f.exists():
        name = json.loads(f.read_text()).get("data", {}).get("name")
        if name:
            return short_name(html.unescape(re.sub(r"<[^>]+>", "", name)))
    return chebi_id


def build_graph(kb_dir):
    """Return (edges, undrawn). edges: list of dict(src, dst, steps, tier, reaction, ec). undrawn: list of (reaction, steps, why)."""
    kb = Path(kb_dir)
    rx = {r["reaction_id"]: r for r in read_tsv(kb / "reaction.tsv")}
    part = {}
    for r in read_tsv(kb / "reaction_participant.tsv"):
        part.setdefault(r["reaction_id"], []).append((r["compound"], r["role"]))
    links = {}
    for r in read_tsv(kb / "step_reaction.tsv"):
        links.setdefault(r["reaction_id"], []).append((r["step_id"], r["tier"]))
    edges, undrawn = [], []
    for rid in rx:
        steps = sorted({s for s, _ in links.get(rid, [])})
        tier = best_tier([t for _, t in links.get(rid, [])])
        subs = [c for c, role in part.get(rid, []) if role == "substrate"]
        prods = [c for c, role in part.get(rid, []) if role == "product" and c not in SIDE_PRODUCTS]
        if not subs or not prods:
            undrawn.append((rid, steps, "no substrate or no product compound in the tables"))
            continue
        if set(subs) <= FEEDERS:
            subs = ["FEED"]
        for s in subs:
            for p in prods:
                edges.append(dict(src=s, dst=p, steps=steps, tier=tier, reaction=rid, ec=rx[rid]["ec"]))
    return edges, undrawn


def dot_text(edges, undrawn, cache_dir, title, rankdir="TB"):
    nodes = {}
    for e in edges:
        for c in (e["src"], e["dst"]):
            nodes[c] = "IPP + FPP" if c == "FEED" else compound_name(c, cache_dir)
    out = ['digraph G {', f'  rankdir={rankdir}; bgcolor="#fcfcfb"; nodesep=0.3; ranksep=0.55; pad=0.3;',
           f'  labelloc=t; labeljust=l; fontname="Helvetica-Bold"; fontsize=15; label="{title}";',
           '  node [fontname="Helvetica", fontsize=11, shape=box, style="rounded,filled", fillcolor="#ffffff", color="#8a8985", penwidth=1.2];',
           '  edge [fontname="Helvetica", fontsize=9, arrowsize=0.8];']
    for c, name in nodes.items():
        nid = c.replace(":", "_")
        if c == "FEED":
            out.append(f'  {nid} [label="{name}", shape=ellipse, fillcolor="#efeeea", color="#b9b8b3", fontcolor="#52514e"];')
        else:
            out.append(f'  {nid} [label="{name}"];')
    for e in edges:
        colour, width, style = TIER_STYLE.get(e["tier"], TIER_STYLE["T4"])
        label = " / ".join(e["steps"]) + (f"\\n{e['tier']}" if e["tier"] else "")
        out.append(f'  {e["src"].replace(":", "_")} -> {e["dst"].replace(":", "_")} '
                   f'[label="{label}", color="{colour}", fontcolor="#0b0b0b", penwidth={width}, style={style}];')
    # legend
    rows = "".join(
        f'<TR><TD ALIGN="LEFT"><FONT COLOR="{TIER_STYLE[t][0]}"><B>{"&#9472;&#9472;" if TIER_STYLE[t][2]=="solid" else ("- - -" if TIER_STYLE[t][2]=="dashed" else "&#183; &#183; &#183;")}</B></FONT> {t}</TD></TR>'
        for t in ("T1a", "T1b", "T1c", "T3", "T4"))
    out.append(f'  legend [shape=plain, label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2"><TR><TD><B>Best evidence tier of the reaction</B></TD></TR>{rows}</TABLE>>];')
    if undrawn:
        lines = "<BR ALIGN=\"LEFT\"/>".join(f'{html.escape(", ".join(steps) or "-")}: {html.escape(r.split(":")[-1])}' for r, steps, _ in undrawn)
        out.append(f'  note [shape=plain, label=<<B>Steps with no drawable edge</B><BR ALIGN="LEFT"/>(a substrate or product is not in the tables)<BR ALIGN="LEFT"/>{lines}<BR ALIGN="LEFT"/>>, fontsize=9];')
    out.append("}")
    return "\n".join(out)


def render(dot_source, out_dir, name):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{name}.dot").write_text(dot_source)
    if shutil.which("dot") is None:
        raise RuntimeError("graphviz `dot` is not on PATH")
    for fmt, extra in (("png", ["-Gdpi=150"]), ("pdf", []), ("svg", [])):
        subprocess.run(["dot", f"-T{fmt}", *extra, str(out_dir / f"{name}.dot"), "-o", str(out_dir / f"{name}.{fmt}")], check=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kb-dir", required=True, type=Path)
    ap.add_argument("--cache-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--name", default="fig6_pathway_diagram")
    ap.add_argument("--rankdir", default="TB", choices=["TB", "LR"])
    ap.add_argument("--title", default="Fungal carotenoid biosynthesis as curated in the knowledge base")
    a = ap.parse_args(argv)
    edges, undrawn = build_graph(a.kb_dir)
    render(dot_text(edges, undrawn, a.cache_dir, a.title, a.rankdir), a.out_dir, a.name)
    print(f"drew {len(edges)} edges, {len(undrawn)} reactions listed without an edge")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
