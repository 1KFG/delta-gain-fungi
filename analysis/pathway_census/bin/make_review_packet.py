#!/usr/bin/env python3
"""Generate the stage 1 review packet (Markdown, Quarto format) for the project owner.

Usage: /usr/bin/python3.12 bin/make_review_packet.py [--out FILE]
Default output: <repo>/docs/superpowers/plans/2026-11-21-stage1-review-packet.md

Inputs (read only, no network): pathways/<name>.yaml, kb/<name>/*.tsv, refs/<name>/provenance.tsv,
refs/<name>/decoys/provenance.tsv, refs/<name>/UNVERIFIED.md, kb/VERIFICATION_REPORT.md and the
tier definitions of the design spec (section 6.2). The output is deterministic.

Sources of numbers. ``build_packet`` returns ``(text, numbers)``. ``numbers`` maps every number
string that the generator wrote to its source. A number is written by one of three helpers:
``Ctx.n`` (a value the generator computed), ``Ctx.c`` (text copied from a table, a YAML field,
the report or the spec) and ``Ctx.verbatim`` (a block copied from UNVERIFIED.md, not scanned).
Prose typed in this script contains no number except through ``Ctx.n``.

Decision items. A decision item is a block of UNVERIFIED.md (after its first ``##`` heading) that
contains ``USER DECISION`` or ``DECISION FOR THE USER``. A block is a list item (with its
continuation lines and a table that follows it), a paragraph, or a whole section when the heading
of the section contains a marker. Blocks are copied character for character.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG))

from pathway_census.kb import TIERS, formula_mass, load_kb, path_possible, species_key  # noqa: E402

PATHWAYS = ("carotenoid", "melanin_dhn", "melanin_ldopa", "mycosporine")
LETTER = {"carotenoid": "C", "melanin_dhn": "D", "melanin_ldopa": "L", "mycosporine": "M"}
MARKERS = ("USER DECISION", "DECISION FOR THE USER")
REVIEW_DATE = "2026-11-21"
SPEC_REL = Path("docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-design.md")
OUT_REL = Path("docs/superpowers/plans/2026-11-21-stage1-review-packet.md")
TOKEN_SPLIT = r"[\s|`(),;:\[\]*_/=<>]+"
VERBATIM_RE = r"^::: \{\.verbatim\}\n(.*?)\n:::$"
ITEM_RE = re.compile(r"(\d+\.|[-*]) ")
HEAD_RE = re.compile(r"(#+) ")


def strip_verbatim(text: str) -> str:
    return re.sub(VERBATIM_RE, "", text, flags=re.S | re.M)


def number_tokens(text: str) -> list[str]:
    out = []
    for tok in re.split(TOKEN_SPLIT, text):
        tok = tok.strip(".'\"")
        if re.fullmatch(r"\d+(\.\d+)?", tok):
            out.append(tok)
    return out


class Ctx:
    """Collects the numbers that the generator writes, with their source."""

    def __init__(self) -> None:
        self.numbers: dict[str, str] = {}

    def n(self, value, source: str) -> str:
        s = str(value)
        self.numbers.setdefault(s, source)
        return s

    def c(self, text, source: str) -> str:
        s = str(text)
        for tok in number_tokens(s):
            self.numbers.setdefault(tok, source)
        return s

    def cell(self, text, source: str) -> str:
        return self.c(str(text).replace("\\|", "\x00").replace("|", "\\|").replace("\x00", "\\|")
                      .replace("\n", " "), source)

    def verbatim(self, text: str) -> str:
        return "::: {.verbatim}\n" + text + "\n:::"


def table(header: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


# ----------------------------------------------------------------------------- UNVERIFIED.md

def decision_blocks(text: str) -> list[dict]:
    """Decision blocks of an UNVERIFIED.md, in file order. Each has text, first, last (1-based)."""
    lines = text.split("\n")
    start = next(i for i, ln in enumerate(lines) if ln.startswith("## "))
    out = []
    i = start
    while i < len(lines):
        ln = lines[i]
        m = HEAD_RE.match(ln)
        if m and any(k in ln for k in MARKERS):
            level = len(m.group(1))
            j = i + 1
            while j < len(lines):
                m2 = HEAD_RE.match(lines[j])
                if m2 and len(m2.group(1)) <= level:
                    break
                j += 1
            out.append(_block(lines, i, j, "section"))
            i = j
            continue
        if not ln.strip() or m or ln.strip() == "---":
            i += 1
            continue
        is_item = bool(ITEM_RE.match(ln))
        j = i + 1
        while j < len(lines):
            nl = lines[j]
            if HEAD_RE.match(nl) or nl.strip() == "---":
                break
            if ITEM_RE.match(nl):
                break
            if not nl.strip():
                k = j
                while k < len(lines) and not lines[k].strip():
                    k += 1
                if (is_item and k < len(lines)
                        and (lines[k].startswith("|") or lines[k].startswith((" ", "\t")))):
                    j = k
                    continue
                break
            j += 1
        blk = _block(lines, i, j, "item" if is_item else "paragraph")
        if any(k in blk["text"] for k in MARKERS):
            out.append(blk)
        i = j
    return out


def _block(lines: list[str], i: int, j: int, kind: str) -> dict:
    seg = lines[i:j]
    while seg and not seg[-1].strip():
        seg.pop()
    return {"text": "\n".join(seg), "first": i + 1, "last": i + len(seg), "kind": kind}


def section_blocks(text: str, pattern: str) -> list[dict]:
    """Sections whose heading matches ``pattern``, copied whole (heading to the next heading of
    the same or a higher level)."""
    lines = text.split("\n")
    out = []
    i = 0
    while i < len(lines):
        m = HEAD_RE.match(lines[i])
        if m and re.search(pattern, lines[i]):
            level = len(m.group(1))
            j = i + 1
            while j < len(lines):
                m2 = HEAD_RE.match(lines[j])
                if m2 and len(m2.group(1)) <= level:
                    break
                j += 1
            out.append(_block(lines, i, j, "section"))
            i = j
        else:
            i += 1
    return out


def owner_sections(text: str) -> list[dict]:
    """Subsections (###) of the hand-written 'Items for a person to check' block, in file order,
    except those whose heading holds a decision marker (those are numbered decisions)."""
    lines = text.split("\n")
    start = next(i for i, ln in enumerate(lines) if ln.startswith("## Items for a person to check"))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    out = []
    i = start + 1
    while i < end:
        if lines[i].startswith("### "):
            j = i + 1
            while j < end and not (HEAD_RE.match(lines[j]) and len(HEAD_RE.match(lines[j]).group(1)) <= 3):
                j += 1
            if not any(k in lines[i] for k in MARKERS):
                out.append(_block(lines, i, j, "section"))
            i = j
        else:
            i += 1
    return out


def context_numbers(text: str) -> str:
    """The route-tier lines of the decision section that give the context reading. Copied
    verbatim: the table after 'Route tiers, MEASURED' (carotenoid) or the numbered item that
    starts with it (melanin_dhn)."""
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        if ln.startswith("Route tiers, MEASURED"):
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            k = j
            while k < len(lines) and lines[k].startswith("|"):
                k += 1
            return "\n".join(lines[i:k])
        if re.match(r"\d+\. \*\*Route tiers, MEASURED", ln):
            return ln
    raise ValueError("no route tier statement found in the decision section")


def md_tables(block_text: str) -> list[tuple[list[str], list[list[str]]]]:
    out = []
    cur: list[str] = []
    for ln in block_text.split("\n") + [""]:
        if ln.startswith("|"):
            cur.append(ln)
        elif cur:
            cells = [[c.strip() for c in re.split(r"(?<!\\)\|", r.strip())[1:-1]] for r in cur]
            if len(cells) >= 2:
                out.append((cells[0], cells[2:]))
            cur = []
    return out


# ----------------------------------------------------------------------------- report

def parse_report(text: str) -> dict:
    def sec(start: str, end: str | None) -> list[str]:
        body = text.split(start, 1)[1]
        if end:
            body = body.split(end, 1)[0]
        return body.split("\n")

    def rows(lines: list[str]) -> list[str]:
        return [ln for ln in lines if ln.startswith("| ") and not ln.startswith("| ---")
                and not ln.startswith("|---")][1:]  # drop the header row

    def split(raw: str) -> list[str]:
        return [c.strip() for c in re.split(r"(?<!\\)\|", raw)[1:-1]]

    per_path = rows(sec("## Counts per pathway", "## Counts per check"))
    per_check = rows(sec("## Counts per check", "## Gene-name passes by kind"))
    kinds = rows(sec("## Gene-name passes by kind", "## Failures"))
    fails = rows(sec("## Failures", "## Unchecked"))
    unchecked = rows(sec("## Unchecked", None))
    return {"per_path": per_path, "per_check": per_check, "kinds": kinds, "fails": fails,
            "unchecked": unchecked, "split": split}


def normalise_reason(reason: str) -> str:
    r = reason.split(":")[0]
    r = re.sub(r"\s*\(flag for review, tier T\d\w*\)", "", r)
    return r.strip()


# ----------------------------------------------------------------------------- per pathway

def load_pathway_data(root: Path, name: str) -> dict:
    y = yaml.safe_load((root / "pathways" / f"{name}.yaml").read_text(encoding="utf-8"))
    kb = load_kb(root / "kb" / name)
    unv = (root / "refs" / name / "UNVERIFIED.md").read_text(encoding="utf-8")
    dprov = pd.read_csv(root / "refs" / name / "decoys" / "provenance.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    return {"y": y, "kb": kb, "unv": unv, "dprov": dprov}


def compact(v) -> str:
    return "-" if v is None or str(v) == "" else str(v)


def route_rows(ctx: Ctx, name: str, d: dict) -> tuple[list[list[str]], list[list[str]]]:
    """Rows of the routes table, and the routes with no product (id, why)."""
    kb = d["kb"]
    formulas = dict(zip(kb["compound"]["chebi_id"], kb["compound"]["formula"]))
    rows, none = [], []
    for r in d["y"]["routes"]:
        try:
            df = path_possible(kb, r["step_ids"])
        except ValueError as e:
            none.append([f"`{r['id']}`", ctx.cell(f"path_possible cannot run: {e}", "path_possible error")])
            continue
        if not len(df):
            none.append([f"`{r['id']}`", "path_possible returns no compound"])
            continue
        for p in df.itertuples():
            rows.append([f"`{r['id']}`", p.compound, ctx.cell(formulas.get(p.compound, "-"), "compound.tsv"),
                         p.min_tier, f"`{p.organisms}`", str(bool(p.single_organism)),
                         ctx.n(p.n_equal_paths, "path_possible n_equal_paths")])
    return rows, none


ROUTE_HEAD = ["route", "product", "formula", "min tier", "organisms", "single_organism", "n_equal_paths"]


def steps_section(ctx: Ctx, name: str, d: dict) -> list[str]:
    sr = d["kb"]["step_reaction"]
    out = ["### Steps and links", "",
           "One row per link (step, reaction, protein). Links at tier T4 are in the T4 section.", ""]
    rows = []
    for st in d["y"]["steps"]:
        links = sr[(sr["step_id"] == st["id"]) & (sr["tier"] != "T4")]
        enzyme = ctx.cell(st["name"], f"{name}.yaml step name")
        if not len(links):
            rows.append([compact(st.get("route")), f"`{st['id']}`", enzyme, "-", "-", "-", "-", "-", "-"])
        for lk in links.itertuples():
            pm = f"PMID {lk.pmid}" if str(lk.pmid) not in ("", "nan") else "-"
            rows.append([compact(st.get("route")), f"`{st['id']}`", enzyme,
                         ctx.cell(species_key(lk.organism), "step_reaction.tsv"),
                         ctx.cell(compact(lk.gene_name), "step_reaction.tsv"),
                         ctx.cell(compact(lk.protein_accession), "step_reaction.tsv"),
                         lk.tier, lk.reaction_id.replace("|", "\\|"), pm])
    out += table(["route", "step", "enzyme", "organism", "gene", "accession", "tier", "reaction", "PMID"],
                 [[ctx.c(c, "step_reaction.tsv") for c in r] for r in rows])
    return out + [""]


def decoys_section(ctx: Ctx, name: str, d: dict) -> list[str]:
    dp = d["dprov"]
    rows = []
    for st in d["y"]["steps"]:
        for dc in st.get("decoys") or []:
            hit = dp[(dp["step_id"] == st["id"]) & (dp["accession"] == dc["accession"])]
            reason = hit.iloc[0]["reason"] if len(hit) else "-"
            rows.append([f"`{st['id']}`", ctx.cell(dc["accession"], "yaml decoy"),
                         ctx.cell(species_key(dc["species"]), "yaml decoy"),
                         ctx.cell(compact(dc.get("gene_name")), "yaml decoy"),
                         ctx.cell(reason, "decoys/provenance.tsv")])
    out = ["### Decoys", "",
           "A decoy is a protein that shares the family of the step and has another function. "
           "The list is PROPOSED by the curator.", ""]
    if rows:
        out += table(["step", "accession", "species", "gene", "reason"], rows)
    else:
        out.append("None.")
    return out + [""]


def compounds_section(ctx: Ctx, name: str, d: dict, blocks: list[dict]) -> list[str]:
    comp = d["kb"]["compound"]
    rows, differ = [], 0
    for r in comp.itertuples():
        calc = f"{formula_mass(r.formula):.4f}"
        if calc != f"{float(r.neutral_exact_mass):.4f}":
            differ += 1
        rows.append([r.chebi_id, ctx.cell(r.formula, "compound.tsv"), ctx.n(calc, "kb.formula_mass"),
                     ctx.cell(r.neutral_exact_mass, "compound.tsv"), ctx.cell(compact(r.kegg_id), "compound.tsv"),
                     ctx.cell(r.inchikey, "compound.tsv")])
    out = ["### Compounds", "",
           "Mass is the neutral exact mass from `kb.formula_mass`. The table mass is the stored value. "
           "ChEBI is the key. KEGG and the InChIKey are cross-reference columns. "
           "PubChem is looked up by InChIKey (the table stores no PubChem id).", ""]
    out += table(["ChEBI", "formula", "mass (formula_mass)", "mass (table)", "KEGG", "InChIKey"], rows)
    out += ["", f"Compounds: {ctx.n(len(comp), 'compound.tsv rows')}. "
            f"Rows where the two masses differ at four decimals: {ctx.n(differ, 'formula_mass vs table')}.", ""]
    # no ChEBI id
    out += ["### No ChEBI id (owner decision)", ""]
    part = d["kb"]["reaction_participant"]
    prod = set(part[part["role"] == "product"]["reaction_id"])
    rx = d["kb"]["reaction"]
    noprod = rx[~rx["reaction_id"].isin(prod)]
    found = False
    for _hdr, trows in [(h, r) for blk in blocks for h, r in md_tables(blk["text"]) if h and h[0] == "Compound"]:
        found = True
        out += ["The decision item lists these compounds (copied from UNVERIFIED.md). "
                "None has a row in `compound.tsv`.", ""]
        out += table(["Compound", "Formula (PubChem)", "PubChem CID", "Effect"],
                     [[ctx.cell(c, "UNVERIFIED.md table") for c in r] for r in trows])
        out.append("")
    if not found:
        out += ["No table of compounds without a ChEBI id is in the decision items of this pathway. "
                "The compounds are named in the prose of the decision items and of the reaction notes.", ""]
    out += [f"Reactions with no product participant in the tables: {ctx.n(len(noprod), 'reaction_participant.tsv')}.", ""]
    if len(noprod):
        out += table(["reaction", "curated", "sentence of the note that names ChEBI"],
                     [[r.reaction_id, r.curated,
                       ctx.cell(next((s for s in re.split(r"(?<=\.) ", r.note) if "ChEBI" in s), "-"),
                                "reaction.tsv note")] for r in noprod.itertuples()])
        out.append("")
    return out


def t4_section(ctx: Ctx, name: str, d: dict) -> list[str]:
    sr = d["kb"]["step_reaction"]
    t4 = sr[sr["tier"] == "T4"]
    out = ["### T4 and excluded items", "",
           "A path through a T4 link is excluded from `path_possible`.", ""]
    out.append(f"T4 links in `step_reaction.tsv`: {ctx.n(len(t4), 'step_reaction.tsv tier T4')}.")
    out.append("")
    if len(t4):
        out += table(["step", "reaction", "organism", "PMID"],
                     [[f"`{r.step_id}`", ctx.c(r.reaction_id, "step_reaction.tsv"),
                       ctx.cell(species_key(r.organism), "step_reaction.tsv"),
                       f"PMID {r.pmid}" if str(r.pmid) not in ("", "nan") else "-"] for r in t4.itertuples()])
        out.append("")
    ystep = {s["id"] for s in d["y"]["steps"]}
    only_kb = [s for s in d["kb"]["step"]["step_id"] if s not in ystep]
    out.append("Steps in `step.tsv` and not in the YAML: " +
               (", ".join(f"`{s}`" for s in only_kb) if only_kb else "none") + ".")
    noref = [s["id"] for s in d["y"]["steps"] if not s.get("refs")]
    out.append("")
    out.append("Steps in the YAML with no reference protein: " +
               (", ".join(f"`{s}`" for s in noref) if noref else "none") + ".")
    notes = [(r["id"], r.get("note") or "") for r in d["y"]["routes"]
             if re.search(r"T4|excluded", r.get("note") or "")]
    out.append("")
    if notes:
        out.append("Route notes that mention T4 or exclusion:")
        out.append("")
        out += [f"- `{rid}`: {ctx.c(note, f'{name}.yaml route note')}" for rid, note in notes]
    else:
        out.append("No route note mentions T4 or exclusion.")
    secs = section_blocks(d["unv"], r"T4|[Ee]xcluded")
    out.append("")
    if secs:
        out += ["The section of `UNVERIFIED.md` on these items, copied verbatim:", ""]
        for s in secs:
            out += [ctx.verbatim(s["text"]), ""]
    else:
        out += ["`UNVERIFIED.md` has no section on T4 or excluded items.", ""]
    return out + [""]


def unverified_refs_section(ctx: Ctx, name: str, d: dict) -> list[str]:
    out = ["### Unverified references (`verified: false`)", "",
           "`verified: true` means database checks passed on the review date. It never means a person reviewed the reference.", ""]
    rows = [(s["id"], "reference", r) for s in d["y"]["steps"] for r in s.get("refs") or [] if not r.get("verified")]
    rows += [(s["id"], "decoy", r) for s in d["y"]["steps"] for r in s.get("decoys") or [] if not r.get("verified")]
    out.append(f"References and decoys with `verified: false`: {ctx.n(len(rows), f'{name}.yaml')}.")
    out.append("")
    if not rows:
        return out + [""]
    out += table(["step", "kind", "accession", "species", "gene", "PMID"],
                 [[f"`{sid}`", kind, ctx.cell(compact(r.get("accession")), "yaml"),
                   ctx.cell(species_key(r["species"]), "yaml"), ctx.cell(compact(r.get("gene_name")), "yaml"),
                   f"PMID {r['pmid']}" if r.get("pmid") else "-"] for sid, kind, r in rows])
    out += ["", "Reason, as written by the curator in `UNVERIFIED.md` (verbatim):", ""]
    unv = d["unv"]
    top = unv.split("\n## Refs with verified: false", 1)[1].split("\n## ", 1)[0] if "\n## Refs with verified: false" in unv else ""
    paras = [p for p in re.split(r"\n\s*\n", top) if p.strip() and not p.lstrip().startswith("|")]
    for sid, kind, r in rows:
        acc = r.get("accession")
        hit = [p.strip() for p in paras if acc and acc in p]
        if not hit:
            hit = [p.strip() for p in re.split(r"\n\s*\n", unv) if acc and acc in p and not p.lstrip().startswith("|")][:1]
        out.append(f"**`{sid}`, {acc}**")
        out.append("")
        out += [ctx.verbatim(hit[0]) if hit else "No note for this accession in `UNVERIFIED.md`.", ""]
    return out


def link_counts_section(ctx: Ctx, name: str, d: dict) -> list[str]:
    sr = d["kb"]["step_reaction"]
    vc = sr["tier"].value_counts()
    rows = [[name, t, ctx.n(int(vc[t]), "step_reaction.tsv")] for t in TIERS if t in vc.index]
    return ["### Links per tier", ""] + table(["pathway", "tier", "links"], rows) + [""]


def owner_items_section(ctx: Ctx, name: str, secs: list[dict]) -> list[str]:
    letters = ", ".join(s["text"].split("\n", 1)[0][4:].split(".")[0] for s in secs)
    out = [f"### Items for the owner ({name}, hand-written sections {letters}, not numbered decisions)", "",
           f"Copied character for character from `refs/{name}/UNVERIFIED.md`, block 'Items for a person to check'. "
           "The sections with a decision marker in their heading are in the numbered decisions below. "
           "An item that holds a marker appears in both places.", ""]
    for s in secs:
        out += [ctx.verbatim(s["text"]), ""]
    return out


def decisions_section(ctx: Ctx, name: str, blocks: list[dict]) -> list[str]:
    L = LETTER[name]
    out = ["### Decisions for the owner", "",
           f"Each item is copied character for character from `refs/{name}/UNVERIFIED.md`. "
           f"The number is `D{L}` plus the position in the file.", ""]
    for i, b in enumerate(blocks, 1):
        out += [f"**D{L}{i}** (`refs/{name}/UNVERIFIED.md`, lines {ctx.n(b['first'], 'UNVERIFIED.md line')} "
                f"to {ctx.n(b['last'], 'UNVERIFIED.md line')}, {b['kind']})", "", ctx.verbatim(b["text"]), ""]
    return out


def verification_section(ctx: Ctx, name: str, rep: dict) -> list[str]:
    split = rep["split"]
    out = ["### Verification", ""]
    row = next(r for r in rep["per_path"] if split(r)[0] == name)
    out += table(["pathway", "checks", "pass", "fail", "unchecked"],
                 [[ctx.c(c, "VERIFICATION_REPORT.md") for c in split(row)]])
    out.append("")
    fails = [r for r in rep["fails"] if split(r)[0] == name]
    out.append(f"Failures in this pathway: {ctx.n(len(fails), 'VERIFICATION_REPORT.md failures')}. "
               "Every failure row is in the verification section at the end of the packet.")
    unc = [split(r) for r in rep["unchecked"] if split(r)[0] == name]
    counts: dict[str, int] = {}
    for c in unc:
        counts[c[1]] = counts.get(c[1], 0) + 1
    out.append("")
    if counts:
        out += table(["pathway", "unchecked check", "rows"],
                     [[name, k, ctx.n(v, "VERIFICATION_REPORT.md unchecked")] for k, v in sorted(counts.items())])
    else:
        out.append("No unchecked rows.")
    out.append("")
    out.append("The report gives pass counts per check for all pathways together only. "
               "It gives no pass count per check per pathway. This packet does not recompute them.")
    return out + [""]


# ----------------------------------------------------------------------------- packet

def build_packet(root: Path = PKG) -> tuple[str, dict]:
    root = Path(root)
    repo = root.parents[1]
    ctx = Ctx()
    data = {n: load_pathway_data(root, n) for n in PATHWAYS}
    blocks = {n: decision_blocks(data[n]["unv"]) for n in PATHWAYS}
    rep = parse_report((root / "kb" / "VERIFICATION_REPORT.md").read_text(encoding="utf-8"))
    spec = (repo / SPEC_REL).read_text(encoding="utf-8")
    tier_lines = [ln for ln in spec.split("### 6.2 ", 1)[1].split("### 6.3", 1)[0].split("\n")[1:] if ln.strip()]

    out: list[str] = ["---", 'title: "Stage 1 review packet"',
                      'subtitle: "Verified step sets, references, decoys and flat tables: four pathways"',
                      f"date: {REVIEW_DATE}", "toc: true", "toc-depth: 3", "---", ""]
    out += ["::: {.callout-important}", "## Limits of this packet", "",
            "- This file is generated by `bin/make_review_packet.py` from the committed YAML files, tables, "
            "`UNVERIFIED.md` files and the verification report. Do not edit it by hand.",
            "- Nothing here is a biological finding.",
            "- No real-data run has been done. No genome was searched.",
            "- A tier is the curator's reading of an abstract. A verification pass means a public record "
            "agrees with a file. It does not mean the biology is right.",
            "- This packet is the review gate. The owner decides before any real-data run.", ":::", ""]

    # ---- rule 1 decision at the top
    out += ["## Decision on rule 1: strict or context reading", "",
            "Rule 1 sets when a link is T1a. The strict reading needs the quoted sentence to name the protein "
            "of the link. The context reading accepts a protein that is named earlier in the abstract. "
            "The tables use the strict reading. The owner decides between the two readings.", ""]
    for n in ("carotenoid", "melanin_dhn"):
        ids = [f"D{LETTER[n]}{i}" for i, b in enumerate(blocks[n], 1)
               if b["kind"] == "section" and "strict or context reading" in b["text"].split("\n", 1)[0]]
        out.append(f"- `{n}`: the full text is in decision {', '.join(ids)}.")
    out += ["- `melanin_ldopa` and `mycosporine` have no section on this reading. The `melanin_dhn` text says "
            "that one decision of the owner should cover the three pathways `carotenoid`, `melanin_dhn` and "
            "`melanin_ldopa`.", ""]
    for n in ("carotenoid", "melanin_dhn"):
        sec = next(b for b in blocks[n] if b["kind"] == "section" and "strict or context reading" in b["text"].split("\n", 1)[0])
        rows, _ = route_rows(ctx, n, data[n])
        out += [f"### Route tiers under each reading: `{n}`", "",
                "Strict reading. Computed with `path_possible` on the committed tables.", ""]
        out += table(ROUTE_HEAD, rows)
        out += ["", "Context reading. As stated in UNVERIFIED.md, not recomputed. Copied verbatim:", "",
                ctx.verbatim(context_numbers(sec["text"])), ""]

    # ---- how to read
    labels = next(ln for ln in data["carotenoid"]["unv"].split("\n") if ln.startswith("Labels:"))
    out += ["## How to read this packet", "",
            "This packet is for review. It has one section per pathway and a verification section at the end.", "",
            f"Tier definitions. Copied from the design spec (section {ctx.n('6.2', 'spec section label')}), "
            f"`{SPEC_REL}`:", ""]
    out += [ctx.c(ln, "spec 6.2") for ln in tier_lines]
    out += ["", "Labels used in the curator's notes (copied from `UNVERIFIED.md`):", "",
            ctx.c(labels, "carotenoid UNVERIFIED.md Labels line"), "",
            "Reading the route tables:", "",
            "- `min tier` is the weakest tier along the best path.",
            "- The links of that path can belong to different organisms. `organisms` lists them.",
            "- When `single_organism` is False, the tier does not hold for any one organism. "
            "Read it only together with the organism list.",
            "- `n_equal_paths` counts organism lists with the best tier and the fewest organisms. "
            "It does not count routes through the reactions.",
            "- A route with no product row cannot be followed by the tables. The reason is in the "
            "decision items of the pathway.", "",
            "Decision numbers: `DC` carotenoid, `DD` melanin_dhn, `DL` melanin_ldopa, `DM` mycosporine. "
            "The letter is followed by the position of the item in `UNVERIFIED.md`.", ""]

    # ---- overview
    ov = []
    for n in PATHWAYS:
        d = data[n]
        sr = d["kb"]["step_reaction"]
        nunver = sum(1 for s in d["y"]["steps"] for r in (s.get("refs") or []) if not r.get("verified"))
        ndec = sum(len(s.get("decoys") or []) for s in d["y"]["steps"])
        ov.append([n, ctx.n(len(d["y"]["steps"]), "yaml steps"), ctx.n(len(d["y"]["routes"]), "yaml routes"),
                   ctx.n(len(sr), "step_reaction.tsv"), ctx.n(ndec, "yaml decoys"),
                   ctx.n(len(d["kb"]["compound"]), "compound.tsv"), ctx.n(nunver, "yaml refs verified false"),
                   ctx.n(len(blocks[n]), "decision blocks")])
    out += ["## Overview", ""]
    out += table(["pathway", "steps", "routes", "links", "decoys", "compounds", "unverified refs", "decision items"], ov)
    out += [""]

    # ---- per pathway
    for n in PATHWAYS:
        d = data[n]
        out += [f"## Pathway `{n}`", ""]
        out += owner_items_section(ctx, n, owner_sections(d["unv"]))
        out += decisions_section(ctx, n, blocks[n])
        rows, none = route_rows(ctx, n, d)
        out += ["### Routes", "", "From `path_possible` on the committed tables.", ""]
        out += table(ROUTE_HEAD, rows) + [""]
        if none:
            out += ["Routes with no product row:", ""]
            out += table(["route", "reason"], none) + [""]
        out += steps_section(ctx, n, d)
        out += decoys_section(ctx, n, d)
        out += compounds_section(ctx, n, d, blocks[n])
        out += t4_section(ctx, n, d)
        out += unverified_refs_section(ctx, n, d)
        out += link_counts_section(ctx, n, d)
        out += verification_section(ctx, n, rep)

    # ---- verification, all pathways
    split = rep["split"]
    out += ["## Verification, all pathways", "",
            "Copied from `kb/VERIFICATION_REPORT.md`. No check was re-run for this packet. "
            "The report states its own limits. A pass is not a finding.", "",
            "### Counts per pathway", ""]
    out += ["| pathway | checks | pass | fail | unchecked |", "|---|---|---|---|---|"]
    out += [ctx.c(r, "VERIFICATION_REPORT.md") for r in rep["per_path"]]
    out += ["", "### Counts per check", "", "| check | pass | fail | unchecked |", "|---|---|---|---|"]
    out += [ctx.c(r, "VERIFICATION_REPORT.md") for r in rep["per_check"]]
    out += ["", "### Gene-name passes by kind", "",
            "| pathway | gene name | gene name (synonym) | ORF name | ordered locus name |", "|---|---|---|---|---|"]
    out += [ctx.c(r, "VERIFICATION_REPORT.md") for r in rep["kinds"]]
    out += ["", f"### Failures ({ctx.n(len(rep['fails']), 'VERIFICATION_REPORT.md failures')})", "",
            "| pathway | check | subject | expected | found |", "|---|---|---|---|---|"]
    out += [ctx.c(r, "VERIFICATION_REPORT.md") for r in rep["fails"]]
    groups: dict[tuple, int] = {}
    for r in rep["unchecked"]:
        c = split(r)
        key = (c[0], c[1], normalise_reason(c[3]))
        groups[key] = groups.get(key, 0) + 1
    out += ["", "### Unchecked, by check and reason", "",
            f"Unchecked rows: {ctx.n(len(rep['unchecked']), 'VERIFICATION_REPORT.md unchecked')}. "
            "Each row below groups rows with the same check and the same reason up to the first colon.", "",
            "| pathway | check | reason | rows |", "|---|---|---|---|"]
    out += [f"| {p} | {ck} | {ctx.cell(rs, 'VERIFICATION_REPORT.md')} | {ctx.n(v, 'VERIFICATION_REPORT.md unchecked')} |"
            for (p, ck, rs), v in sorted(groups.items())]
    out.append("")
    return "\n".join(out), ctx.numbers


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=PKG.parents[1] / OUT_REL)
    a = ap.parse_args(argv)
    text, numbers = build_packet(PKG)
    a.out.write_text(text, encoding="utf-8")
    print(f"wrote {a.out}: {len(text.splitlines())} lines, {len(numbers)} distinct numbers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
