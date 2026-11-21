"""Flat knowledge tables for pathway curation (spec section 6.1) and their validator.

One directory per pathway under ``kb/<pathway>/`` holds five TSV files (header required, UTF-8):
compound, reaction, reaction_participant, step, step_reaction. ``taxon_capability`` is added by the
census later and is not part of this module.

``validate_kb`` returns error strings. Each names the table, the row key and the rule. The exact
formats are the ERR_* templates below (tests assert them).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

TIERS = ("T1a", "T1b", "T1c", "T2", "T3", "T4")
EVIDENCE_TIERS = ("T1a", "T1b", "T1c", "T2")  # these need pmid and quote
ROLES = ("substrate", "product")
MASS_TOL = 0.0005

COLUMNS = {
    "compound": ["chebi_id", "inchikey", "inchikey_block1", "smiles", "formula",
                 "neutral_exact_mass", "kegg_id", "np_class"],
    "reaction": ["reaction_id", "ec", "kegg_rn", "direction", "curated", "note", "rhea_id"],
    "reaction_participant": ["reaction_id", "compound", "role", "stoich"],
    "step": ["step_id", "route", "ko", "hmm_id", "family", "optional"],
    "step_reaction": ["step_id", "reaction_id", "organism", "gene_name", "protein_accession",
                      "tier", "assay", "pmid", "quote"],
}

# Error templates (single definition; tests assert these exact strings).
ERR_MISSING_COL = "{table}: missing column '{col}'"
ERR_EXTRA_COL = "{table}: unexpected column '{col}' (kegg_id is informational only)"
ERR_TIER = "step_reaction[{step},{rxn}]: tier '{tier}' not in " + ",".join(TIERS)
ERR_EMPTY = "step_reaction[{step},{rxn}]: {field} is empty"
ERR_EVIDENCE = "step_reaction[{step},{rxn}]: tier '{tier}' requires {field}"
ERR_PART_COMPOUND = "reaction_participant[{rxn},{cpd}]: compound not in compound table"
ERR_ROLE = "reaction_participant[{rxn},{cpd}]: role '{role}' not in " + ",".join(ROLES)
ERR_BLOCK1 = "compound[{cpd}]: inchikey_block1 '{got}' != first 14 characters of inchikey '{want}'"
ERR_CHARGED = "compound[{cpd}]: charged species not allowed, use the neutral parent"
ERR_MASS = ("compound[{cpd}]: neutral_exact_mass {got} differs from formula mass {want:.4f} "
            "by more than 0.0005")
ERR_MASS_NUM = "compound[{cpd}]: neutral_exact_mass '{got}' is not a number"
ERR_FORMULA = "compound[{cpd}]: formula '{formula}' cannot be parsed ({why})"
ERR_CURATED = "reaction[{rxn}]: curated=false requires rhea_id"

# Monoisotopic masses (mass of the most abundant isotope, unified atomic mass units).
# Source: NIST Atomic Weights and Isotopic Compositions for All Elements
# (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl), isotopes 1H, 12C, 14N, 16O,
# 23Na, 31P, 32S, 35Cl, 39K, 56Fe.
MONO_MASS = {
    "H": 1.00782503207,
    "C": 12.0,
    "N": 14.0030740048,
    "O": 15.99491461956,
    "Na": 22.9897692809,
    "P": 30.97376163,
    "S": 31.97207100,
    "Cl": 34.96885268,
    "K": 38.96370668,
    "Fe": 55.9349375,
}

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_SMILES_CHARGE = re.compile(r"\[[^\]]*[+-][^\]]*\]")


def formula_mass(formula: str) -> float:
    """Monoisotopic mass of a neutral formula such as ``C40H52O4``.

    Supports C, H, N, O, P, S, Na, K, Cl, Fe. Raises ValueError on any other element or on text
    that is not a plain formula (no charges, no parentheses, no hydrates).
    """
    formula = formula.strip()
    if not formula:
        raise ValueError("empty formula")
    pos, total = 0, 0.0
    while pos < len(formula):
        m = _TOKEN.match(formula, pos)
        if not m:
            raise ValueError(f"cannot parse formula {formula!r} at position {pos}")
        el, n = m.group(1), m.group(2)
        if el not in MONO_MASS:
            raise ValueError(f"unsupported element {el!r}")
        total += MONO_MASS[el] * (int(n) if n else 1)
        pos = m.end()
    return total


def load_kb(kb_dir: Path) -> dict[str, pd.DataFrame]:
    """Read the five TSV tables of one pathway. All cells are strings; empty cells are ''."""
    kb_dir = Path(kb_dir)
    out = {}
    for name in COLUMNS:
        out[name] = pd.read_csv(kb_dir / f"{name}.tsv", sep="\t", dtype=str,
                                keep_default_na=False, encoding="utf-8")
    return out


def _blank(v) -> bool:
    return v is None or (isinstance(v, str) and not v.strip()) or (not isinstance(v, str) and pd.isna(v))


def validate_kb(tables: dict[str, pd.DataFrame]) -> list[str]:
    """Return a list of error strings (empty when valid). See the ERR_* templates."""
    errs: list[str] = []
    for name, cols in COLUMNS.items():
        if name not in tables:
            errs.append(f"{name}: table missing")
            continue
        have = list(tables[name].columns)
        for c in cols:
            if c not in have:
                errs.append(ERR_MISSING_COL.format(table=name, col=c))
        for c in have:
            if c not in cols:
                errs.append(ERR_EXTRA_COL.format(table=name, col=c))
    if errs:
        return errs

    comp, rxn = tables["compound"], tables["reaction"]
    part, sr = tables["reaction_participant"], tables["step_reaction"]

    for _, r in comp.iterrows():
        cpd = r["chebi_id"]
        ik, b1 = str(r["inchikey"]), str(r["inchikey_block1"])
        if b1 != ik[:14]:
            errs.append(ERR_BLOCK1.format(cpd=cpd, got=b1, want=ik[:14]))
        formula = str(r["formula"])
        charged = bool(re.search(r"[+-]", formula)) or bool(_SMILES_CHARGE.search(str(r["smiles"])))
        if charged:
            errs.append(ERR_CHARGED.format(cpd=cpd))
            continue
        try:
            want = formula_mass(formula)
        except ValueError as e:
            why = str(e)
            errs.append(ERR_FORMULA.format(cpd=cpd, formula=formula, why=why))
            continue
        try:
            got = float(r["neutral_exact_mass"])
        except (TypeError, ValueError):
            errs.append(ERR_MASS_NUM.format(cpd=cpd, got=r["neutral_exact_mass"]))
            continue
        if abs(got - want) > MASS_TOL:
            errs.append(ERR_MASS.format(cpd=cpd, got=r["neutral_exact_mass"], want=want))

    known = set(comp["chebi_id"])
    for _, r in part.iterrows():
        if r["compound"] not in known:
            errs.append(ERR_PART_COMPOUND.format(rxn=r["reaction_id"], cpd=r["compound"]))
        if r["role"] not in ROLES:
            errs.append(ERR_ROLE.format(rxn=r["reaction_id"], cpd=r["compound"], role=r["role"]))

    for _, r in rxn.iterrows():
        if str(r["curated"]).strip().lower() in ("false", "0", "no") and _blank(r["rhea_id"]):
            errs.append(ERR_CURATED.format(rxn=r["reaction_id"]))

    for _, r in sr.iterrows():
        s, x, t = r["step_id"], r["reaction_id"], r["tier"]
        if t not in TIERS:
            errs.append(ERR_TIER.format(step=s, rxn=x, tier=t))
        for f in ("organism", "gene_name"):
            if _blank(r[f]):
                errs.append(ERR_EMPTY.format(step=s, rxn=x, field=f))
        if t in EVIDENCE_TIERS:
            for f in ("pmid", "quote"):
                if _blank(r[f]):
                    errs.append(ERR_EVIDENCE.format(step=s, rxn=x, tier=t, field=f))
    return errs


def path_possible(tables: dict[str, pd.DataFrame], step_ids: list[str]) -> pd.DataFrame:
    """Compounds reachable by following ``step_reaction`` links in the given step order.

    Returns DataFrame ``compound, min_tier, steps`` for the products of the LAST step.
    A link (step, reaction) from step k feeds step k+1 when a product of the first link is a
    substrate of the next link (``reaction_participant``). ``min_tier`` is the weakest tier along
    the path (order T1a < T1b < T1c < T2 < T3 < T4, weakest = highest index). Links with tier T4
    are dropped, so any path through a T4 link is excluded. When several paths reach the same
    compound, the one with the strongest weakest-tier is kept (ties: first found).

    Simplifications: only the declared step order is followed (no branching enumeration, no
    skipping of optional steps, no loops); the ``direction`` column is ignored (substrates feed
    forward as stored); the starting compounds are the substrates of the first step's reactions;
    cofactors listed as participants count like any other substrate or product.
    """
    sr, part = tables["step_reaction"], tables["reaction_participant"]
    subs: dict[str, set[str]] = {}
    prods: dict[str, set[str]] = {}
    for _, r in part.iterrows():
        d = subs if r["role"] == "substrate" else prods
        d.setdefault(r["reaction_id"], set()).add(r["compound"])

    # reach: compound -> (tier index, [step ids])
    reach: dict[str, tuple[int, list[str]]] | None = None
    for sid in step_ids:
        links = sr[sr["step_id"] == sid]
        nxt: dict[str, tuple[int, list[str]]] = {}
        for _, l in links.iterrows():
            if l["tier"] not in TIERS or l["tier"] == "T4":
                continue
            ti = TIERS.index(l["tier"])
            rid = l["reaction_id"]
            if reach is None:
                base, steps = ti, [sid]
            else:
                feeders = [reach[c] for c in subs.get(rid, ()) if c in reach]
                if not feeders:
                    continue
                best = min(feeders, key=lambda x: x[0])
                base, steps = max(best[0], ti), best[1] + [sid]
            for p in prods.get(rid, ()):
                if p not in nxt or base < nxt[p][0]:
                    nxt[p] = (base, steps)
        reach = nxt
        if not reach:
            break
    rows = []
    for c, (ti, steps) in sorted((reach or {}).items()):
        if len(steps) == len(step_ids):
            rows.append({"compound": c, "min_tier": TIERS[ti], "steps": ",".join(steps)})
    return pd.DataFrame(rows, columns=["compound", "min_tier", "steps"])
