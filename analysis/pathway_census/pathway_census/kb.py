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
ROLES = ("substrate", "product", "cofactor")
DIRECTIONS = ("forward", "reversible", "unknown")
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
ERR_CURATED_VAL = "reaction[{rxn}]: curated '{v}' not in true,false"
ERR_DIRECTION = "reaction[{rxn}]: direction '{v}' not in " + ",".join(DIRECTIONS)
ERR_STOICH = "reaction_participant[{rxn},{cpd}]: stoich '{v}' is not a positive number"
ERR_FK_STEP = "step_reaction[{step},{rxn}]: step_id not in step table"
ERR_FK_RXN = "step_reaction[{step},{rxn}]: reaction_id not in reaction table"
ERR_PART_RXN = "reaction_participant[{rxn},{cpd}]: reaction_id not in reaction table"
ERR_DUP = "{table}[{key}]: duplicate key"
ERR_PMID = "step_reaction[{step},{rxn}]: pmid '{v}' is not digits only"
ERR_INCHIKEY_EMPTY = "compound[{cpd}]: inchikey is empty but inchikey_block1 is '{b1}'"

# Monoisotopic masses (mass of the most abundant isotope, unified atomic mass units).
# Source: NIST Atomic Weights and Isotopic Compositions for All Elements
# (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl), fetched 2026-11-21, isotopes
# 1H, 12C, 14N, 16O, 23Na, 31P, 32S, 35Cl, 39K, 56Fe.
MONO_MASS = {
    "H": 1.00782503223,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "Na": 22.9897692820,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "Cl": 34.968852682,
    "K": 38.9637064864,
    "Fe": 55.93493633,
}

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_BRACKET_ATOM = re.compile(r"\[([^\]]*)\]")
_CHARGE = re.compile(r"(\++|-+)$|([+-])(\d+)$")


def net_charge(smiles: str) -> int:
    """Sum of the charges of all bracket atoms in a SMILES (``+``, ``++``, ``+2``, ``-``, ``-2``)."""
    total = 0
    for inner in _BRACKET_ATOM.findall(smiles):
        inner = re.sub(r":\d+$", "", inner)  # atom-map number
        m = _CHARGE.search(inner)
        if not m:
            continue
        if m.group(1):
            total += len(m.group(1)) * (1 if m.group(1)[0] == "+" else -1)
        else:
            total += int(m.group(3)) * (1 if m.group(2) == "+" else -1)
    return total


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
        if not ik.strip() and b1.strip():
            errs.append(ERR_INCHIKEY_EMPTY.format(cpd=cpd, b1=b1))
        elif b1 != ik[:14]:
            errs.append(ERR_BLOCK1.format(cpd=cpd, got=b1, want=ik[:14]))
        formula = str(r["formula"])
        charged = bool(re.search(r"[+-]", formula)) or net_charge(str(r["smiles"])) != 0
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
    rxn_ids, step_ids = set(rxn["reaction_id"]), set(tables["step"]["step_id"])
    for _, r in part.iterrows():
        k = dict(rxn=r["reaction_id"], cpd=r["compound"])
        if r["compound"] not in known:
            errs.append(ERR_PART_COMPOUND.format(**k))
        if r["reaction_id"] not in rxn_ids:
            errs.append(ERR_PART_RXN.format(**k))
        if r["role"] not in ROLES:
            errs.append(ERR_ROLE.format(role=r["role"], **k))
        try:
            ok = float(r["stoich"]) > 0
        except (TypeError, ValueError):
            ok = False
        if not ok:
            errs.append(ERR_STOICH.format(v=r["stoich"], **k))

    for _, r in rxn.iterrows():
        cur = str(r["curated"])
        if cur not in ("true", "false"):
            errs.append(ERR_CURATED_VAL.format(rxn=r["reaction_id"], v=cur))
        elif cur == "false" and _blank(r["rhea_id"]):
            errs.append(ERR_CURATED.format(rxn=r["reaction_id"]))
        if r["direction"] not in DIRECTIONS:
            errs.append(ERR_DIRECTION.format(rxn=r["reaction_id"], v=r["direction"]))

    for _, r in sr.iterrows():
        s, x, t = r["step_id"], r["reaction_id"], r["tier"]
        if s not in step_ids:
            errs.append(ERR_FK_STEP.format(step=s, rxn=x))
        if x not in rxn_ids:
            errs.append(ERR_FK_RXN.format(step=s, rxn=x))
        if t not in TIERS:
            errs.append(ERR_TIER.format(step=s, rxn=x, tier=t))
        for f in ("organism", "gene_name"):
            if _blank(r[f]):
                errs.append(ERR_EMPTY.format(step=s, rxn=x, field=f))
        if t in EVIDENCE_TIERS:
            for f in ("pmid", "quote"):
                if _blank(r[f]):
                    errs.append(ERR_EVIDENCE.format(step=s, rxn=x, tier=t, field=f))
        if not _blank(r["pmid"]) and not str(r["pmid"]).strip().isdigit():
            errs.append(ERR_PMID.format(step=s, rxn=x, v=r["pmid"]))

    keys = {"compound": ["chebi_id"], "reaction": ["reaction_id"],
            "reaction_participant": ["reaction_id", "compound", "role"], "step": ["step_id"],
            "step_reaction": ["step_id", "reaction_id", "organism", "gene_name"]}
    for name, cols in keys.items():
        dup = tables[name][tables[name].duplicated(cols, keep="first")]
        for _, r in dup.iterrows():
            errs.append(ERR_DUP.format(table=name, key=",".join(str(r[c]) for c in cols)))
    return errs


def path_possible(tables: dict[str, pd.DataFrame], step_ids: list[str]) -> pd.DataFrame:
    """Compounds reachable by following ``step_reaction`` links in the given step order.

    Returns DataFrame ``compound, min_tier, steps`` for the products of the LAST step.
    A link (step, reaction) of step k+1 is usable when ALL non-cofactor substrates of its reaction
    are products of the previous step's usable links (``reaction_participant``). Cofactor rows are
    ignored on both sides. ``min_tier`` is the weakest tier along the path (order
    T1a < T1b < T1c < T2 < T3 < T4, weakest = highest index; across the substrates of one link
    the weakest also counts). Links with tier T4 are dropped, so a path through a T4 link is
    excluded. When several links give the same product, the strongest weakest-tier is kept.
    ``steps`` is the comma-joined declared order. Raises ValueError for a step_id absent from the
    ``step`` table.

    Simplifications: only the declared step order is followed (no branching enumeration, no
    skipping of optional steps, no loops). Roles must be stored in the direction the step uses;
    ``direction`` is informational and ignored here. The starting compounds are the substrates
    of the first step's reactions, with no requirement on them.
    """
    known = set(tables["step"]["step_id"])
    for sid in step_ids:
        if sid not in known:
            raise ValueError(f"unknown step_id {sid!r}")
    sr, part = tables["step_reaction"], tables["reaction_participant"]
    subs: dict[str, set[str]] = {}
    prods: dict[str, set[str]] = {}
    for _, r in part.iterrows():
        if r["role"] == "substrate":
            subs.setdefault(r["reaction_id"], set()).add(r["compound"])
        elif r["role"] == "product":
            prods.setdefault(r["reaction_id"], set()).add(r["compound"])

    reach: dict[str, int] | None = None  # compound -> weakest tier index so far
    for sid in step_ids:
        nxt: dict[str, int] = {}
        for _, l in sr[sr["step_id"] == sid].iterrows():
            if l["tier"] not in TIERS or l["tier"] == "T4":
                continue
            base = TIERS.index(l["tier"])
            rid = l["reaction_id"]
            if reach is not None:
                need = subs.get(rid, set())
                if not need or not need <= reach.keys():
                    continue
                base = max([base] + [reach[c] for c in need])
            for p in prods.get(rid, ()):
                if p not in nxt or base < nxt[p]:
                    nxt[p] = base
        reach = nxt
        if not reach:
            break
    steps = ",".join(step_ids)
    rows = [{"compound": c, "min_tier": TIERS[ti], "steps": steps}
            for c, ti in sorted((reach or {}).items())]
    return pd.DataFrame(rows, columns=["compound", "min_tier", "steps"])
