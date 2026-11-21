"""Tests for the L-DOPA melanin / other melanin flat knowledge tables (stage 1, task S5)."""
import re
from pathlib import Path

import pandas as pd
import pytest

from pathway_census.defs import load_pathway
from pathway_census.kb import (EVIDENCE_TIERS, formula_mass, load_kb, path_possible,
                               species_key, validate_kb)

ROOT = Path(__file__).resolve().parent.parent
KB = ROOT / "kb" / "melanin_ldopa"
YAML = ROOT / "pathways" / "melanin_ldopa.yaml"
QUOTES = ROOT / "tests" / "data" / "kb" / "melanin_ldopa_quotes.tsv"
C = "CUR:melanin_ldopa:"
TYR, DOPA, DQ, DC = "CHEBI:17895", "CHEBI:15765", "CHEBI:16852", "CHEBI:15772"
HGA, HPP, ASPE, MAA = "CHEBI:44747", "CHEBI:15999", "CHEBI:17704", "CHEBI:47904"
DA, DACHR, DHICA = "CHEBI:18243", "CHEBI:189002", "CHEBI:2003"
T4_STEP = "dopachrome_onward"
NO_ENZYME = "not assigned (no enzyme identified)"


@pytest.fixture(scope="module")
def kb():
    return load_kb(KB)


def _parts(kb, rid, role):
    p = kb["reaction_participant"]
    return set(p[(p["reaction_id"] == rid) & (p["role"] == role)]["compound"])


def _routes():
    return {r.id: r.step_ids for r in load_pathway(YAML).routes}


def _links(kb, step):
    sr = kb["step_reaction"]
    sr = sr[sr["step_id"] == step]
    return {(r, a, p): t for r, a, p, t in zip(sr["reaction_id"], sr["protein_accession"],
                                              sr["pmid"], sr["tier"])}


def test_tables_validate(kb):
    assert validate_kb(kb) == []


def test_tables_not_empty(kb):
    for name, df in kb.items():
        assert len(df) > 0, name


def test_every_t1_t2_link_has_pmid_quote_and_assay(kb):
    sr = kb["step_reaction"]
    ev = sr[sr["tier"].isin(EVIDENCE_TIERS)]
    assert len(ev) > 0
    for _, r in ev.iterrows():
        assert r["pmid"].strip().isdigit(), (r["step_id"], r["reaction_id"])
        assert r["quote"].strip(), (r["step_id"], r["reaction_id"])
        assert r["assay"].strip(), (r["step_id"], r["reaction_id"])


def test_tiers_used(kb):
    assert set(kb["step_reaction"]["tier"]) == {"T1a", "T1b", "T1c", "T3", "T4"}


def test_one_quoted_sentence_per_link(kb):
    sr = kb["step_reaction"]
    assert not sr.duplicated(["step_id", "reaction_id", "protein_accession", "pmid"]).any()
    for _, r in sr.iterrows():
        q = r["quote"]
        if q:
            assert r["pmid"], (r["step_id"], r["protein_accession"])
            assert len(q) <= 400 and q == q.strip() and q.endswith(".")
            assert not re.search(r"[.!?] [A-Z]", q[:-1]), (r["pmid"], q)
    # PMID 15643074 gives two sentences: one names lac2, the other the lac1 mutants
    by = sr[sr["pmid"] == "15643074"]
    quotes = dict(zip(by["protein_accession"], by["quote"]))
    assert quotes["J9VQZ4"] != quotes["J9VY90"] and set(by["tier"]) == {"T3"}
    assert "lac2 mutation" in quotes["J9VQZ4"] and "lac1 mutants" in quotes["J9VY90"]


def test_every_link_names_organism_gene_and_accession_of_a_yaml_ref(kb):
    refs = {}
    for st in load_pathway(YAML).steps:
        for r in st.refs:
            refs[(st.id, r.accession)] = r
    sr = kb["step_reaction"]
    real = sr[sr["step_id"] != T4_STEP]
    for _, r in real.iterrows():
        ref = refs[(r["step_id"], r["protein_accession"])]
        assert r["organism"] == ref.species, (r["step_id"], r["protein_accession"])
        assert r["gene_name"] == ref.gene_name, (r["step_id"], r["protein_accession"])
    # every reference has at least one link
    assert {(s, a) for s, a in zip(real["step_id"], real["protein_accession"])} == set(refs)


def test_kb_steps_are_the_yaml_steps_plus_the_t4_step(kb):
    assert set(kb["step"]["step_id"]) == {s.id for s in load_pathway(YAML).steps} | {T4_STEP}


def test_yaml_step_reactions_are_the_kb_reactions(kb):
    sr = kb["step_reaction"]
    for st in load_pathway(YAML).steps:
        assert set(st.reactions) == set(sr[sr["step_id"] == st.id]["reaction_id"]), st.id


# ---- requirement 5: compounds ----

def test_compounds_of_requirement_5_with_kegg_cross_references(kb):
    c = kb["compound"].set_index("chebi_id")
    want = {TYR: ("C9H11NO3", "C00082"), DOPA: ("C9H11NO4", "C00355"),
            DQ: ("C9H9NO4", "C00822"), DC: ("C9H7NO4", "C01693"),
            HGA: ("C8H8O4", "C00544"), ASPE: ("C17H12O5", "C02006"),
            HPP: ("C9H8O4", "C01179"), MAA: ("C8H8O6", "C01036"),
            DA: ("C8H11NO2", "C03758"), DACHR: ("C8H7NO2", ""), DHICA: ("C9H7NO4", "C04185")}
    assert set(c.index) == set(want)
    for cid, (formula, kegg) in want.items():
        assert c.loc[cid, "formula"] == formula, cid
        assert c.loc[cid, "kegg_id"] == kegg, cid


def test_compounds_are_neutral_parents_not_the_rhea_species(kb):
    """Rhea lists zwitterions and anions (CHEBI:58315, 57504, 57924, 36242, 16169, 17105,
    58240, 57509, 16875). The table holds the neutral parents."""
    ids = set(kb["compound"]["chebi_id"])
    for charged in ("58315", "57504", "57924", "36242", "16169", "17105", "58240", "57509",
                    "16875", "193114", "193115"):
        assert f"CHEBI:{charged}" not in ids
    for _, r in kb["compound"].iterrows():
        assert not re.search(r"[+-]", r["formula"]), r["chebi_id"]
        assert "+" not in r["smiles"] and "-]" not in r["smiles"], r["chebi_id"]


def test_masses_are_computed_from_formulas(kb):
    for _, r in kb["compound"].iterrows():
        assert r["neutral_exact_mass"] == f"{formula_mass(r['formula']):.4f}", r["chebi_id"]
        assert r["inchikey"].strip() and r["inchikey_block1"] == r["inchikey"][:14]
    c = kb["compound"].set_index("chebi_id")
    # isomers: L-dopachrome and 5,6-dihydroxyindole-2-carboxylic acid have one formula
    assert c.loc[DC, "formula"] == c.loc[DHICA, "formula"]
    assert c.loc[DC, "inchikey_block1"] != c.loc[DHICA, "inchikey_block1"]


def test_kegg_ids_are_cross_references_only(kb):
    ko = dict(zip(kb["step"]["step_id"], kb["step"]["ko"]))
    assert ko == {"tyrosinase": "K00505", "laccase_diphenol": "K00421;K05909",
                  "hppd": "K00457", "hmga": "K00451", "aspulvinone_e_synthase": "",
                  "asp_melanin_tyrosinase": "", T4_STEP: ""}
    assert (kb["reaction"]["kegg_rn"] == "").all()
    assert list(kb["compound"].columns) == ["chebi_id", "inchikey", "inchikey_block1", "smiles",
                                            "formula", "neutral_exact_mass", "kegg_id",
                                            "np_class"]


def test_curated_reaction_ids(kb):
    for _, r in kb["reaction"].iterrows():
        if r["curated"] == "true":
            assert r["reaction_id"].startswith(C), r["reaction_id"]
            assert r["rhea_id"] == ""
        else:
            assert r["reaction_id"] == r["rhea_id"] and r["rhea_id"].startswith("RHEA:")
    assert set(kb["reaction"]["reaction_id"]) == {
        "RHEA:18117", "RHEA:34287", C + "ldopa_oxidation_laccase",
        C + "dopamine_oxidation_laccase", C + "hga_oxidation_laccase",
        C + "dopaquinone_cyclisation", "RHEA:13041", "RHEA:16189", "RHEA:15449",
        C + "hga_autooxidation", "RHEA:63824", "RHEA:74195", "RHEA:74199"}
    ec = dict(zip(kb["reaction"]["reaction_id"], kb["reaction"]["ec"]))
    assert ec["RHEA:18117"] == ec["RHEA:34287"] == "1.14.18.1"
    assert ec[C + "ldopa_oxidation_laccase"] == ec[C + "dopamine_oxidation_laccase"] == "1.10.3.2"
    assert ec["RHEA:16189"] == "1.13.11.27" and ec["RHEA:15449"] == "1.13.11.5"
    assert ec["RHEA:63824"] == ec["RHEA:74195"] == ec["RHEA:74199"] == ""
    assert ec["RHEA:13041"] == "5.3.3.12" and ec[C + "dopaquinone_cyclisation"] == ""


# ---- requirement 1: tyrosinase and laccase are separate steps ----

def test_tyrosinase_reactions(kb):
    assert _parts(kb, "RHEA:18117", "substrate") == {TYR}
    assert _parts(kb, "RHEA:18117", "product") == {DQ}
    assert _parts(kb, "RHEA:34287", "substrate") == {DOPA}
    assert _parts(kb, "RHEA:34287", "product") == {DQ}
    # free L-DOPA is not a product of any reaction of the table (pathway review, O7)
    p = kb["reaction_participant"]
    assert DOPA not in set(p[p["role"] == "product"]["compound"])


def test_laccase_step_has_no_reaction_whose_substrate_is_tyrosine(kb):
    sr = kb["step_reaction"]
    rids = set(sr[sr["step_id"] == "laccase_diphenol"]["reaction_id"])
    assert rids == {C + "ldopa_oxidation_laccase", C + "dopamine_oxidation_laccase",
                    C + "hga_oxidation_laccase"}
    for rid in rids:
        assert TYR not in _parts(kb, rid, "substrate"), rid
        assert TYR not in _parts(kb, rid, "product"), rid
    assert _parts(kb, C + "ldopa_oxidation_laccase", "substrate") == {DOPA}
    assert _parts(kb, C + "ldopa_oxidation_laccase", "product") == {DQ}
    assert _parts(kb, C + "dopamine_oxidation_laccase", "substrate") == {DA}
    assert _parts(kb, C + "dopamine_oxidation_laccase", "product") == {DACHR}
    assert _parts(kb, C + "hga_oxidation_laccase", "substrate") == {HGA}
    assert _parts(kb, C + "hga_oxidation_laccase", "product") == set()
    # tyrosine is a substrate of the tyrosinase reaction only
    p = kb["reaction_participant"]
    assert set(p[p["compound"] == TYR]["reaction_id"]) == {"RHEA:18117"}
    # no link joins a laccase to a tyrosinase reaction or a tyrosinase to a laccase reaction
    tyr_rx = set(sr[sr["step_id"] == "tyrosinase"]["reaction_id"])
    assert tyr_rx == {"RHEA:18117", "RHEA:34287"} and not (tyr_rx & rids)


def test_laccase_reactions_are_curated_and_say_why(kb):
    rx = kb["reaction"].set_index("reaction_id")
    for name in ("ldopa_oxidation_laccase", "dopamine_oxidation_laccase", "hga_oxidation_laccase"):
        assert rx.loc[C + name, "curated"] == "true"
        assert "RHEA:11276" in rx.loc[C + name, "note"], name
    assert "RHEA:34287" in rx.loc[C + "ldopa_oxidation_laccase", "note"]
    assert "exogenous" in rx.loc[C + "ldopa_oxidation_laccase", "note"]


def test_tyrosinase_tiers_are_bound_to_protein_and_organism(kb):
    assert _links(kb, "tyrosinase") == {
        ("RHEA:18117", "P00440", "2147580"): "T1c",
        ("RHEA:34287", "P00440", ""): "T3",
        ("RHEA:18117", "A0A1S9DK56", "23053534"): "T3",
        ("RHEA:34287", "A0A1S9DK56", "22213164"): "T3",
        ("RHEA:18117", "Q00234", "7893753"): "T3",
        ("RHEA:34287", "Q00234", "7893753"): "T3",
        ("RHEA:18117", "C7FF05", "28500345"): "T3",
        ("RHEA:34287", "C7FF05", "28500345"): "T3",
        ("RHEA:18117", "C7FF04", "21598903"): "T3",
        ("RHEA:34287", "C7FF04", "21598903"): "T3",
        ("RHEA:18117", "A7BHQ9", "17617709"): "T3",
        ("RHEA:34287", "A7BHQ9", "17617709"): "T3",
        ("RHEA:18117", "Q2TL94", "16151802"): "T3",
        ("RHEA:34287", "Q2TL94", "16151802"): "T3",
    }


# Binding rule 1, applied as written (fix round 1): a T1 link needs a quoted sentence that names
# the protein of the link. The names below are the words of the sentence that name it.
NAMED = {
    ("P00440", "2147580"): "tyrosinase",        # T1c: transformant with high tyrosinase levels
    ("Q55P57", "8760791"): "CNLAC1",            # T1c
    ("Q4WHU1", "19028908"): "hppD",             # T1b
    ("Q4WHT9", "19028908"): "hmgA",             # T1b
    ("Q00667", "7673153"): "HmgA",              # T1b
    ("Q0CRX1", "27133313"): "MelA",             # T1a
    ("Q0CRX1", "28791090"): "atmelA",           # T1a
}


def test_every_t1_quote_names_the_protein_of_its_link(kb):
    sr = kb["step_reaction"]
    ev = sr[sr["tier"].isin(EVIDENCE_TIERS)]
    assert {(a, p) for a, p in zip(ev["protein_accession"], ev["pmid"])} == set(NAMED)
    for _, r in ev.iterrows():
        assert NAMED[(r["protein_accession"], r["pmid"])] in r["quote"], r["protein_accession"]
    assert dict(ev["tier"].value_counts()) == {"T1a": 2, "T1b": 3, "T1c": 2}
    # T1a: the sentence states the product for the named protein
    for q in ev[ev["tier"] == "T1a"]["quote"]:
        assert "aspulvinone E" in q
    # T1b: the sentence states the intermediate
    for q in ev[ev["tier"] == "T1b"]["quote"]:
        assert "homogentisic acid" in q or "homogentisate" in q


def test_no_tyrosinase_link_is_above_t1c_and_held_links_say_why(kb):
    """No quoted tyrosinase sentence names the protein together with L-tyrosine or L-DOPA.
    The substrate-class rule of the first version is withdrawn."""
    sr = kb["step_reaction"]
    tyr = sr[sr["step_id"] == "tyrosinase"]
    assert dict(tyr["tier"].value_counts()) == {"T3": 13, "T1c": 1}
    held = tyr[(tyr["tier"] == "T3") & (tyr["quote"] != "")]
    assert len(held) == 12
    for _, r in held.iterrows():
        assert r["assay"].startswith("Held at T3."), r["protein_accession"]
    readme = (ROOT / "kb" / "README.md").read_text(encoding="utf-8")
    assert "substrate class" not in readme


def test_laccase_tiers(kb):
    assert _links(kb, "laccase_diphenol") == {
        (C + "ldopa_oxidation_laccase", "Q55P57", "8760791"): "T1c",
        (C + "dopamine_oxidation_laccase", "Q55P57", "8300520"): "T3",
        (C + "hga_oxidation_laccase", "Q55P57", "17098915"): "T3",
        (C + "ldopa_oxidation_laccase", "J9VY90", "15643074"): "T3",
        (C + "hga_oxidation_laccase", "J9VY90", "17098915"): "T3",
        (C + "ldopa_oxidation_laccase", "J9VQZ4", "15643074"): "T3",
    }
    sr = kb["step_reaction"]
    dop = sr[sr["reaction_id"] == C + "dopamine_oxidation_laccase"].iloc[0]
    assert dop["quote"] == ("Dopamine was shown by mass spectroscopy to be oxidized to decarboxy "
                            "dopachrome, an intermediate of melanin synthesis.")
    assert dop["assay"].startswith("Held at T3.") and "does not name the enzyme" in dop["assay"]
    # the strain of PMID 15643074 is not confirmed: the H99 links are held at T3 (I3)
    h99 = sr[sr["pmid"] == "15643074"]
    assert len(h99) == 2
    for a in h99["assay"]:
        assert a.startswith("Held at T3.") and "names no strain" in a
    hga = sr[sr["reaction_id"] == C + "hga_oxidation_laccase"]
    assert set(hga["quote"]) == {
        "HGA melanization was laccase dependent and repressed by high concentrations of glucose."}
    for a in hga["assay"]:
        assert a.startswith("Held at T3.") and "does not name the laccase gene" in a


def test_pyomelanin_and_asp_melanin_tiers(kb):
    assert _links(kb, "hppd") == {("RHEA:16189", "Q4WHU1", "19028908"): "T1b"}
    assert _links(kb, "hmga") == {("RHEA:15449", "Q4WHT9", "19028908"): "T1b",
                                  ("RHEA:15449", "Q00667", "7673153"): "T1b"}
    q00667 = kb["step_reaction"][kb["step_reaction"]["protein_accession"] == "Q00667"].iloc[0]
    assert q00667["quote"] == (
        "The HmgA protein catalyzes an essential step in phenylalanine catabolism, and disruption "
        "of the gene results in accumulation of homogentisate in broths containing phenylalanine.")
    assert _links(kb, "aspulvinone_e_synthase") == {
        ("RHEA:63824", "Q0CRX1", "27133313"): "T1a",
        ("RHEA:63824", "Q0CRX1", "28791090"): "T1a"}
    # The abstract of PMID 27133313 states that TyrP activates aspulvinone E for polymerisation;
    # it states no assay and no product of TyrP. Held at T3 with the sentence as context.
    assert _links(kb, "asp_melanin_tyrosinase") == {
        ("RHEA:74195", "Q0CRX0", "27133313"): "T3",
        ("RHEA:74199", "Q0CRX0", "27133313"): "T3"}
    sr = kb["step_reaction"]
    tyrp = sr[sr["step_id"] == "asp_melanin_tyrosinase"]
    assert set(tyrp["quote"]) == {"We demonstrate that MelA produces aspulvinone E, which is "
                                  "activated for polymerization by TyrP."}
    for a in tyrp["assay"]:
        assert a.startswith("Held at T3.")
    assert _parts(kb, "RHEA:16189", "substrate") == {HPP}
    assert _parts(kb, "RHEA:16189", "product") == {HGA}
    assert _parts(kb, "RHEA:15449", "substrate") == {HGA}
    assert _parts(kb, "RHEA:15449", "product") == {MAA}
    assert _parts(kb, "RHEA:63824", "substrate") == {HPP}
    assert _parts(kb, "RHEA:63824", "product") == {ASPE}
    # the TyrP products have charged ChEBI entries only: no product participant
    for rid in ("RHEA:74195", "RHEA:74199"):
        assert _parts(kb, rid, "substrate") == {ASPE} and _parts(kb, rid, "product") == set()
        note = kb["reaction"].set_index("reaction_id").loc[rid, "note"]
        assert "CHEBI:19311" in note and "no neutral" in note


def test_dopaquinone_onward_is_non_enzymatic_and_t4_for_any_fungal_enzyme(kb):
    rid = C + "dopaquinone_cyclisation"
    assert _parts(kb, rid, "substrate") == {DQ} and _parts(kb, rid, "product") == {DC}
    assert _parts(kb, "RHEA:13041", "substrate") == {DC}
    assert _parts(kb, "RHEA:13041", "product") == {DHICA}
    rx = kb["reaction"].set_index("reaction_id")
    assert "non-enzymatic" in rx.loc[rid, "note"]
    sr = kb["step_reaction"]
    onward = sr[sr["reaction_id"].isin([rid, "RHEA:13041"])]
    assert set(onward["step_id"]) == {T4_STEP} and set(onward["tier"]) == {"T4"}
    assert set(onward["organism"]) == {NO_ENZYME} and set(onward["protein_accession"]) == {""}
    assert set(sr[sr["tier"] == "T4"]["step_id"]) == {T4_STEP}
    # a path through a T4 link is never reported
    assert path_possible(kb, ["tyrosinase", T4_STEP]).empty
    # homogentisate autooxidation: curated, no enzyme link
    assert C + "hga_autooxidation" not in set(sr["reaction_id"])
    assert _parts(kb, C + "hga_autooxidation", "substrate") == {HGA}
    assert _parts(kb, C + "hga_autooxidation", "product") == set()


# Expected final products of every declared route (path_possible over the route's step order):
# compound -> (weakest tier along the best path, organisms of the links on that path, number of
# organism lists with that tier and that size; the list shown is the alphabetically first one).
# An empty dict means: the route cannot be followed to a product in the tables (UNVERIFIED.md, D).
AB, AO, TS = "agaricus bisporus", "aspergillus oryzae", "trametes sanguinea"
CD, CN = "cryptococcus deneoformans", "cryptococcus neoformans"
AF, AN, AT = "aspergillus fumigatus", "aspergillus nidulans", "aspergillus terreus"
ROUTE_PRODUCTS = {
    # Fix round 1 (binding rule 1 as written). T1c for Neurospora crassa only (colour phenotype
    # of an overexpression strain); every other tyrosinase link is T3.
    "tyrosinase_dopa": {DQ: ("T1c", ["neurospora crassa"], 1)},
    # L-dopaquinone at T1c (colour phenotype) for Cryptococcus deneoformans only; the H99 links
    # are T3. Dopaminechrome at T3 (the sentence does not name the enzyme).
    "laccase_exogenous": {DQ: ("T1c", [CD], 1), DACHR: ("T3", [CD], 1)},
    "pyomelanin": {HGA: ("T1b", [AF], 1)},
    "tyrosine_degradation": {MAA: ("T1b", [AF], 1)},
    # the TyrP products have no neutral ChEBI entry
    "asp_melanin": {},
}


def _have(got):
    return {c: (t, o.split(","), single, n) for c, t, o, single, n in
            zip(got["compound"], got["min_tier"], got["organisms"], got["single_organism"],
                got["n_equal_paths"])}


def test_every_declared_route_gives_its_expected_products(kb):
    routes = _routes()
    assert set(routes) == set(ROUTE_PRODUCTS)
    for rid, want in ROUTE_PRODUCTS.items():
        got = path_possible(kb, routes[rid])
        assert _have(got) == {c: (t, o, len(o) == 1, n) for c, (t, o, n) in want.items()}, rid
        # every route tier of this pathway holds for one organism
        assert got["single_organism"].all(), rid


def test_organisms_of_the_tables(kb):
    sr = kb["step_reaction"]
    organisms = sorted({species_key(o) for o in sr[sr["tier"] != "T4"]["organism"]})
    assert organisms == sorted([AB, AF, AN, AO, AT, CD, CN, "neurospora crassa",
                                "pholiota microspora", TS])


def test_routes_one_organism_at_a_time(kb):
    """What each organism reaches with its own links (the filter is applied before the search)."""
    routes = _routes()
    want = {
        "tyrosinase_dopa": {AB: ("T3", DQ), AO: ("T3", DQ), TS: ("T3", DQ),
                            "neurospora crassa": ("T1c", DQ), "pholiota microspora": ("T3", DQ)},
        "laccase_exogenous": {CD: ("T1c", DQ), CN: ("T3", DQ)},
        "pyomelanin": {AF: ("T1b", HGA)},
        "tyrosine_degradation": {AF: ("T1b", MAA)},
        "asp_melanin": {},
    }
    sr = kb["step_reaction"]
    organisms = sorted({species_key(o) for o in sr[sr["tier"] != "T4"]["organism"]})
    for rid, steps in routes.items():
        got = {}
        for org in organisms:
            out = path_possible(kb, steps, organism=org)
            if not out.empty:
                tiers = dict(zip(out["compound"], out["min_tier"]))
                cpd = want[rid][org][1] if org in want[rid] else sorted(tiers)[0]
                got[org] = (tiers[cpd], cpd)
                assert out["single_organism"].all()
        assert got == want[rid], rid
    # dopaminechrome: Cryptococcus deneoformans only
    out = path_possible(kb, routes["laccase_exogenous"], organism="Cryptococcus neoformans")
    assert set(out["compound"]) == {DQ}


def test_asp_melanin_route_reaches_aspulvinone_e_then_stops(kb):
    """Route asp_melanin gives no product: the TyrP products have charged ChEBI entries only.
    The first step alone gives aspulvinone E at T1a for Aspergillus terreus."""
    steps = _routes()["asp_melanin"]
    got = path_possible(kb, steps[:1])
    assert _have(got) == {ASPE: ("T1a", [AT], True, 1)}
    assert path_possible(kb, steps).empty


def test_no_route_makes_free_l_dopa_and_cryptococcus_has_no_tyrosine_reaction(kb):
    for rid, steps in _routes().items():
        for n in range(1, len(steps) + 1):
            assert DOPA not in set(path_possible(kb, steps[:n])["compound"]), rid
    sr = kb["step_reaction"]
    crypto = sr[sr["organism"].str.startswith("Cryptococcus")]
    assert len(crypto) == 6 and set(crypto["step_id"]) == {"laccase_diphenol"}
    for rid in set(crypto["reaction_id"]):
        assert TYR not in _parts(kb, rid, "substrate")


def test_organism_filter_accepts_the_names_used_in_the_papers(kb):
    # serotype D: NCBI gives "Cryptococcus neoformans var. neoformans" the taxid of
    # Cryptococcus deneoformans (40410)
    assert species_key("Cryptococcus neoformans var. neoformans") == CD
    assert species_key("Cryptococcus deneoformans (strain B-3501A)") == CD
    assert species_key("Cryptococcus neoformans var. grubii") == CN
    assert species_key("Cryptococcus neoformans (strain H99 / ATCC 208821 / CBS 10515 / "
                       "FGSC 9487)") == CN
    assert CD != CN
    assert species_key("Pholiota nameko") == species_key("Pholiota microspora")
    # "Pycnoporus sanguineus" (the name in PMID 16151802) has no NCBI taxid hit: no alias
    assert species_key("Pycnoporus sanguineus") != species_key("Trametes sanguinea")


def test_quotes_equal_the_committed_fixture(kb):
    """Change detector: any edit of a quote (or of the links that carry one) fails here.
    The fixture is not a check against PubMed; Task 7 does that."""
    cols = ["step_id", "reaction_id", "protein_accession", "pmid", "quote", "tier", "assay"]
    fix = pd.read_csv(QUOTES, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8",
                      quoting=3)
    assert list(fix.columns) == cols
    cur = kb["step_reaction"][cols].reset_index(drop=True)
    assert len(cur) == 29 and (fix["quote"] != "").sum() == 26
    assert cur.to_dict("records") == fix.to_dict("records")


def test_kb_readme_documents_the_ldopa_tables():
    text = (ROOT / "kb" / "README.md").read_text(encoding="utf-8")
    assert "kb/melanin_ldopa/" in text and "exogenous_substrate" in text
    assert "`cryptococcus deneoformans`" in text and "40410" in text
