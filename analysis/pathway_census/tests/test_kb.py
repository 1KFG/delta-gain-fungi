from pathlib import Path

import pandas as pd
import pytest

from pathway_census import kb

KB = Path(__file__).parent / "data" / "kb"


@pytest.fixture
def tables():
    return kb.load_kb(KB)


def _set(tables, name, mask_col, mask_val, col, val):
    df = tables[name].copy()
    df.loc[df[mask_col] == mask_val, col] = val
    tables[name] = df
    return tables


def test_load_kb_has_five_tables(tables):
    assert set(tables) == {"compound", "reaction", "reaction_participant", "step", "step_reaction"}
    assert list(tables["compound"].columns) == kb.COLUMNS["compound"]


def test_valid_fixture_passes(tables):
    assert kb.validate_kb(tables) == []


def test_tiers_constant():
    assert kb.TIERS == ("T1a", "T1b", "T1c", "T2", "T3", "T4")


def test_bad_tier(tables):
    _set(tables, "step_reaction", "step_id", "s1", "tier", "T9")
    assert kb.validate_kb(tables) == ["step_reaction[s1,R1]: tier 'T9' not in T1a,T1b,T1c,T2,T3,T4"]


def test_empty_organism(tables):
    _set(tables, "step_reaction", "step_id", "s1", "organism", "")
    assert kb.validate_kb(tables) == ["step_reaction[s1,R1]: organism is empty"]


def test_empty_gene_name(tables):
    _set(tables, "step_reaction", "step_id", "s2", "gene_name", "")
    assert kb.validate_kb(tables) == ["step_reaction[s2,R2]: gene_name is empty"]


@pytest.mark.parametrize("tier", ["T1a", "T1b", "T1c", "T2"])
def test_pmid_quote_required(tables, tier):
    _set(tables, "step_reaction", "step_id", "s1", "tier", tier)
    _set(tables, "step_reaction", "step_id", "s1", "pmid", "")
    _set(tables, "step_reaction", "step_id", "s1", "quote", "")
    assert kb.validate_kb(tables) == [
        f"step_reaction[s1,R1]: tier '{tier}' requires pmid",
        f"step_reaction[s1,R1]: tier '{tier}' requires quote",
    ]


def test_t3_needs_no_pmid(tables):
    assert kb.validate_kb(tables) == []  # s3 is T3 with empty pmid and quote


def test_participant_unknown_compound(tables):
    _set(tables, "reaction_participant", "reaction_id", "R1", "compound", "CHEBI:ZZ")
    errs = kb.validate_kb(tables)
    assert "reaction_participant[R1,CHEBI:ZZ]: compound not in compound table" in errs


def test_participant_bad_role(tables):
    _set(tables, "reaction_participant", "compound", "CHEBI:B", "role", "catalyst")
    errs = kb.validate_kb(tables)
    assert errs.count("reaction_participant[R1,CHEBI:B]: role 'catalyst' not in substrate,product,cofactor") == 1


def test_block1_mismatch(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "inchikey_block1", "XXXXXXXXXXXXXX")
    assert kb.validate_kb(tables) == [
        "compound[CHEBI:A]: inchikey_block1 'XXXXXXXXXXXXXX' != first 14 characters of inchikey 'AAAAAAAAAAAAAA'"
    ]


def test_charged_formula(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "formula", "C10H7O4-")
    assert kb.validate_kb(tables) == ["compound[CHEBI:A]: charged species not allowed, use the neutral parent"]


def test_charged_smiles(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "smiles", "CC([O-])=O")
    assert kb.validate_kb(tables) == ["compound[CHEBI:A]: charged species not allowed, use the neutral parent"]


def test_mass_disagrees(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "neutral_exact_mass", "192.1")
    errs = kb.validate_kb(tables)
    assert errs == ["compound[CHEBI:A]: neutral_exact_mass 192.1 differs from formula mass 192.0423 by more than 0.0005"]


def test_unparsable_formula(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "formula", "C10Xx2")
    errs = kb.validate_kb(tables)
    assert errs == ["compound[CHEBI:A]: formula 'C10Xx2' cannot be parsed (unsupported element 'Xx')"]


def test_extra_compound_column(tables):
    tables["compound"] = tables["compound"].assign(kegg_pubchem_cid="1")
    assert kb.validate_kb(tables) == ["compound: unexpected column 'kegg_pubchem_cid' (kegg_id is informational only)"]


def test_missing_column(tables):
    tables["reaction"] = tables["reaction"].drop(columns=["rhea_id"])
    assert "reaction: missing column 'rhea_id'" in kb.validate_kb(tables)


def test_uncurated_needs_rhea(tables):
    _set(tables, "reaction", "reaction_id", "R3", "rhea_id", "")
    assert kb.validate_kb(tables) == ["reaction[R3]: curated=false requires rhea_id"]


def test_formula_mass_values():
    assert kb.formula_mass("C40H52O4") == pytest.approx(596.3866, abs=0.0005)
    assert kb.formula_mass("C10H8O4") == pytest.approx(192.0423, abs=0.0005)
    assert kb.formula_mass("C9H11NO4") == pytest.approx(197.0688, abs=0.0005)


def test_formula_mass_other_elements():
    # NaCl: 22.98977 + 34.96885
    assert kb.formula_mass("NaCl") == pytest.approx(57.95862, abs=0.0005)


def test_formula_mass_unknown_element():
    with pytest.raises(ValueError, match="unsupported element 'Br'"):
        kb.formula_mass("C2H5Br")


def test_formula_mass_garbage():
    with pytest.raises(ValueError):
        kb.formula_mass("c2h5")


def test_path_possible_three_steps(tables):
    out = kb.path_possible(tables, ["s1", "s2", "s3"])
    assert list(out.columns) == ["compound", "min_tier", "steps"]
    assert out["compound"].tolist() == ["CHEBI:D"]
    assert out["min_tier"].tolist() == ["T3"]
    assert out["steps"].tolist() == ["s1,s2,s3"]


def test_path_possible_two_steps_weakest(tables):
    out = kb.path_possible(tables, ["s1", "s2"])
    assert out["compound"].tolist() == ["CHEBI:C"]
    assert out["min_tier"].tolist() == ["T1c"]


def test_path_possible_excludes_t4(tables):
    _set(tables, "step_reaction", "step_id", "s2", "tier", "T4")
    assert kb.path_possible(tables, ["s1", "s2"]).empty
    assert kb.path_possible(tables, ["s1", "s2", "s3"]).empty


def test_path_possible_broken_chain(tables):
    assert kb.path_possible(tables, ["s1", "s3"]).empty


def test_path_possible_alternative_links_pick_strongest(tables):
    extra = pd.DataFrame([dict(step_id="s2", reaction_id="R2", organism="Synthetica three",
                               gene_name="geneB2", protein_accession="ACC9", tier="T1a",
                               assay="in vitro", pmid="3", quote="synthetic quote")])
    tables["step_reaction"] = pd.concat([tables["step_reaction"], extra], ignore_index=True)
    out = kb.path_possible(tables, ["s1", "s2"])
    assert out["min_tier"].tolist() == ["T1a"]


# ---- fix round 1 -------------------------------------------------------------------------
import json  # noqa: E402

CHARGE_MSG = "compound[CHEBI:A]: charged species not allowed, use the neutral parent"


def _smiles_fixture():
    return json.loads((KB / "chebi_smiles_charge.json").read_text())


def _with_smiles(tables, smiles):
    return _set(tables, "compound", "chebi_id", "CHEBI:A", "smiles", smiles)


@pytest.mark.parametrize("cid", ["CHEBI:57504", "CHEBI:27798"])
def test_neutral_zwitterion_and_nitro_from_chebi_pass(tables, cid):
    _with_smiles(tables, _smiles_fixture()[cid]["smiles"])
    assert kb.validate_kb(tables) == []


def test_literal_nitro_passes(tables):
    _with_smiles(tables, "[N+](=O)[O-]")
    assert kb.validate_kb(tables) == []


@pytest.mark.parametrize("smiles", ["CC([O-])=O", "[Na+]", "[NH3+]C", "[Fe+2]", "[O--]"])
def test_net_charged_smiles_fail(tables, smiles):
    _with_smiles(tables, smiles)
    assert kb.validate_kb(tables) == [CHARGE_MSG]


def test_net_charge_parser():
    assert kb.net_charge("[NH3+][C@@H](C)C(=O)[O-]") == 0
    assert kb.net_charge("[N++]") == 2
    assert kb.net_charge("[Fe+2]") == 2
    assert kb.net_charge("[O-2]") == -2
    assert kb.net_charge("[C@@H]") == 0


def test_plus_in_formula_fails(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "formula", "C10H9O4+")
    assert kb.validate_kb(tables) == [CHARGE_MSG]


@pytest.mark.parametrize("formula,chebi_num", [("H3O4P", "26078"), ("H2O4S", "26836"),
                                               ("C10H10Fe", "30672"), ("KCl", "32588")])
def test_formula_mass_vs_stored_chebi(formula, chebi_num):
    # Stored ChEBI responses (phosphoric acid, sulfuric acid, ferrocene, KCl), captured 2026-11-21.
    d = json.loads((KB / f"chebi_{chebi_num}.json").read_text())
    mono = float(d["chemical_data"]["monoisotopic_mass"])
    assert kb.formula_mass(formula) == pytest.approx(mono, abs=0.0005)


def test_nist_masses_exact():
    # NIST Atomic Weights and Isotopic Compositions, fetched 2026-11-21.
    nist = {"H": 1.00782503223, "C": 12.0, "N": 14.00307400443, "O": 15.99491461957,
            "Na": 22.9897692820, "P": 30.97376199842, "S": 31.9720711744, "Cl": 34.968852682,
            "K": 38.9637064864, "Fe": 55.93493633}
    for el, m in nist.items():
        assert kb.MONO_MASS[el] == pytest.approx(m, abs=1e-6), el


def test_fk_step_reaction_step(tables):
    _set(tables, "step_reaction", "step_id", "s1", "step_id", "sX")
    assert "step_reaction[sX,R1]: step_id not in step table" in kb.validate_kb(tables)


def test_fk_step_reaction_reaction(tables):
    _set(tables, "step_reaction", "step_id", "s1", "reaction_id", "RX")
    assert "step_reaction[s1,RX]: reaction_id not in reaction table" in kb.validate_kb(tables)


def test_fk_participant_reaction(tables):
    _set(tables, "reaction_participant", "compound", "CHEBI:A", "reaction_id", "RX")
    assert "reaction_participant[RX,CHEBI:A]: reaction_id not in reaction table" in kb.validate_kb(tables)


def test_stoich_numeric(tables):
    _set(tables, "reaction_participant", "compound", "CHEBI:A", "stoich", "two")
    assert kb.validate_kb(tables) == ["reaction_participant[R1,CHEBI:A]: stoich 'two' is not a positive number"]


def test_direction_vocab(tables):
    _set(tables, "reaction", "reaction_id", "R1", "direction", "backward")
    assert kb.validate_kb(tables) == ["reaction[R1]: direction 'backward' not in forward,reversible,unknown"]


def test_curated_vocab(tables):
    _set(tables, "reaction", "reaction_id", "R1", "curated", "yes")
    assert kb.validate_kb(tables) == ["reaction[R1]: curated 'yes' not in true,false"]


@pytest.mark.parametrize("name,row_idx,msg", [
    ("compound", 0, "compound[CHEBI:A]: duplicate key"),
    ("reaction", 0, "reaction[R1]: duplicate key"),
    ("reaction_participant", 0, "reaction_participant[R1,CHEBI:A,substrate]: duplicate key"),
    ("step", 0, "step[s1]: duplicate key"),
    ("step_reaction", 0, "step_reaction[s1,R1,Synthetica one,geneA]: duplicate key"),
])
def test_duplicate_keys(tables, name, row_idx, msg):
    df = tables[name]
    tables[name] = pd.concat([df, df.iloc[[row_idx]]], ignore_index=True)
    assert msg in kb.validate_kb(tables)


def test_pmid_digits(tables):
    _set(tables, "step_reaction", "step_id", "s1", "pmid", "PMC123")
    assert kb.validate_kb(tables) == ["step_reaction[s1,R1]: pmid 'PMC123' is not digits only"]


def test_blank_inchikey_with_block1(tables):
    _set(tables, "compound", "chebi_id", "CHEBI:A", "inchikey", "")
    assert kb.validate_kb(tables) == [
        "compound[CHEBI:A]: inchikey is empty but inchikey_block1 is 'AAAAAAAAAAAAAA'"
    ]


def test_cofactor_role_is_valid(tables):
    extra = pd.DataFrame([dict(reaction_id="R1", compound="CHEBI:D", role="cofactor", stoich="1")])
    tables["reaction_participant"] = pd.concat([tables["reaction_participant"], extra], ignore_index=True)
    assert kb.validate_kb(tables) == []


def _add_part(tables, rid, cpd, role):
    extra = pd.DataFrame([dict(reaction_id=rid, compound=cpd, role=role, stoich="1")])
    tables["reaction_participant"] = pd.concat([tables["reaction_participant"], extra], ignore_index=True)


def test_path_cofactor_product_does_not_feed(tables):
    # R2 consumes CHEBI:D as a cofactor only; R1 makes CHEBI:D as a cofactor product. Ignored.
    p = tables["reaction_participant"]
    tables["reaction_participant"] = p[~((p.reaction_id == "R2") & (p.compound == "CHEBI:B"))].reset_index(drop=True)
    _add_part(tables, "R1", "CHEBI:D", "cofactor")
    _add_part(tables, "R2", "CHEBI:D", "cofactor")
    assert kb.path_possible(tables, ["s1", "s2"]).empty


def test_path_cofactor_substrate_not_required(tables):
    _add_part(tables, "R2", "CHEBI:D", "cofactor")  # not produced by step 1, but a cofactor
    out = kb.path_possible(tables, ["s1", "s2"])
    assert out["compound"].tolist() == ["CHEBI:C"]


def test_path_requires_all_substrates(tables):
    _add_part(tables, "R2", "CHEBI:D", "substrate")  # CHEBI:D is not made by step 1
    assert kb.path_possible(tables, ["s1", "s2"]).empty


def test_path_all_substrates_met_takes_weakest(tables):
    _add_part(tables, "R1", "CHEBI:D", "product")
    _add_part(tables, "R2", "CHEBI:D", "substrate")
    out = kb.path_possible(tables, ["s1", "s2"])
    assert out["compound"].tolist() == ["CHEBI:C"]
    assert out["min_tier"].tolist() == ["T1c"]


def test_path_unknown_step_raises(tables):
    with pytest.raises(ValueError, match="unknown step_id 'sX'"):
        kb.path_possible(tables, ["s1", "sX"])
