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
    _set(tables, "reaction_participant", "compound", "CHEBI:B", "role", "cofactor")
    errs = kb.validate_kb(tables)
    assert errs.count("reaction_participant[R1,CHEBI:B]: role 'cofactor' not in substrate,product") == 1


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
    assert kb.formula_mass("FeSPK") > 100


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
