# Flat knowledge tables

One directory per pathway: `kb/<pathway>/` with five TSV files (UTF-8, header required). No pathway tables exist yet; fixtures are in `tests/data/kb/`.

| File | Columns |
|---|---|
| `compound.tsv` | chebi_id, inchikey, inchikey_block1, smiles, formula, neutral_exact_mass, kegg_id, np_class |
| `reaction.tsv` | reaction_id, ec, kegg_rn, direction, curated, note, rhea_id |
| `reaction_participant.tsv` | reaction_id, compound (a chebi_id), role (substrate or product), stoich |
| `step.tsv` | step_id, route, ko, hmm_id, family, optional |
| `step_reaction.tsv` | step_id, reaction_id, organism, gene_name, protein_accession, tier, assay, pmid, quote |

Tiers: T1a, T1b, T1c, T2, T3, T4 (strongest to weakest). T1* and T2 rows need `pmid` and `quote`. Compounds are neutral parents only; the curator maps Rhea charged species to the neutral parent. `kegg_id` is informational.

Code: `pathway_census/kb.py` (`load_kb`, `validate_kb`, `formula_mass`, `path_possible`) and `pathway_census/kbsources.py` (ChEBI, PubChem, UniProt, PubMed fetchers with a JSON cache in `analysis/pathway_census/.cache/`, 3 requests per second at most).
