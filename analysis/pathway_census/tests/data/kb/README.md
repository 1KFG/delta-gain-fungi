# Test fixtures for kb and kbsources

## Captured from public APIs (real data)

Captured once on 2026-11-21 from the HPCC login node with curl, User-Agent `pathway-census/0.1 (BFD-HPCC-master agent)`. Each file is the raw response (all under 30 KB). Tests never use the network.

| File | Source URL |
|---|---|
| `chebi_40968.json` | https://www.ebi.ac.uk/chebi/backend/api/public/compound/40968/ (astaxanthin) |
| `pubchem_MQZIGYBFDRPAKN.json` | https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/MQZIGYBFDRPAKN-UWFIBFSHSA-N/property/InChIKey,MolecularFormula,ExactMass/JSON |
| `uniprot_Q3HR17.json` | https://rest.uniprot.org/uniprotkb/Q3HR17.json (length 557, unreviewed) |
| `uniprot_P21334.json` | https://rest.uniprot.org/uniprotkb/P21334.json (Neurospora al-1) |
| `pubmed_16416328.xml` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=16416328&retmode=xml |

## Synthetic (invented for validator tests, not chemistry)

`compound.tsv`, `reaction.tsv`, `reaction_participant.tsv`, `step.tsv`, `step_reaction.tsv` form a small valid knowledge base. Compound ids, InChIKeys, SMILES, organisms, genes, accessions and quotes are placeholders. Formulas and masses are real values for the stated formulas. The tests copy these tables and break one rule at a time.

## Added in fix round 1 (captured 2026-11-21)

- `chebi_26078.json` (phosphoric acid), `chebi_26836.json` (sulfuric acid), `chebi_32588.json` (potassium chloride), `chebi_30672.json` (ferrocene): raw ChEBI public API responses, used as an independent check of `formula_mass`.
- `chebi_smiles_charge.json`: SMILES of L-dopa zwitterion (CHEBI:57504) and nitrobenzene (CHEBI:27798) copied from the ChEBI `es_search` API.

## Added in task 7 (captured 2026-11-21 from the cached real responses)

| File | Source |
|---|---|
| `uniprot_B2J6X9.json` | https://rest.uniprot.org/uniprotkb/B2J6X9.json (Nostoc punctiforme; ordered locus name Npun_R5600 only) |
| `ncbi_protein_WEV88679.gp` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=WEV88679&rettype=gp&retmode=text (GenBank protein, Phaffia rhodozyma DDGS) |

`tests/test_verify_kb.py` builds its small pathways in a temporary directory and writes PubMed XML inline.
