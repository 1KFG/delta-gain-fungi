# Flat knowledge tables

One directory per pathway: `kb/<pathway>/` with five TSV files (UTF-8, header required). Tables exist for `carotenoid`, `melanin_dhn`, `melanin_ldopa` and `mycosporine` (stage 1, 2026-11-21); fixtures are in `tests/data/kb/`.

| File | Columns |
|---|---|
| `compound.tsv` | chebi_id, inchikey, inchikey_block1, smiles, formula, neutral_exact_mass, kegg_id, np_class |
| `reaction.tsv` | reaction_id, ec, kegg_rn, direction, curated, note, rhea_id |
| `reaction_participant.tsv` | reaction_id, compound (a chebi_id), role (substrate or product), stoich |
| `step.tsv` | step_id, route, ko, hmm_id, family, optional |
| `step_reaction.tsv` | step_id, reaction_id, organism, gene_name, protein_accession, tier, assay, pmid, quote |

Tiers: T1a, T1b, T1c, T2, T3, T4 (strongest to weakest). T1* and T2 rows need `pmid` and `quote`. Compounds are neutral parents only; the curator maps Rhea charged species to the neutral parent. `kegg_id` is informational.

Code: `pathway_census/kb.py` (`load_kb`, `validate_kb`, `formula_mass`, `path_possible`) and `pathway_census/kbsources.py` (ChEBI, PubChem, UniProt, PubMed fetchers with a JSON cache in `analysis/pathway_census/.cache/`, 3 requests per second at most).

## Vocabularies and rules

- `reaction.direction`: `forward`, `reversible`, `unknown`. It is informational. `path_possible` ignores it.
- `reaction.curated`: `true` or `false`. A `false` row needs a `rhea_id`.
- `reaction_participant.role`: `substrate`, `product`, `cofactor`. Store the roles in the direction the step uses (the stored substrates are consumed and the stored products are made by the step). `path_possible` ignores cofactor rows. A reaction of a step is applied only if ALL its non-cofactor substrates are available: products of the previous step, or compounds made earlier inside the same step. Within one step the function takes the closure over all reactions of the step, so one step can run a chain or a grid of reactions. T4 links are never used. The result is the set of compounds made by the last step, with columns `compound, min_tier, steps, organisms, single_organism, n_equal_paths`. A reaction with zero non-cofactor substrates after the first step is skipped (it can never be reached).
- `stoich` is a positive number. `pmid` is digits only when present.
- Foreign keys: `step_reaction.step_id` is in `step`; `step_reaction.reaction_id` and `reaction_participant.reaction_id` are in `reaction`; `reaction_participant.compound` is in `compound`.
- Unique keys: compound `chebi_id`; reaction `reaction_id`; step `step_id`; participant (`reaction_id`, `compound`, `role`); step_reaction (`step_id`, `reaction_id`, `organism`, `gene_name`, `protein_accession`, `pmid`). Rows that differ in any of these six fields are distinct evidence.
- Charge: `formula` may not contain `+` or `-`. SMILES may contain charged atoms (zwitterions, nitro) only when the net charge is 0.
- Fetchers raise `KbFetchError` on an HTTP error or an incomplete response, and cache only complete responses. Each cache entry carries a version tag and the request URL; a mismatch is a miss. The 3 requests per second limit is per process.

## Conventions used by the carotenoid tables

- `reaction_id` is the Rhea master id (`RHEA:<n>`) or a curated id `CUR:<pathway>:<name>` where Rhea has no reaction.
- Small molecules (O2, H2O, acceptor, NAD(P), diphosphate) are not listed as participants. A curated reaction may have no participant when ChEBI has no entry for the compound; the `note` says so.
- `step_reaction.organism` is the UniProt organism name of `protein_accession`. One paper gives one quoted sentence, used for all links from that paper. A T3 link may carry a `pmid` (the paper that was read) without a quote.
- A step that exists only at tier T4 and has no reference protein is in `step.tsv` but not in the pathway YAML (carotenoid: `torularhodin_formation`).

## Route tiers and organisms (`path_possible`)

- `min_tier` is the weakest tier along the best path. The links of that path can belong to different organisms.
- `organisms` lists the species keys of the links on the best path (sorted, comma-joined). `single_organism` is True when the list has one name.
- **When `single_organism` is False, the tier does NOT hold for any one organism.** Report such a tier only together with the organism list. Example (carotenoid, `core` route): beta-carotene is T1a with *Mucor lusitanicus* (the UniProt name is Mucor circinelloides f. lusitanicus) and *Xanthophyllomyces dendrorhous* links; no single organism reaches T1a.
- `path_possible(tables, step_ids, organism="Neurospora crassa")` uses only the links of that organism. The filter is applied before the search. In the carotenoid tables only *Neurospora crassa* completes a route alone (`core` and `neurosporaxanthin`, at T3, because its GGPP synthase link is T3).
- Best path: strongest tier; among equal tiers the fewest organisms; then the alphabetically first organism list.
- **Ties are broken alphabetically.** The organism list shown is one of `n_equal_paths` choices with the same tier and the same number of organisms. `n_equal_paths` counts organism lists, not routes through the reactions. Example (carotenoid, `core`): gamma-carotene is T1a with 3 equal lists (Mucor lusitanicus with Neurospora crassa, with Phycomyces blakesleeanus, or with Xanthophyllomyces dendrorhous); the first is shown.
- `species_key(organism)` makes the key: lowercase, text in parentheses or square brackets dropped, then the alias table `SPECIES_ALIASES` in `kb.py`, else the first two words.
- **Alias rule: a name is an alias of a key only when NCBI Taxonomy gives both names the same taxon id.** Names with different taxids are never merged. Checked on 2026-11-21 (esearch, db=taxonomy, field All Names; stored in `tests/data/kb/ncbi_taxonomy_names.tsv`; a test compares the table with it):

| Name | NCBI taxid | NCBI scientific name | Key |
|---|---|---|---|
| Gibberella fujikuroi | 5127 | Fusarium fujikuroi | `fusarium fujikuroi` |
| Fusarium fujikuroi | 5127 | Fusarium fujikuroi | `fusarium fujikuroi` |
| Phaffia rhodozyma | 264483 | Phaffia rhodozyma | `xanthophyllomyces dendrorhous` |
| Xanthophyllomyces dendrorhous | 264483 | Phaffia rhodozyma | `xanthophyllomyces dendrorhous` |
| Mucor circinelloides f. lusitanicus | 29924 | Mucor lusitanicus | `mucor lusitanicus` |
| Mucor racemosus var. lusitanicus | 29924 | Mucor lusitanicus | `mucor lusitanicus` |
| Mucor lusitanicus | 29924 | Mucor lusitanicus | `mucor lusitanicus` |
| Mucor circinelloides | 36080 | Mucor circinelloides | `mucor circinelloides` (not merged with 29924) |
| Mucor racemosus | 4841 | Mucor racemosus | `mucor racemosus` (not merged with 29924) |
| Mycosarcoma maydis | 5270 | Mycosarcoma maydis | `ustilago maydis` |
| Ustilago maydis | 5270 | Mycosarcoma maydis | `ustilago maydis` |
| Wangiella dermatitidis | 5970 | Exophiala dermatitidis | `exophiala dermatitidis` |
| Exophiala dermatitidis | 5970 | Exophiala dermatitidis | `exophiala dermatitidis` |
| Magnaporthe oryzae | 318829 | Pyricularia oryzae | `pyricularia oryzae` |
| Pyricularia oryzae | 318829 | Pyricularia oryzae | `pyricularia oryzae` |
| Magnaporthe grisea | 148305 | Pyricularia grisea | `pyricularia grisea` (not merged with 318829) |
| Pyricularia grisea | 148305 | Pyricularia grisea | `pyricularia grisea` (not merged with 318829) |
| Colletotrichum lagenaria | 5462 | Colletotrichum lagenaria | `colletotrichum lagenaria` (not merged with 5465) |
| Colletotrichum orbiculare | 5465 | Colletotrichum orbiculare | `colletotrichum orbiculare` |
| Cryptococcus neoformans | 5207 | Cryptococcus neoformans | `cryptococcus neoformans` |
| Cryptococcus neoformans var. grubii | 5207 | Cryptococcus neoformans | `cryptococcus neoformans` (first two words; no alias needed) |
| Cryptococcus deneoformans | 40410 | Cryptococcus deneoformans | `cryptococcus deneoformans` (not merged with 5207) |
| Cryptococcus neoformans var. neoformans | 40410 | Cryptococcus deneoformans | `cryptococcus deneoformans` (alias) |
| Pholiota microspora | 1538424 | Pholiota microspora | `pholiota microspora` |
| Pholiota nameko | 1538424 | Pholiota microspora | `pholiota microspora` (alias) |
| Trametes sanguinea | 158606 | Trametes sanguinea | `trametes sanguinea` |
| Nostoc punctiforme | 272131 | Nostoc punctiforme | `nostoc punctiforme` |
| Trichormus variabilis | 264691 | Trichormus variabilis | `trichormus variabilis` |
| Anabaena variabilis | 264691 | Trichormus variabilis | `trichormus variabilis` (alias) |

- Key choice: for taxid 29924 the key is the NCBI scientific name, because the old names point to two other species. For taxids 264483 and 5270 the key is the name used in the cited papers, not the current NCBI scientific name.
- Limit: a name without an alias is keyed by its first two words. A forma or variety that NCBI holds as its own taxon is kept apart only if it is in the alias table. Check the taxid before adding organisms from a new pathway.

## Conventions used by the DHN-melanin tables (`kb/melanin_dhn/`)

- Three PKS routes (`pks_thn`, `pks_hexaketide`, `pks_heptaketide`) share the steps from 1,3,6,8-tetrahydroxynaphthalene to 1,8-dihydroxynaphthalene. `scytalone_dehydratase` is listed twice in each route (scytalone, then vermelone). The oxidases are in the route `polymerisation`, which has no product compound.
- The two reductase steps each hold the main reaction of their type and the side reaction of the other type (substrate overlap, PMID 10956664, 12967213). The reference sets of the two steps are disjoint.
- Every link is in `tests/data/kb/melanin_dhn_quotes.tsv` (change detector), T3 links included.
- A reaction with no step link is allowed (`CUR:melanin_dhn:t4hn_autooxidation`, the flaviolin shunt).
- `pks_hexaketide` gives no product: 2-acetyl-1,3,6,8-tetrahydroxynaphthalene has no ChEBI entry. `refs/melanin_dhn/UNVERIFIED.md`, section D.
- Organism key of the P79068 links (`thn_synthase`): `glomerella lagenarium`, the first two words of the UniProt organism name "Glomerella lagenarium" (UniProt taxid 5462). `species_key("Colletotrichum lagenarium")`, the name in the papers, gives a third key, `colletotrichum lagenarium`. NCBI Taxonomy gives neither "Glomerella lagenarium" nor "Colletotrichum lagenarium" a taxid hit (esearch, All Names, 2026-11-21), so no alias is allowed under the taxid rule. The row `Colletotrichum lagenaria` (taxid 5462) in the table above is the NCBI scientific name of the UniProt taxid; no link has the key `colletotrichum lagenaria`. None of these is merged with Colletotrichum orbiculare (taxid 5465).
- Quotes: at most one quoted sentence per link. One paper may give different sentences to different links (PMID 10515939: one sentence for the arp1 and arp2 links, another for the abr1 link). A T3 link may carry a sentence as context (abr1).
- No organism completes a DHN-melanin route with its own links. Route tiers hold only together with the organism list.

## Conventions used by the L-DOPA / other melanin tables (`kb/melanin_ldopa/`)

- Five routes, all alternatives: `tyrosinase_dopa`, `laccase_exogenous`, `pyomelanin`, `tyrosine_degradation`, `asp_melanin`. Tyrosinase (EC 1.14.18.1) and laccase (EC 1.10.3.2) are separate steps in separate routes.
- **`exogenous_substrate`.** The schema has no field for a substrate that the organism does not make. The flag is text: each route `note` in `pathways/melanin_ldopa.yaml` starts with `exogenous_substrate: true` or `exogenous_substrate: false`, and the `lineage_note` of step `laccase_diphenol` carries `exogenous_substrate: true`. Only `laccase_exogenous` is true (Cryptococcus neoformans makes melanin only from exogenous substrates, PMID 17098915). `path_possible` does not read the flag: in the first step of a route it applies every link with no requirement on the substrates. So the products of `laccase_exogenous` mean "possible if the diphenol is supplied".
- The laccase reactions are curated (`CUR:melanin_ldopa:ldopa_oxidation_laccase`, `dopamine_oxidation_laccase`, `hga_oxidation_laccase`) with EC 1.10.3.2. Rhea has only the generic laccase reaction RHEA:11276. No laccase link points to a tyrosinase reaction (RHEA:18117, RHEA:34287), and no laccase reaction has L-tyrosine as a participant.
- Free L-DOPA is a substrate only. Rhea writes L-tyrosine to L-dopaquinone as one reaction, so no reaction of the table makes L-DOPA.
- Dopaquinone onward is non-enzymatic in the table. The step `dopachrome_onward` is in `step.tsv` only (no reference protein) and has two T4 links (`CUR:melanin_ldopa:dopaquinone_cyclisation`, RHEA:13041) with the organism text `not assigned (no enzyme identified)`. `path_possible` never uses a T4 link. Exclude T4 rows before listing organisms.
- Reactions that end in a polymer have no product participant (`hga_oxidation_laccase`, `hga_autooxidation`). `hga_autooxidation` has no step link.
- `asp_melanin` gives no product: the TyrP products have charged ChEBI entries only (CHEBI:193114, CHEBI:193115). The first step alone gives aspulvinone E. `refs/melanin_ldopa/UNVERIFIED.md`, section D.
- `tyrosine_degradation` is not a pigment route. It is declared so that `hmga` (removal of homogentisic acid) is counted against `pyomelanin`.
- Tier rule (binding rule 1, applied as written): a T1 link needs a quoted sentence that names the protein of the link. T1a also needs the assay or the product for the reaction in that sentence. A sentence that does not name the protein, or that names no substrate or product of the reaction, is held at T3 as context (the `assay` text starts with `Held at T3.`). Result: no tyrosinase link is above T1c, and the only T1a links are the two MelA links.
- Every link is in `tests/data/kb/melanin_ldopa_quotes.tsv` (change detector; columns step_id, reaction_id, protein_accession, pmid, quote, tier, assay).
- Organism keys: serotype D references have the key `cryptococcus deneoformans` (NCBI taxid 40410; alias `cryptococcus neoformans var. neoformans`), serotype A references `cryptococcus neoformans` (taxid 5207). The papers call both Cryptococcus neoformans; the keys are not merged. `Pycnoporus sanguineus` (the name in PMID 16151802) has no NCBI taxid hit, so it is not an alias of `trametes sanguinea`.
- Every route result of this pathway holds for one organism (`single_organism` True).

## Conventions used by the mycosporine tables (`kb/mycosporine/`)

- **Source of the step set.** The steps come from the fungal literature: DDGS, OMT and ATPG of one cluster in Phaffia rhodozyma UCD 67-385 (PMID 36983003). The set is not the cyanobacterial gene list. It has no second-amine step (MysD or an NRPS-like enzyme), and its ATP-grasp step makes mycosporine-glutaminol.
- **Fungal references are GenBank proteins.** UniProtKB has no entry for the three Phaffia proteins, so the accessions are GenBank protein ids (WEV88679.1, WEV88680.1, WEV88681.1; CDS of the nucleotide record OQ547787). `organism` and `species` are the GenBank organism with the strain: `Phaffia rhodozyma (strain UCD 67-385)`. `provenance.tsv` names GenBank in `source_database`. These references are `verified: false`, because the rule for `verified: true` needs a UniProt entry.
- **Bacterial references and the T2 tier.** Cyanobacterial proteins (Nostoc punctiforme, Trichormus variabilis) are references for tree placement. Their links are T2 or T3 and are never fungal evidence. A T2 link needs a quoted abstract sentence that names the gene of the protein, and either the product of that reaction in the sentence (mysC, B2J6X7: mycosporine-glycine) or a UniProt catalytic activity record with experimental evidence for the protein (mysA, B2J6X9). PMID 21890703 names NpR5600, NpR5599 and NpR5598. mysB (B2J6X8) has neither and is held at T3. A sentence that names only the enzyme class is held at T3 (Q3M6C3). A T1 link needs a fungal protein and the fungal paper; a test checks both.
- **The bacterial ATP-grasp link is on another reaction.** MysC ligates glycine. Its link points to `CUR:mycosporine:mycosporine_glycine_formation`, not to the fungal reaction `CUR:mycosporine:mycosporine_glutaminol_formation`.
- **Gene names of bacterial entries.** UniProt gives these entries an ordered locus name only (for example Npun_R5600). `kbsources.uniprot_entry` does not read ordered locus names; the curation scripts read them from the full UniProt record. The papers write the short form (NpR5600).
- **4-deoxygadusol has no ChEBI entry.** It is the product of `omt` and the substrate of `atp_grasp_ligase`. It is not in `compound.tsv`, so `ddg_o_methylation` has no product and `mycosporine_glutaminol_formation` has no substrate. Result: **no declared route gives a product** in `path_possible`. The prefix `ddgs` gives demethyl-4-deoxygadusol; the step `atp_grasp_ligase` alone (as the first step of a call) gives mycosporine-glutaminol. Other compounds with no ChEBI entry: mycosporine-glycine, mycosporine-glutaminol-glucoside, mycosporine-glutamicol-glucoside, asterina-330. `refs/mycosporine/UNVERIFIED.md`, section D.
- **Route flag.** Each route `note` starts with `fungal_evidence: true` or `fungal_evidence: partial` (text, not a schema field). `partial` means that the last step is T4.
- **T4 steps are in the YAML.** `mg_glucosyltransferase` and `glutamicol_formation` have `refs: []`, no decoy and one T4 link each (organism text `not assigned (no enzyme identified)`). They are in the YAML, unlike the T4 steps of the carotenoid and melanin tables, because the route variants need the step ids. `path_possible` never uses a T4 link.
- **Reactions with no step link:** `shinorine_formation` and `porphyra334_formation` (bacterial context, no fungal enzyme). **Compounds that no reaction makes:** gadusol, mycosporine-serinol (CHEBI:7039, ChEBI name "Mycosporine") and palythine.
- **KEGG ids of compounds.** Where neither ChEBI nor KEGG links the other record, the KEGG id is matched by KEGG name and formula (C21281, C22846).
- **Domains.** `hmm_id` of `atp_grasp_ligase` is empty: the Phaffia protein has no Pfam family at the gathering cutoff. `omt` lists two Pfam ids (any-of): PF13578 for the fungal reference and PF01596 for the bacterial one.
- The glucosylation step is T4 for the protein link. An enzyme activity without a sequence is recorded in the notes (Ascochyta fabae, PMID 24258137).
- `tests/data/kb/mycosporine_sequences.tsv` holds the length and SHA-256 of every reference and decoy sequence (written by `build.py` from the fetched records); a test compares the FASTA files with it.
- Every link is in `tests/data/kb/mycosporine_quotes.tsv` (change detector; columns step_id, reaction_id, protein_accession, pmid, quote, tier, assay).
- Reported producers are in `refs/mycosporine/fungal_lineages.tsv` (organism, group, compound, pmid, evidence, note). It is a list of reports with the PMID of each; it is not a table of the knowledge base and `load_kb` does not read it.
- Organism keys: `xanthophyllomyces dendrorhous` (Phaffia rhodozyma, taxid 264483, alias from the carotenoid tables), `nostoc punctiforme`, `trichormus variabilis` (alias `anabaena variabilis`, the name in the papers; one NCBI taxid, 264691).
