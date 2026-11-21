# Flat knowledge tables

One directory per pathway: `kb/<pathway>/` with five TSV files (UTF-8, header required). Tables exist for `carotenoid` (stage 1, 2026-11-21); fixtures are in `tests/data/kb/`.

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
- `reaction_participant.role`: `substrate`, `product`, `cofactor`. Store the roles in the direction the step uses (the stored substrates are consumed and the stored products are made by the step). `path_possible` ignores cofactor rows. A reaction of a step is applied only if ALL its non-cofactor substrates are available: products of the previous step, or compounds made earlier inside the same step. Within one step the function takes the closure over all reactions of the step, so one step can run a chain or a grid of reactions. T4 links are never used. The result is the set of compounds made by the last step, with columns `compound, min_tier, steps, organisms, single_organism`. A reaction with zero non-cofactor substrates after the first step is skipped (it can never be reached).
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
- **When `single_organism` is False, the tier does NOT hold for any one organism.** Report such a tier only together with the organism list. Example (carotenoid, `core` route): beta-carotene is T1a with *Mucor circinelloides* and *Xanthophyllomyces dendrorhous* links; no single organism reaches T1a.
- `path_possible(tables, step_ids, organism="Neurospora crassa")` uses only the links of that organism. The filter is applied before the search. In the carotenoid tables only *Neurospora crassa* completes a route alone (`core` and `neurosporaxanthin`, at T3, because its GGPP synthase link is T3).
- Best path: strongest tier; among equal tiers the fewest organisms; then the alphabetically first organism list.
- `species_key(organism)` makes the key: lowercase, text in parentheses or square brackets dropped, first two words, after the alias table `SPECIES_ALIASES` in `kb.py`. Each alias names its source (NCBI Taxonomy record, fetched 2026-11-21): Gibberella fujikuroi = Fusarium fujikuroi (taxid 5127); Phaffia rhodozyma = Xanthophyllomyces dendrorhous (taxid 264483); Mucor circinelloides f. lusitanicus = Mucor racemosus var. lusitanicus = Mucor lusitanicus (taxid 29924, key `mucor circinelloides`); Mycosarcoma maydis = Ustilago maydis (taxid 5270). The key is the name used in the cited papers; it is not always the current NCBI scientific name (NCBI uses Phaffia rhodozyma, Mucor lusitanicus and Mycosarcoma maydis).
