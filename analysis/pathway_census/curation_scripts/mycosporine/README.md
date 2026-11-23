# One-off build helpers for the mycosporine curation (stage 1, task S6)

These scripts are one-off build helpers. They fetch records from public APIs (UniProt, UniParc, NCBI protein, PubMed, PMC, Europe PMC, ChEBI, PubChem, Rhea, KEGG, NCBI Taxonomy) through `pathway_census.kbsources` (rate limit, user agent, cache in `analysis/pathway_census/.cache/`). They are not run by the tests and they are not tested code. Task 7 (S7) replaces them with a tested verifier (`bin/verify_kb.py`).

Run from the script directory (the scripts import `helpers.py` from it): `cd analysis/pathway_census/curation_scripts/mycosporine`, then `/usr/bin/python3.12 <script>.py`. A run needs network access unless the cache is full. `pairwise_identity.py` needs Biopython. All paths are relative to the script file.

| Script | What it does |
|---|---|
| `helpers.py` | Shared helper. `raw()` fetches a URL that has no `kbsources` fetcher, through `kbsources._cached` and `default_fetch`. |
| `genbank_proteins.py` | Reads the three GenBank proteins of the Phaffia rhodozyma UCD 67-385 cluster OQ547787 (WEV88679.1 DDGS, WEV88680.1 OMT, WEV88681.1 ATPG). UniProtKB has no entry for them. `--fasta FILE` writes a FASTA file. |
| `records.py` | One interface, `rec(acc)`, for UniProtKB entries and the GenBank proteins: organism, gene names, sequence, EC, Rhea ids, Pfam / NCBIfam / PROSITE ids, cited PMIDs. It adds the ordered locus names, which `kbsources.uniprot_entry` does not read. |
| `build.py` | Holds the curated data (references, decoys, tiers, quote positions, reactions, routes, lineage rows). It checks each fact against the fetched records with an assert and writes `pathways/mycosporine.yaml`, `refs/mycosporine/` (FASTA, provenance, decoys, `fungal_lineages.tsv`), `kb/mycosporine/*.tsv`, `tests/data/kb/mycosporine_quotes.tsv` and `tests/data/kb/mycosporine_sequences.tsv` (length and SHA-256 of each fetched sequence). It does not write `UNVERIFIED.md`. |
| `verify.py` | Second pass from the written files: organism, gene name, sequence, PMID, quote in abstract, whether the record cites the paper of each link, decoy EC and Rhea ids, compound records, route table, lineage PMIDs. Prints Markdown tables. Exit code 1 on a mismatch. |
| `show_uniprot.py` | Prints the checked fields of UniProt entries, with PROSITE ids and ordered locus names (arguments: accessions). |
| `search_uniprot.py` | UniProt searches (arguments: queries). |
| `show_pubmed.py` | PubMed abstracts, paced at 2 requests per second or less (arguments: PMIDs). |
| `search_pubmed.py` | PubMed search (esearch, then esummary titles), paced at 2 requests per second or less (arguments: queries). |
| `fetch_fulltext.py` | Open-access full text from Europe PMC, with a word search (arguments: PMID, words). Prints `NO FULL TEXT` when Europe PMC has none or the fetch fails. |
| `search_chebi_rhea.py` | ChEBI name search, ChEBI record, Rhea search. |
| `find_compounds.py` | For a compound name: PubChem compounds of that name, then a ChEBI search by InChIKey and by its first block. A name with no ChEBI hit is a USER DECISION. |
| `pairwise_identity.py` | Pairwise global identity between references and decoys (arguments: accessions; UniProt or the WEV proteins). Source of the MEASURED identities in the lineage notes and in `UNVERIFIED.md`. |
| `ncbi_taxonomy.py` | NCBI Taxonomy taxid of the organism names. With `--write` it appends new rows to `tests/data/kb/ncbi_taxonomy_names.tsv`. |

Quotes: `build.py` cuts each quoted sentence out of the fetched abstract by a start and an end fragment, so the quote text is never typed by hand. A link has at most one sentence.

Full text: `build.py` asserts six sentences of the full text of PMID 36983003 (Europe PMC, PMC10057719) and three sentences of the full texts of PMID 22741921 and 20813918 (NCBI PMC efetch, PMC3405187 and PMC3116657). Full text is used for notes and ties only. A quote in `step_reaction.tsv` is always a sentence of the abstract.

Pfam domains of the Phaffia proteins: the proteins have no UniProt cross-references, so the domains were computed once on 2026-11-21 on the UCR HPCC and the result is stored in `build.py` (`HMMSCAN_GA`):

```
/usr/bin/python3.12 genbank_proteins.py --fasta phaffia.faa
module load hmmer/3
hmmscan --cut_ga --noali --domtblout phaffia.domtbl <Pfam 38.2>/Pfam-A.hmm phaffia.faa
```

Result (gathering cutoff): WEV88679.1 PF01761, PF24621, PF13685; WEV88680.1 PF13578; WEV88681.1 no family. Without the cutoff (`-E 1`): WEV88680.1 also PF01596 (sequence E-value 1.4e-07); WEV88681.1 PF02786 (0.054) and PF02655 (0.13). The Pfam database was `/srv/projects/db/pfam/2026-01-27-Pfam38.2` on that cluster; the path is not used by any script.

KEGG: the compound entries C21281, C22846, C22843 and C10607 and the KO search for K21342 are read through the cache by `build.py`.
