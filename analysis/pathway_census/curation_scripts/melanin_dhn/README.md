# One-off build helpers for the DHN-melanin curation (stage 1, task S4)

These scripts are one-off build helpers. They fetch records from public APIs (UniProt, PubMed, ChEBI, PubChem, Rhea, KEGG `list`, NCBI Taxonomy) through `pathway_census.kbsources` (rate limit, user agent, cache in `analysis/pathway_census/.cache/`). They are not run by the tests and they are not tested code. Task 7 (S7) replaces them with a tested verifier (`bin/verify_kb.py`).

Run from any directory with `/usr/bin/python3.12 analysis/pathway_census/curation_scripts/melanin_dhn/<script>.py`. A run needs network access unless the cache is full. `a5.py` needs Biopython.

| Script | What it does |
|---|---|
| `h.py` | Shared helper (copy of the carotenoid helper). `raw()` fetches a URL that has no `kbsources` fetcher, through `kbsources._cached` and `default_fetch`. |
| `build.py` | Holds the curated data (references, decoys, tiers, quote positions, reactions, routes). It checks each fact against the fetched records and writes `pathways/melanin_dhn.yaml`, `refs/melanin_dhn/` (FASTA, provenance, decoys), `kb/melanin_dhn/*.tsv` and `tests/data/kb/melanin_dhn_quotes.tsv`. It does not write `UNVERIFIED.md`. |
| `verify.py` | Second pass from the written files: organism, gene name, sequence, PMID, quote in abstract, compound records, route table. Prints Markdown tables. Exit code 1 on a mismatch. |
| `a1.py` | Prints the checked fields of UniProt entries (arguments: accessions). |
| `a2.py` | UniProt searches (arguments: queries). |
| `a3.py` | PubMed abstracts, paced at 2 requests per second or less (arguments: PMIDs). |
| `a4.py` | ChEBI name search, ChEBI record, Rhea search. |
| `a5.py` | Pairwise global identity between references and decoys (arguments: accessions). Source of the MEASURED identities in `UNVERIFIED.md`. |
| `b1.py` | NCBI Taxonomy taxid of the organism names. With `--write` it appends new rows to `tests/data/kb/ncbi_taxonomy_names.tsv` (the source of `kb.SPECIES_ALIASES`). |

Quotes: `build.py` cuts each quoted sentence out of the fetched abstract by a start and an end fragment, so the quote text is never typed by hand. One paper gives one sentence.
