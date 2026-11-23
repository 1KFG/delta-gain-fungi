# One-off build helpers for the carotenoid curation (stage 1, task S3)

These scripts are one-off build helpers. They fetch records from public APIs (UniProt, PubMed, ChEBI, PubChem, Rhea, InterPro, NCBI protein) through `pathway_census.kbsources` (rate limit, user agent, cache in `analysis/pathway_census/.cache/`). They are not run by the tests and they are not tested code. Task 7 (S7) replaces them with a tested verifier (`bin/verify_kb.py`).

Run from any directory with `/usr/bin/python3.12 analysis/pathway_census/curation_scripts/carotenoid/<script>.py`. A run needs network access unless the cache is full.

| Script | What it does |
|---|---|
| `h.py` | Shared helper. `raw()` fetches a URL that has no `kbsources` fetcher, through `kbsources._cached` and `default_fetch`. |
| `build.py` | Holds the curated data (references, decoys, tiers, quote positions, reactions, routes). It checks each fact against the fetched records and writes `pathways/carotenoid.yaml`, `refs/carotenoid/` (FASTA, provenance, decoys), `kb/carotenoid/*.tsv` and `tests/data/kb/carotenoid_quotes.tsv`. It does not write `UNVERIFIED.md`. |
| `verify.py` | Second pass from the written files: organism, gene name, sequence, PMID, quote in abstract, compound records. |
| `rep.py` | Writes the tables of the task report (`tables.md` in the current directory). |
| `a1.py` to `a11.py`, `b1.py`, `b2.py` | Exploration scripts, kept as a record of the queries: UniProt entries and searches, PubMed abstracts, domain coordinates, Rhea and ChEBI searches, NCBI protein records, NCBI Taxonomy names (`b2.py`). |
| `b3.py` | Looks up the NCBI taxid of every organism name used by `kb.species_key` and writes `tests/data/kb/ncbi_taxonomy_names.tsv` (the source of `kb.SPECIES_ALIASES`). `a3.py` takes an output file and PMIDs as arguments. |

Quotes: `build.py` cuts each quoted sentence out of the fetched abstract by a start and an end fragment, so the quote text is never typed by hand. One paper gives one sentence.
