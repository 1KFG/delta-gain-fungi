# One-off build helpers for the L-DOPA / other melanin curation (stage 1, task S5)

These scripts are one-off build helpers. They fetch records from public APIs (UniProt, PubMed, ChEBI, PubChem, Rhea, NCBI Taxonomy) through `pathway_census.kbsources` (rate limit, user agent, cache in `analysis/pathway_census/.cache/`). They are not run by the tests and they are not tested code. Task 7 (S7) replaces them with a tested verifier (`bin/verify_kb.py`).

Run from the script directory (the scripts import `h.py` from it): `cd analysis/pathway_census/curation_scripts/melanin_ldopa`, then `/usr/bin/python3.12 <script>.py`. A run needs network access unless the cache is full. `a5.py` needs Biopython. All paths are relative to the script file.

| Script | What it does |
|---|---|
| `h.py` | Shared helper (copy of the DHN-melanin helper). `raw()` fetches a URL that has no `kbsources` fetcher, through `kbsources._cached` and `default_fetch`. |
| `build.py` | Holds the curated data (references, decoys, tiers, quote positions, reactions, routes). It checks each fact against the fetched records with an assert and writes `pathways/melanin_ldopa.yaml`, `refs/melanin_ldopa/` (FASTA, provenance, decoys), `kb/melanin_ldopa/*.tsv` and `tests/data/kb/melanin_ldopa_quotes.tsv`. It does not write `UNVERIFIED.md`. |
| `verify.py` | Second pass from the written files: organism, gene name, sequence, PMID, quote in abstract, whether the UniProt entry cites the paper of each link, compound records, route table. Prints Markdown tables. Exit code 1 on a mismatch. |
| `a1.py` | Prints the checked fields of UniProt entries (arguments: accessions). |
| `a2.py` | UniProt searches (arguments: queries). |
| `a3.py` | PubMed abstracts, paced at 2 requests per second or less (arguments: PMIDs). |
| `a4.py` | ChEBI name search, ChEBI record, Rhea search. |
| `a5.py` | Pairwise global identity between references and decoys (arguments: accessions). Source of the MEASURED identities in the lineage notes and in `UNVERIFIED.md`. |
| `a6.py` | PubMed search (esearch, then esummary titles), paced at 2 requests per second or less (arguments: queries). |
| `b1.py` | NCBI Taxonomy taxid of the organism names. With `--write` it appends new rows to `tests/data/kb/ncbi_taxonomy_names.tsv` (the source of `kb.SPECIES_ALIASES`). |

Quotes: `build.py` cuts each quoted sentence out of the fetched abstract by a start and an end fragment, so the quote text is never typed by hand. A link has at most one sentence; a paper can give a second sentence to another link (key `<pmid>:<name>` in `Q`).

KEGG: the five KO names (K00505, K00421, K05909, K00457, K00451) were read once with `curl https://rest.kegg.jp/list/ko:<id>` on 2026-11-21. No script stores KEGG content.
