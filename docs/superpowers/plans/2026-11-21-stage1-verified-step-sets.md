# Stage 1: Verified Step Sets, References, Decoys and Flat Tables Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce verified, route-aware step sets for carotenoid, DHN-melanin and L-DOPA/other melanin pathways, with per-domain references, decoys, and flat knowledge tables (compounds, reactions, steps, evidence), and hand the user a review packet before any real-data run.

**Architecture:** Extend the pathway definition schema (routes, any-of domains, organism-keyed references, nullable EC). Add a flat-table module with validators and cached public-database fetchers (ChEBI, PubChem, UniProt, PubMed). Curate each pathway into YAML plus TSV tables. Verify every accession and PMID against its source. End with a review packet.

**Tech Stack:** Python 3.12 (`/usr/bin/python3.12`), pandas, pyyaml, pytest, the `pathway_census` package on `main`, UniProt/ChEBI/PubChem/PubMed/Rhea REST (read-only, cached).

**Spec:** `docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-design.md` (version 2, approved 2026-11-21). Read it first. Source material for curation, to be treated as leads and never as proof: `docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-REVIEW_PATHWAYS.md` (REVIEWED corrections, PMIDs, tiers) and `...-REVIEW_FABLE.md`.

**Out of scope here:** the real-data census, tree placement, SIRIUS, the EBC comparison, and the species tree. The earlier `2026-11-21-pathway-census-carotenoid-real-data.md` plan stays valid for stage 2 once the step sets here are approved; its R1 is superseded by this plan.

## Global Constraints

- Run all Python with `/usr/bin/python3.12`. Never the default `python3`.
- Public-database access is read-only, rate-limited (KEGG at most 3 requests/s; others polite), and cached under `analysis/pathway_census/.cache/` (git-ignored). Tests use stored fixtures and never hit the network.
- **Never invent an accession, PMID, EC number, ChEBI id or InChIKey.** If a fact cannot be checked, leave the field null and add it to the task report's "unverified" list.
- `verified: true` on a reference means: accession fetched from UniProt (or NCBI) on the run date, organism and gene name match, the entry or the abstract supports the step's reaction, and the PMID exists. It never means a person reviewed it.
- Gene names are always keyed by organism (`carS` means two different genes in *Phycomyces* and *Fusarium fujikuroi*).
- A tier belongs to the named organism and protein only. Tiers: T1a (fungal enzyme assayed or product identified by heterologous expression), T1b (deletion or complementation with intermediate analysis in the fungus), T1c (deletion, colour phenotype only), T2 (shown in another organism), T3 (inferred), T4 (hypothesis). T1 and T2 links need a PMID and a quoted sentence.
- Compounds are neutral parents. Rhea charged species map to the neutral compound. Match structures on the first block of the InChIKey. Never join through KEGG's PubChem field (it is a substance ID). KEGG ids are cross-reference columns only; copy no KEGG pathway or reaction content.
- Do not copy publisher text beyond one quoted sentence per link.
- Label numbers MEASURED, REVIEWED, REPORTED, PROPOSED or UNVERIFIED in docs.
- Do not commit data caches, PDFs or large files. Keep curated TSV and YAML small.
- All work in `analysis/pathway_census/`; docs in `docs/superpowers/`. The worktree is `DeltaGain_Fungi_wt_realdata` (branch `pathway-census-realdata`).
- This is computational annotation and database curation. No wet-lab work.

## Review Focus

1. **A reference verified against the wrong protein** (accession belongs to another gene or species). Task S6 re-fetches every accession and compares organism and gene name.
2. **A tier that travels between orthologues.** Validators reject a tier without organism and gene name.
3. **A charged or non-neutral species in the compound table.** Validators check formula, mass and charge.
4. **A route variant collapsed into one linear list** (DHN-melanin). Schema tests require routes.
5. **A curated quote that is not in the cited abstract.** S6 checks quote text against the fetched PubMed abstract.

---

### Task S1: Schema v2 for pathway definitions

**Files:**
- Modify: `analysis/pathway_census/pathway_census/defs.py`, `analysis/pathway_census/pathway_census/quality.py` (only the any-of handling), `analysis/pathway_census/pathways/carotenoid.yaml` (rename keys only)
- Test: `tests/test_defs.py`, `tests/test_quality.py`, fixtures under `tests/data/`

**Interfaces:**
- Produces:
  - `Ref(accession: str | None, species: str, pmid: str | None, verified: bool, gene_name: str | None, region: tuple[int, int] | None)`; `region` is 1-based inclusive and must lie inside the sequence when given.
  - `Step(..., domains_any: list[str], ec: str | None, route: str | None, reactions: list[str], decoys: list[Ref], lineage_note: str | None)`; `copy_expect` and `role` stay required.
  - `Route(id: str, name: str, step_ids: list[str], note: str | None)`; `Pathway(id, steps, verified, routes: list[Route])`.
  - `load_pathway(path)` rejects the old key `domains` with the message "use domains_any"; `validate(p)` reports: a route naming an unknown step, a step with a `route` not declared, a reaction id duplicated within a step, a decoy accession that is also a reference of the same step.
  - `quality.quality_flags` treats a step's expected domains as **any-of**: `True` if at least one listed domain is present, `False` if none, `NA` if the list is empty.

- [ ] **Step 1: Write failing tests** with exact assertions: loading a v2 fixture with two routes; the old `domains` key raises with that message; route naming an unknown step is a validate error; a decoy accession equal to a reference accession is a validate error; `region` outside the sequence length is rejected (give a sequence length via the reference FASTA when present, otherwise skip the check); any-of domain flag True/False/NA cases.
- [ ] **Step 2: Run the tests; expect FAIL.**
- [ ] **Step 3: Implement.** Keep `copy_expect`/`role` required. Keep v1 loader behaviour for every other key. Rename `domains` to `domains_any` in `carotenoid.yaml` and in existing fixtures and tests.
- [ ] **Step 4: Run the full suite** `/usr/bin/python3.12 -m pytest analysis/pathway_census/tests -q`; expect all pass.
- [ ] **Step 5: Commit** `feat(pathway_census): pathway schema v2 (routes, any-of domains, decoys)`.

---

### Task S2: Flat knowledge tables and public-database fetchers

**Files:**
- Create: `analysis/pathway_census/pathway_census/kb.py`, `analysis/pathway_census/pathway_census/kbsources.py`, `analysis/pathway_census/kb/README.md`
- Test: `tests/test_kb.py`, `tests/test_kbsources.py`, fixtures `tests/data/kb/`

**Interfaces:**
- Produces (all under `kb.py`):
  - `TIERS = ("T1a","T1b","T1c","T2","T3","T4")`.
  - `load_kb(kb_dir: Path) -> dict[str, pd.DataFrame]` reading `compound.tsv`, `reaction.tsv`, `reaction_participant.tsv`, `step.tsv`, `step_reaction.tsv` with the columns in spec section 6.1.
  - `validate_kb(tables: dict[str, pd.DataFrame]) -> list[str]` (empty when valid), rules:
    - `tier` in `TIERS`; `step_reaction` needs non-empty `organism` and `gene_name`; T1* and T2 need `pmid` and `quote`.
    - every `reaction_participant.compound` exists in `compound`; roles are `substrate` or `product`.
    - `inchikey_block1` equals the first 14 characters of `inchikey`; `formula` has no charge sign; `neutral_exact_mass` agrees with the monoisotopic mass computed from `formula` within 0.0005.
    - `kegg_id` present only as an informational column; no column named like a PubChem field from KEGG.
    - a `reaction` marked `curated=false` must have a `rhea_id`.
  - `formula_mass(formula: str) -> float` (elements C, H, N, O, P, S, Na, K, Cl, Fe; raises on others).
  - `path_possible(tables, step_ids: list[str]) -> pd.DataFrame` returning per compound reachable by consecutive `step_reaction` links with columns `compound, min_tier, steps`; tier order T1a < T1b < T1c < T2 < T3 < T4 (lower is stronger); a path containing T4 is excluded.
- Produces (under `kbsources.py`), each with an on-disk JSON cache and a `fetch=` hook so tests inject stored responses:
  - `chebi_compound(chebi_id: str) -> dict` (name, formula, monoisotopic mass, InChIKey, SMILES; read from the public ChEBI API).
  - `pubchem_by_inchikey(inchikey: str) -> dict` (CID, InChIKey, formula, exact mass).
  - `uniprot_entry(accession: str) -> dict` (gene names, organism, length, sequence, EC list, cited PMIDs, reviewed flag).
  - `pubmed_abstract(pmid: str) -> dict` (title, abstract text, year).
  - Politeness: at most 3 requests per second; a descriptive User-Agent that carries no personal information.

- [ ] **Step 1: Write failing tests** with exact values: a valid fixture passes; each rule above has one failing fixture with the exact error text; `formula_mass("C40H52O4")` equals 596.3866 within 0.0005 and `formula_mass("C10H8O4")` equals 192.0423; `path_possible` on a three-step fixture returns the weakest tier and excludes a path with a T4 step; fetchers parse stored response fixtures (ChEBI 40968 gives InChIKey `MQZIGYBFDRPAKN-UWFIBFSHSA-N`; UniProt Q3HR17 gives length 557); a cache hit makes no second call.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** Use stored fixtures captured once from the real APIs (record the capture date in a README under `tests/data/kb/`).
- [ ] **Step 4: Run the full suite; expect all pass.**
- [ ] **Step 5: Commit** `feat(pathway_census): flat knowledge tables, validators and cached database fetchers`.

---

### Task S3: Carotenoid curation (model: curation tier)

**Files:**
- Modify: `analysis/pathway_census/pathways/carotenoid.yaml`; `refs/carotenoid/*.faa`; `refs/carotenoid/provenance.tsv`; `refs/carotenoid/UNVERIFIED.md`
- Create: `analysis/pathway_census/refs/carotenoid/decoys/<step>.faa`, `refs/carotenoid/decoys/provenance.tsv`; `analysis/pathway_census/kb/carotenoid/` TSVs (the five tables, carotenoid rows only)
- Test: `tests/test_carotenoid_def.py` (update), `tests/test_kb_carotenoid.py`

**Requirements (from the spec's section 6.3, REVIEWED corrections):**
1. Phytoene desaturase: two reactions (four-step lycopene-forming; five-step 3,4-didehydrolycopene-forming); `ec` null on the step; the five-step reaction T1a only for *Neurospora* AL-1 (PMID 11017770), the four-step reaction T1a only for *Phycomyces* carB (PMID 9079885). Verify both PMIDs and both UniProt accessions (P21334, P54982) before use.
2. Astaxanthin synthase: reference accession `Q3HR17` (fetch it; expect 557 residues and identity with CAC16261); keep `ec` null; CrtR (PMID 18837978) as a separate step; model the astaxanthin route as a reaction grid (ketolation first or hydroxylation first, both to phoenicoxanthin) with tiers T1a for 4-ketolation and T1b for 3-hydroxylation; torulene to 4-ketotorulene to 3-hydroxy-4-ketotorulene by the same enzyme.
3. Fused phytoene synthase and lycopene cyclase references split **per domain**: cut each reference sequence by UniProt/Pfam/InterPro domain coordinates (PF00494 phytoene synthase; cyclase domain by PF18916 or TIGR03462 where annotated). Record `region`. A reference without usable coordinates is dropped from the split files and listed as unverified. The two split FASTA files must differ.
4. Add steps with evidence from the pathway review (verify each): *Phycomyces* carS and acaA (PMID 21854466, 22445959); tsp3 (PMID 17822929); TSP2 (PMID 18931040); tdh (PMID 8828234); *Ustilago* Cco1 as a reference of the carX step (PMID 19584000); keep carX, carT/cao-2, ylo-1/carD with their existing verified references.
5. Mark as T4 or exclude: torulene to torularhodin (no enzyme found), carotenoid glycosylation, later trisporoid oxidations.
6. Decoys per step (3-8, UniProt, verified, with reason): farnesyl/geranyl diphosphate synthases for GGPP synthase; squalene synthase for phytoene synthase; other FAD amine oxidases for phytoene desaturase; non-carotenoid carotenoid-cleavage-family members including *Neurospora* CAO-1 (PMID 23893079) for carX/carT; other aldehyde dehydrogenases for ylo-1; fungal cytochrome P450s that are not astaxanthin synthase for CrtS. If a step has no suitable decoy, say so.
7. Compounds and reactions into the flat tables: every compound with ChEBI id, InChIKey (from ChEBI, confirmed with PubChem), formula and neutral mass (mass computed, not typed); Rhea ids where Rhea has the reaction and a curated id where it does not. Never use RHEA:56772 for the *X. dendrorhous* route.
8. Gene names carry the organism (use `gene_name` plus species).

- [ ] **Step 1: Write failing tests**: `validate(load_pathway(...)) == []`; `validate_kb` passes on the carotenoid tables; phytoene desaturase step has `ec is None` and two reactions; the astaxanthin step has accession Q3HR17 and CrtR is a separate step; the two split FASTA files differ and every `region` lies inside its sequence; each decoy accession is in `decoys/provenance.tsv` and is not a reference of the same step; every T1/T2 link has a PMID and quote.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Curate** (network use via the `kbsources` fetchers, cached). Write what could not be verified into `UNVERIFIED.md`. The pathway stays `verified: false` until the user's review.
- [ ] **Step 4: Run the full suite; expect all pass.**
- [ ] **Step 5: Commit** `feat(pathway_census): carotenoid step set v2 with per-domain references, decoys and flat tables`.

---

### Task S4: DHN-melanin curation (model: curation tier)

**Files:**
- Create: `analysis/pathway_census/pathways/melanin_dhn.yaml`, `analysis/pathway_census/refs/melanin_dhn/*.faa`, `refs/melanin_dhn/provenance.tsv`, `refs/melanin_dhn/decoys/`, `refs/melanin_dhn/UNVERIFIED.md`, `analysis/pathway_census/kb/melanin_dhn/` TSVs
- Test: `tests/test_melanin_dhn_def.py`, `tests/test_kb_melanin_dhn.py`

**Requirements (spec section 6.3, REVIEWED):**
1. Route variants as `routes`: (a) THN synthase, *Colletotrichum lagenarium* PKS1 (PMID 10501004, 23261597); (b) hexaketide PKS plus Yg1, *Wangiella dermatitidis* WdPKS1 and WdYg1 (PMID 11179356, 18676950); (c) heptaketide PKS plus Ayg1, *Aspergillus fumigatus* Alb1 and Ayg1 (PMID 11350964). Ayg1/Yg1 is required with (b) and (c) only.
2. Two reductase types as separate steps and reference clades: 4HNR/Arp2 type (T4HN to scytalone, EC 1.1.1.252; PMID 10956664) and 3HNR/THR1/BUF1/BRM2 type (T3HN to vermelone; PMID 8672814, 10956664). `lineage_note` says one paralogue present does not mean the other is absent.
3. Scytalone dehydratase (SCD1, Arp1, SDH1; PMID 8953707, 9383199, 10515939). The vermelone dehydration is a **curated reaction with its own evidence record** and tier T3 unless a paper shows it; record the uncertainty about which enzyme acts.
4. Oxidases: laccase (Abr2 type) and multicopper oxidase (Abr1 type) as two families, tier T1c (PMID 10515939, 16988817, 21277986).
5. Lineage restrictions in `lineage_note`: *A. terreus* lacks DHN-melanin and has Asp-melanin; Botrytis and other Leotiomycetes carry two PKS (PMID 26514268, 27133313).
6. Decoys: *A. nidulans* wA (UniProt Q03149) for the PKS clades; other short-chain dehydrogenases for the reductases; other multicopper oxidases and laccases for the oxidase step. Verify each.
7. Add compounds: scytalone, 1,3,6,8-tetrahydroxynaphthalene, 1,3,8-trihydroxynaphthalene, vermelone, 1,8-dihydroxynaphthalene, flaviolin, YWA1/heptaketide intermediates where a ChEBI entry exists.
8. KEGG cross-references only: C00779, C04033, C01624, C18012, and KOs K17740, K17739; do not copy KEGG reaction content. K19580 `rppA` is a bacterial enzyme: do not use it for fungi.

- [ ] **Step 1: Failing tests**: schema validates; three routes exist and the Ayg1/Yg1 step is only in routes (b) and (c); two reductase steps with distinct reference sets; vermelone dehydration is a separate curated reaction with tier T3 or a verified paper; every reference accession is in provenance; flat tables validate; decoys do not overlap references.
- [ ] **Step 2: Run; expect FAIL. Step 3: Curate. Step 4: Run the full suite; expect pass.**
- [ ] **Step 5: Commit** `feat(pathway_census): DHN-melanin step set with route variants`.

---

### Task S5: L-DOPA melanin and other melanin routes (model: curation tier)

**Files:**
- Create: `analysis/pathway_census/pathways/melanin_ldopa.yaml`, `refs/melanin_ldopa/`, `kb/melanin_ldopa/`
- Test: `tests/test_melanin_ldopa_def.py`, `tests/test_kb_melanin_ldopa.py`

**Requirements (spec section 6.3, REVIEWED):**
1. Tyrosinase (EC 1.14.18.1; *Neurospora crassa* P00440; *Aspergillus oryzae* melB; *Agaricus bisporus*) and laccase (EC 1.10.3.2; *Cryptococcus neoformans* LAC1, PMID 8300520) are **separate steps**. Laccase carries diphenol oxidation only, not tyrosine hydroxylation.
2. The *Cryptococcus* route is substrate-dependent: the fungus does not make L-DOPA (PMID 17098915). Record this in `lineage_note` and as a route variant flag `exogenous_substrate`.
3. Other melanin routes as small branches or explicit exclusions with references: pyomelanin (HppD, HmgA; PMID 19028908) and Asp-melanin (MelA NRPS-like, TyrP tyrosinase; PMID 27133313; UniProt Q0CRX1, Q0CRX0; RHEA:63824). Dopaquinone to dopachrome onward is non-enzymatic in the table, T4 for any fungal enzyme.
4. Decoys for the tyrosinase and laccase steps (other type-3 copper proteins, other multicopper oxidases). Verify each.
5. Compounds: L-tyrosine, L-DOPA, dopaquinone, L-dopachrome, homogentisate, aspulvinone E where ChEBI has them; neutral parents only.

- [ ] **Step 1: Failing tests**: tyrosinase and laccase are distinct steps; the laccase step has no reaction whose substrate is tyrosine; `exogenous_substrate` is set on the *Cryptococcus* route; flat tables validate; decoys do not overlap references.
- [ ] **Step 2: Run; expect FAIL. Step 3: Curate. Step 4: Run the full suite; expect pass.**
- [ ] **Step 5: Commit** `feat(pathway_census): L-DOPA melanin and other melanin routes`.

---

### Task S6: Independent reference verification

**Files:**
- Create: `analysis/pathway_census/bin/verify_kb.py`, `analysis/pathway_census/kb/VERIFICATION_REPORT.md` (generated, committed)
- Test: `tests/test_verify_kb.py`, fixtures `tests/data/kb/`

**Interfaces:**
- Consumes: `kbsources` fetchers, all three pathway YAMLs and their `kb/` tables.
- Produces: `verify_pathway(yaml_path: Path, kb_dir: Path, *, fetch=None) -> pd.DataFrame` with columns `check, subject, expected, found, status` where status is `pass`, `fail` or `unchecked`. Checks: every reference accession exists in UniProt and its organism and gene name match the YAML; each `pmid` exists in PubMed; each `quote` occurs (after whitespace normalisation) in the fetched abstract, or status `unchecked` with the reason when the abstract is not available; each ChEBI id gives the formula and InChIKey in `compound.tsv`; the first InChIKey block agrees between ChEBI and PubChem; each decoy accession exists.

- [ ] **Step 1: Failing tests** with stored fixtures: a matching reference passes; a swapped organism fails; a missing PMID fails; a quote not in the abstract fails; an unavailable abstract is `unchecked`.
- [ ] **Step 2: Run; expect FAIL. Step 3: Implement. Step 4: Run the tests; expect pass.**
- [ ] **Step 5: Run `verify_kb.py` over the three pathways** (network, cached, rate-limited). Write the table summary to `kb/VERIFICATION_REPORT.md`: counts of pass, fail and unchecked per pathway and every fail with its detail. Fix curation errors in the YAML and tables found by failures (separate commit) or list them as open. Never flip a `fail` to `pass` without a source.
- [ ] **Step 6: Commit** `feat(pathway_census): independent knowledge-base verification and report`.

---

### Task S7: Review packet for the user

**Files:**
- Create: `docs/superpowers/plans/2026-11-21-stage1-review-packet.md`
- Test: `tests/test_review_packet.py`

**Interfaces:**
- Produces: `bin/make_review_packet.py` writing a markdown packet from the YAMLs and tables: per pathway, a table of steps (route, enzyme, organism, gene, accession, tier, PMID), the decoys, the compounds with neutral mass, the T4 and excluded items, the unverified list, and the verification counts. Every number in the packet comes from a table, with no hand-typed counts.

- [ ] **Step 1: Failing tests**: the packet contains each step id; the verification counts equal those in the verification table; no digit string in the generated text lacks a source table (check by regenerating and comparing; the packet is deterministic).
- [ ] **Step 2: Run; expect FAIL. Step 3: Implement and generate. Step 4: Run the tests; expect pass.**
- [ ] **Step 5: Commit** `docs(pathway_census): stage 1 review packet`. Then STOP: this is the user's review gate before any real-data run.

---

## Self-review

- **Spec coverage:** section 6.1 tables (S2); 6.2 tiers (S2 validators, S3-S5 data); 6.3 step sets and corrections (S3-S5); 6.4 steps 2-5 (S3-S5); 6.5 neutral mass and first InChIKey block (S2); review gate (S7). Not here by design: census, SIRIUS, comparison.
- **Open for the user at the gate:** the step lists, the tier assignments, the decoy lists, and whether to flip any `verified` flag.
- **Known limits:** abstract access only for most papers; the reviewer found nine yaml PMIDs and the RSY1 synonym unchecked (S3 and S6 must check them); decoy lists are PROPOSED.
