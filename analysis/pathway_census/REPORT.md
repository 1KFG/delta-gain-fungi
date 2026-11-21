# Pathway knowledge base and the metabolite corpus: what was tried, what failed, what to do next

Status: in progress. Last updated 2026-11-21. Code: `analysis/pathway_census/`.

::: {.callout-note}
## What this is, what it is not, and where it stops being reliable

This is a working record of explorations done on 2026-11-20 and 2026-11-21. It is written so that a person or an agent can use the corpus, or connect it to other data, without repeating the failures.

**Nothing here is a biological finding.** No bench work was done. The comparison between fungal genes and metabolite observations was not run, because the pre-set detectability gate failed for every pathway tested in fungi. The curation (Stage 1) is not finished: Tasks 1 to 5 are closed, Task 6 (mycosporines) is running, and no step set has had the owner's review.

Every number carries a label. **MEASURED**: this agent ran or fetched it on 2026-11-21. **REPORTED**: taken from a message of the EBC agent (the agent that holds the corpus), not re-measured here. **PROPOSED**: a choice, not a fact. **UNVERIFIED**: not checked. Numbers from the EBC agent use the corpus version `quant_v5`.
:::

## 1. The question and why it changed

The first request was two analysis chapters for the fungal genome census: carotenoid biosynthesis and melanin biosynthesis (DHN and L-DOPA), modelled on the earlier peroxygenase (UPO) project. Structure and affinity prediction was dropped. Gene copies were kept as a ranking for later functional testing.

The request then grew. The owner asked for a public-knowledge **knowledge base**, as flat tables, of what fungal enzymes can make (upstream and downstream), and a comparison with the processed metabolite corpus. Mycosporines and mycosporine-like amino acids (MAAs) were added as a fourth pathway.

The working question was: **if a fungal genome carries the genes for a pathway, does the metabolite corpus show the product more often for that species?**

## 2. The two data sets

| Side | What it holds | Source |
|---|---|---|
| Genome side (BFD) | Fungal genomes only. Counts for taxa named by the EBC notes: *Batrachochytrium* 15 genomes (2 species), *Pleurotus* 53, *Agaricus* 20, *Naganishia* 28, *Rhodotorula* 220, *Phaffia* 12, *Penicillium* 655, *Aspergillus* 1,632 (MEASURED, `species.parquet`, 2026-11-21). | this repository |
| Corpus side | 11,236 acquisition groups and 355.4 M features. 16.2% of features have MS2, 10.3% have an accepted predicted formula, 0.64% have a top-1 library hit. 21.5% of the hits are NIST14 "Spectral Match to X" labels (REPORTED). `FUNGAL_REF` holds 116 groups from 84 accessions; 32 of the 116 are yeast (REPORTED). | EBC agent, `quant_v5` |

A genome count says a genome exists. It says nothing about pathway genes. BFD has no bacterial, plant or animal genomes, so a bacterial skin-microbe pathway cannot be placed on the genome side.

## 3. What was built (genome side)

### 3.1 Tool and tables

- A Python package, `analysis/pathway_census/` (merged to `main`, merge `1122d5f`). It reads a pathway definition (YAML, schema v2: routes, any-of domains, decoys, regions, organism-keyed references) and counts genes per genome.
- Flat TSV tables, not a graph store: `compound`, `reaction`, `reaction_participant`, `step`, `step_reaction`. Each link between a step, a reaction and a protein carries an **evidence tier** bound to the organism and the protein.
- Tiers (spec 6.2): T1a (the fungal enzyme assayed, or its product identified after heterologous expression), T1b (gene deletion or complementation with intermediate analysis in the fungus), T1c (deletion with colour phenotype only), T2 (reaction shown in another organism for the same enzyme family), T3 (inferred from family or orthology), T4 (hypothesis, no direct evidence). Curation also applies a stricter quote rule for T1a (section 5.1).
- Route-aware completeness: `path_possible` closes each step per organism and reports `organisms`, `single_organism` and `n_equal_paths`. A route tier that holds only across several organisms is not a tier for any one organism.

### 3.2 State of curation (MEASURED, review reports in the ledger)

| Pathway | Task | State | Headline |
|---|---|---|---|
| Carotenoid | 3 | closed after 3 fix rounds | 29 verified references, 61 decoy rows, 36 links at the first report; tiers changed during review (section 5) |
| DHN-melanin | 4 | closed after 1 fix round | Review found the two DHN reductases swapped in the pathway review |
| L-DOPA and other melanins | 5 | closed after 1 fix round | 15 references (14 verified), 10 decoy proteins, 29 links: 2 T1a, 3 T1b, 2 T1c, 20 T3, 2 T4 |
| Mycosporines and MAAs | 6 | running | Starts from literature leads; no fungal gene list is assumed |
| Verifier and review packet | 7, 8 | not started | Then a stop for the owner's review |

Each pathway has an `UNVERIFIED.md` under `analysis/pathway_census/refs/<pathway>/` that lists every item a person must decide.

## 4. The comparison design and what happened to it

### 4.1 Four models, three kept

| Model | Idea | Outcome |
|---|---|---|
| M1 | Mass shifts between features (H2, O, H2O) mark a biosynthetic step | **Dropped before any run.** These shifts occur everywhere in LC-MS data. Cyclisation has zero shift. Isotope, adduct and in-source-loss pairs give the same differences (Fable review and EBC review) |
| M2 | Capability (genome has the pathway) against observation (compound seen) | Gate failed for every pathway tested in fungi (section 4.3) |
| M3 | CANOPUS and NPClassifier class anchors for MS2 spectra | Run on 734 spectra; comparison waits for a matched control set (section 4.4) |
| Methods check | Test the method on a pathway that is detectable | Fucose fails, MAAs have no fungal signal; fungal natural products (for example roquefortine, andrastin) not tried |

### 4.2 The detectability gate (set before any counts arrived)

A compound enters M2 only if all three hold (PROPOSED, accepted by the owner):

1. At least 6 capable accessions in the Orbitrap-positive filamentous stratum.
2. At least 3 of them show the compound by library name, or by formula plus MS2.
3. The positive-control cohort shows it by library name in at least 3 accessions.

The minimum effect, also fixed in advance: observation rate among capable accessions at least 45% and at least 40 percentage points above non-capable accessions; permutation p at most 0.01; Benjamini-Hochberg q at most 0.05; same sign when any one accession is dropped. A weaker result is "not testable", not "no effect". Exact one-sided Fisher power (80% power, 100 non-capable accessions, 5% background, alpha 0.05): the capable rate must reach 50%, 39% and 32% with 6, 10 and 16 capable accessions.

The thresholds were not lowered after the counts arrived.

### 4.3 Detectability results (REPORTED by the EBC agent; the gate reading is mine)

| Target | Count | Gate reading |
|---|---|---|
| β-carotene, zeaxanthin, canthaxanthin | 0 of 84 fungal-reference accessions each | Fail |
| Astaxanthin, retinal, torulene | 1, 2 and 0 of 84 | Fail |
| Melanin library names | 5 of 7 names have zero hits; L-DOPA 0 fungal-reference groups; THN 1 | Fail |
| Melanin formulas | Present in 10 to 35% of groups, in every cohort | Background, not evidence (MSI level 4 at best) |
| Free fucose and rhamnose | 17 library features in 13 accessions; none in fungal-reference or milk accessions | Fail; curation not started |
| MAAs in fungi | 0 of 82 fungal-reference accessions | Fail |
| MAAs overall | 36 of 87 marine accessions carry an MAA library name, against 216 of 5,379 others; odds ratio 16.9; permutation p 0.0005 | A robust marine positive control |

The MAA enrichment holds after removing the top three contributing accessions, the largest marine study and all fungal-labelled accessions (odds ratio 15.7 when restricted to precursors within 5 mDa). It shows that MAAs are present in marine samples. It does not show who made them: cyanobacteria, algae, fungi, and animals that take MAAs up all occur there (UNVERIFIED).

::: {.callout-important}
## Two traps in these numbers

The operator note quoted 32 of 88 marine accessions (odds ratio 14.5). The later, filtered measurement is 36 of 87 (odds ratio 16.9). Use the later one. The MAA query counts 82 fungal-reference accessions and the carotenoid query counts 84. The reason for the difference was not stated; it is **UNVERIFIED** and should not be reconciled by guessing.
:::

### 4.4 SIRIUS on the 734 EBC feature spectra (MEASURED)

A formula, fingerprint, CANOPUS and structure-search chain ran in a container (SIRIUS 6.3.12), strictly serially, on 734 spectra chosen by the EBC agent from fungal marker tables.

| Stratum | Spectra | With formula | With structure | With CANOPUS |
|---|---|---|---|---|
| Orbitrap negative | 176 | 176 | 154 | 145 |
| Orbitrap positive | 520 | 520 | 445 | 404 |
| qTOF positive | 38 | 38 | 33 | 25 |
| All | 734 | 734 | 632 (86.1%) | 574 (78.2%) |

Run time was 8,616 s in total (four chunks: 2,917; 3,636; 1,586; 477 s). The counts are counts of output rows. They do not measure correctness.

Only 3 of the 734 features carry a pigment-related class label. Two read "Carotenoids (C40)" at superclass level but have the formula C28H40, which cannot be a C40 skeleton as written; the class label is a model output and should not be read as a carotenoid detection. The third is a quinazoline alkaloid. These are **leads for a person to look at**, not hits.

M3 is valid only if a matched negative-control set goes through SIRIUS in the identical way, because the 734 spectra were chosen from fungal marker tables. That set has not been requested.

### 4.5 Amphibian skin chemistry (REPORTED, census only)

1,778 amphibian-labelled runs in 24 accessions. 1,586 runs in 23 accessions are host samples, mostly laboratory *Xenopus* eye, brain and testis. 99 runs in one accession are *Batrachochytrium* cultures. Skin, mucus, gland and swab labels: 0 runs. A skin-community screen is not possible in this corpus. 12 of the 24 accessions have usable feature tables.

## 5. The negative results and the errors, each with what caught it

### 5.1 Scientific and data errors found by review

| What went wrong | What caught it | Rule that came out of it |
|---|---|---|
| In the pathway review, the two DHN reductases were swapped (Q9HFV6 is the 4HNR, Q12634 the 3HNR) | Task 4 implementer and Fable reviewer, by fetching UniProt records | Every accession is re-fetched; a recalled accession is a lead |
| KEGG "PubChem" field is a substance ID, not a compound ID (astaxanthin SID 10773 is CID 5281224) | Compound check against ChEBI and PubChem | Match on the first block of the InChIKey |
| Rhea reactions use charged ChEBI species; the tables use neutral parents | Table validation | Neutral parents only; masses from formulas |
| An alias rule from memory merged *Mucor lusitanicus* (taxid 29924) with *M. circinelloides* (36080) | Independent NCBI taxid check in review | An alias only where NCBI gives both names one taxid |
| Q2UP46 was carried as the assayed melB protein; the assayed sequence is A0A1S9DK56 (matches PDB 3W6W at 615 of 616 positions; Q2UP46 at 585) | Task 5 reviewer compared with the deposited structure | Check the sequence against a structure when a paper's residue numbers are used |
| 7 of 9 T1a links in Task 5 quoted sentences that did not name the protein or the substrate | Task 5 reviewer | T1a only when the quote names the link's own protein and an assay or product; otherwise T3 |
| Q9Y8G7 (a *Fusarium* P450 fused to a P450 reductase domain) was listed as a decoy of the `crtr` reductase step, but it carries the reductase domains and so can do that reaction | Task 3 review | Removed from the `crtr` decoys; kept as a decoy of astaxanthin synthase only. A decoy must not carry the reference reaction |
| The pathway review read abstracts only (apart from one paper); nine YAML PMIDs and one synonym were not checked | Spec risk section, then full re-fetch | Every PMID and quote is fetched and matched to the abstract text |
| Phytoene desaturase tiers in the review list differed from what the abstracts show by organism (four-step T1a for three organisms; five-step T1a for *Neurospora* AL-1 only, T1b for *Fusarium* carB) | Task 3 | Tier changes written into the spec (commit `d29d73c`); awaiting the owner's decision |

### 5.2 Operational failures

| What went wrong | Effect | Fix |
|---|---|---|
| SIRIUS 6.3.12 exits 0 when its tool chain fails | Two "finished" chunks took 7 and 10 s and held no results | Check for result tables; fail the job when none appear |
| The output directory was created before SIRIUS ran | SIRIUS stops with "is a directory, must be a file" and still exits 0 | Let SIRIUS create the project path |
| Pilot outputs were written to node-local `$SCRATCH` and lost | Pilot results gone | Full run writes to shared storage |
| Module run on a node with old glibc | CLP and CBC solvers failed to load; fell back to GLPK | Run in the container |
| One shared login file | SIRIUS runs strictly serially | Serial chain of four jobs |
| A model safeguard stopped an Opus pathway-biology review | The review was not done by that agent | Split to Fable and the EBC agent; briefs say "database curation, no wet-lab" |
| Usage limit interrupted a fix round | Round resumed from saved work | Saved uncommitted work and resumed by message |

## 6. Hypotheses and their status

| # | Hypothesis | Status | Evidence |
|---|---|---|---|
| H1 | Fungal genomes with the carotenoid pathway show carotenoids more often in the corpus | **Not testable here** | 0 of 84 accessions show the main carotenoids; the gate fails |
| H2 | Same for melanin (DHN and L-DOPA) | **Not testable here** | Section 4.3 |
| H3 | Same for MAAs in fungi | **Not testable here** | 0 of 82 accessions; the marine signal does not name a producer |
| H4 | Mass-shift pairs reveal pathway steps | **Rejected before use** | Shifts are not specific |
| H5 | CANOPUS class anchors match pathway products | **Open** | Needs the matched control set; coarse classes only |
| H6 | The method works on some detectable fungal pathway | **Open** | Fucose fails; fungal natural-product pathways not yet tried |
| H7 | Pathways are often silent in standard culture, so a null means "not seen under these conditions" | **Plausible, not tested** | Reason a null cannot be read as absence |

Competing explanations for the carotenoid and melanin nulls, none resolved: the compounds are not extracted or ionised by the methods in the corpus; the cultures do not express the pathway; the compounds are polymeric or unstable (melanin) so the monomer is absent; the library lacks the names; the fungal reference set is small and dominated by yeasts that do not make these compounds.

## 7. How to connect this corpus to other data

1. **The unit is the accession, not the group.** Hold out whole accessions. Collapse the 174 mirrored duplicate groups first (REPORTED).
2. **Stratify by instrument and polarity.** Orbitrap positive is the primary stratum. Orbitrap negative and qTOF positive are reported separately only when a cell has at least 6 capable and 30 non-capable accessions (PROPOSED).
3. **Separate yeast from filamentous fungi.** Capability, growth form and extraction method are confounded.
4. **Drop genus-only "Fungi" groups** (15 of 116).
5. **Define "observed" two ways and report both:** library name, and accepted formula plus MS2 class. Formula alone is MSI level 4 at best.
6. **Never match on mass alone.** Isobars exist (palythene and usujirene share m/z 285.14; mycosporine-methylamine-threonine and mycosporine glutaminol share 303.154).
7. **Library labels need care.** 21.5% of hits are NIST14 spectral-match labels. Two MAA-named groups with a fungal dominant taxon hold non-fungal runs, so the hit cannot be assigned to the fungus.
8. **Run the gate first.** Get the detectability counts before building a model. This saved a full comparison here.
9. **Run any SIRIUS comparison with a matched control set** processed the same way.
10. **Keep evidence tiers attached to claims.** A T3 link is an annotation, not an assay. A route tier that holds only across several organisms is not a tier for any one organism.
11. **Identifiers:** KEGG ids are cross-references only (academic licence). ChEBI and PubChem licences for a shared table are UNVERIFIED. Release outside the group is an open owner decision.
12. **Genome gaps:** BFD holds fungi only. Bacterial, algal, plant and animal producers need genomes from elsewhere.

## 8. What would change these conclusions

- A compound with at least 6 capable and 3 observing accessions in the primary stratum, and a positive control at 3 accessions, would open M2 for that compound.
- Metabolomics of fungi grown for a pigment (induced carotenoid or melanin conditions, targeted extraction) would test H1 and H2 directly. This is outside the corpus.
- A matched negative-control SIRIUS run would decide whether class anchors (M3) mean anything.
- Fungal genomes with confirmed MAA genes, once Task 6 finishes, would let the marine signal be compared with a fungal capability call.

## 9. Not done here

- No M2 or M3 test was run. No census was run on real data for any curated step set.
- No bench work and no chemical identification.
- No full-text reading by the curation agents except where the report says so; abstracts and database records were read.
- BUDDY against SIRIUS formula comparison; matched control set; the fungal natural-product methods check.
- Tasks 6 to 8 of Stage 1; the owner's review of step sets.
- KEGG KO and compound ids, some ChEBI records and the TIGRFAM ids on two entries were not independently verified in Task 5.

## 10. Decisions waiting for the owner

The full lists are in the `UNVERIFIED.md` file of each pathway. In outline:

- Carotenoid: the acaA verified flag; the phytoene desaturase tier changes; the 3-hydroxylation and grid-branch tiers of astaxanthin synthase.
- DHN-melanin: the Q9HFV6 gene-name rule; vermelone dehydration T1a on P56221; the two routes with no product compound; the *Colletotrichum* and *Glomerella* organism key; abr1 at T3.
- L-DOPA melanin: Q2TL94 (no gene name); the melB reference A0A1S9DK56 (ORF name only, tied by structure); the P00440 tie; the LAC1 T1c sentence without a substrate; TyrP at T3; no product for the Asp-melanin route; the rule for entries with an ORF name only.
- Project: whether to release tables outside the group; the negative-control taxa; whether to run the fungal natural-product methods check.

## 11. Where things are

| What | Where |
|---|---|
| Census spec | `docs/superpowers/specs/2026-11-20-pathway-census-design.md` |
| Knowledge base spec (v2 and addendum 2.1) and its reviews | `docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-design.md`, `...-REVIEW_FABLE.md`, `...-REVIEW_PATHWAYS.md` |
| Stage 1 plan | `docs/superpowers/plans/2026-11-21-stage1-verified-step-sets.md` |
| Exploration log (status, counts, restart information) | `docs/superpowers/notes/2026-11-21-exploration-log.md` |
| Package, pathways, references, tables | `analysis/pathway_census/` (`pathways/`, `refs/`, `kb/`, `curation_scripts/`, `tests/`) |
| SIRIUS runs and merged tables | `analysis/ms_comparison/sirius_full/merged/` on the main checkout |
| Corpus-side messages and count tables | private repository `stajichlab/EBC-KT`, folder `collab/` (counts and small tables only) |
| Branch | `pathway-census-realdata` on the project remote; not merged |
