# Design amendment: a pathway knowledge base for fungal carotenoid and melanin biosynthesis (version 2)

Written 2026-11-21 by the BFD-HPCC-master agent. Version 2 revises version 1 after three reviews. Status: **APPROVED by the user 2026-11-21** (no second reviewer pass requested).
Amends and extends `2026-11-20-pathway-census-design.md` (approved 2026-11-20). The census stays as the survey engine. This amendment adds a compound layer, a bounded corpus comparison, and melanin as a second pathway.

Labels: **MEASURED** = I checked it against a database or file on 2026-11-21. **REVIEWED** = an independent reviewer agent checked it against a source (abstract read or database record read); I have not re-checked it. **REPORTED** = copied from the EBC agent's documents. **PROPOSED** = my choice, no data behind it. **UNVERIFIED** = from memory, not checked.

## What changed from version 1
| Review | File | Effect on this version |
|---|---|---|
| Fable, methods | `2026-11-21-pathway-knowledge-graph-REVIEW_FABLE.md` | Flat tables, not a graph store. M1 dropped. M2 and M3 narrowed. Detectability gate D0 added. Sterol control dropped. Identifier traps added. Tiers extended. |
| Fable, pathways | `2026-11-21-pathway-knowledge-graph-REVIEW_PATHWAYS.md` | Six step-list errors corrected. DHN-melanin rewritten as route variants. L-DOPA steps separated. Evidence tiers split by assay type and bound to organism. Missing steps added. |
| EBC agent, comparison | `EBC-KT/collab/from-ebc/2026-11-21-ebc-004-comparison-design-review.md` | Stratification, accession as the unit, yeast stratum, mirrored duplicates, label leakage, control-set selection bias. |
| User, 2026-11-21 | (conversation) | Flat tables for now; reduced plan accepted; thresholds in 6.6 accepted; melanin included. |
A fourth, domain-expert review could not run: a model safeguard stopped it. Its detectability questions went to the EBC agent instead (request 005).

---

## 1. Intent

### 1.1 The question
Which compounds can fungal enzymes make along carotenoid and melanin pathways, which fungal lineages encode the enzymes, and, where the data can test it, does a metabolite corpus show the compounds where the genomes say they are possible?

### 1.2 Why a compound layer
- The census says which enzymes a genome carries. It does not say what they could make. The compound layer links enzyme steps to reactions and compounds, so a count of enzymes becomes "this lineage could make X, by path P, with evidence tier T".
- Every link carries its source, so a language model can answer questions and cite each answer.
- Compounds carry InChIKey, formula and exact mass, which makes comparison with a metabolite corpus possible.
- The corpus is the GNPS-processed `quant_v5` data held by the EBC agent on another host. The two agents communicate through the private repository `stajichlab/EBC-KT` (`collab/`).

### 1.3 Why melanin is in scope
Melanin intermediates may be easier to detect in LC-MS than carotenoids (UNVERIFIED): scytalone, 1,3,6,8-tetrahydroxynaphthalene and L-DOPA are small and polar. Carotenoids ionise poorly in electrospray. This is untested, and the first step of the comparison is to measure it (D0, section 6.6). The census needs melanin definitions anyway (spec 5.8).

## 2. Goals, non-goals, anti-goals

### Goals
- **G1.** Verified, route-aware step sets for carotenoid, DHN-melanin and L-DOPA melanin.
- **G2.** Flat tables of compounds, reactions, steps, evidence and taxon capability (section 6.1), with an evidence tier and a source on every enzyme-to-reaction link.
- **G3.** The real-data census attached to those tables: which genomes and lineages carry which steps, with copy number, identity to the nearest characterised enzyme, and quality flags.
- **G4.** An export for the comparison: a compound table (ChEBI, InChIKey, formula, neutral exact mass, adduct m/z values) and a capability table per taxon.
- **G5.** The bounded comparison of section 6.6, with a detectability gate, a positive control, and a permutation null.
- **G6.** A written description of the diversity of potential biology across the kingdom, each statement tied to a table.
- **G7.** Everything reproducible from the repository, with labelled numbers.

### Non-goals
- Activity, affinity or kinetics prediction.
- A claim that a compound is produced because genes are present. The tables say "possible", never "produced".
- A graph database or a query service in this version. A graph can be derived from the flat tables later, if the tables prove their worth.
- Wet-lab work.
- Publishing KEGG content (section 3).
- A genome-scale metabolic model.
- Improving spectral identification. Identification levels are reported only.

### Anti-goals (the work fails even if it runs)
- **An edge with no source.** Every link needs a citation or database identifier.
- **Family similarity read as function.** A step is T3 unless a paper shows the reaction for that organism and protein. `carX` and `carT` are close relatives with different substrates (REVIEWED).
- **Gene names without the organism.** `carS` is the beta-carotene cleavage oxygenase in *Phycomyces* and a repressor of carotenoid synthesis in *Fusarium fujikuroi* (REVIEWED, PMID 22445959, 33383912). Every gene name in a table carries its organism.
- **One EC for a step whose reaction differs between species.** Phytoene desaturase does four desaturations in *Phycomyces* carB (lycopene-forming, PMID 9079885) and five in *Neurospora* AL-1 (3,4-didehydrolycopene-forming, PMID 11017770); both REVIEWED. The number of steps cannot be read from sequence. This decides whether a lineage can make torulene and neurosporaxanthin.
- **A single linear list for a pathway with route variants.** DHN-melanin has three PKS routes (6.3). A linear list would mark *Colletotrichum lagenarium* as lacking Ayg1 and *Aspergillus terreus* as defective. Both calls would be wrong.
- **A tier that travels between orthologues.** A tier holds for the named organism and protein only.
- **Mass alone read as identity.** The EBC revalidation found m/z and retention-time bin matching at chance level (REPORTED). Corpus claims need a library name, MS2, or a permutation null.
- **Isomers merged or split wrongly.** beta-carotene, lycopene and gamma-carotene are all C40H56 (MEASURED).
- **KEGG's PubChem field read as a compound ID.** MEASURED: for astaxanthin it is a substance ID; PubChem compound 10773 is sodium hippurate, SID 10773 maps to compound 5281224 (astaxanthin).
- **The UniProt keyword read as evidence.** The keyword "Carotenoid biosynthesis" lists GGPP synthases of *Saccharomyces*, *Schizosaccharomyces* and *Eremothecium* (REVIEWED, E9). Use it to find candidates only.
- **Language-model extraction without verification.** An extracted link needs a quoted sentence, its PMID, and a second check.
- **Presence read as production.** A derived "can produce" statement must give its minimum tier and which steps the census found.

## 3. Constraints

- **KEGG terms (MEASURED, kegg.jp, 2026-11-21).** The REST API is for academic use by academic users. Providing a service built on KEGG needs an academic service licence. Rule (PROPOSED): store KEGG identifiers as cross-references only, and copy no KEGG pathway or reaction content into published output.
- **Open sources.** UniProt and Rhea are CC BY 4.0 (REVIEWED). ChEBI and PubChem licences are UNVERIFIED: check before any release.
- **Identifier traps (REVIEWED unless marked):**
  - Rhea uses charged ChEBI species. Neurosporaxanthin(1-) has a mass off by 1.0073 from the neutral compound. Map each Rhea participant to its neutral parent before computing a mass.
  - The full InChIKey for scytalone differs between ChEBI and PubChem. Match on the first block (connectivity layer).
  - Rhea has no reaction for vermelone, torularhodin or the trisporic acid steps (name search). Curated steps carry their own identifiers.
  - RHEA:56772 (adonixanthin to astaxanthin) is not on the Ojima route. Do not use it as the fungal example. Rhea has beta-carotene to echinenone (RHEA:55660) and echinenone to canthaxanthin (RHEA:55664).
- **SIRIUS.** Run through the container `/bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif` (MEASURED: the file exists), with `singularity exec --bind /bigdata,/scratch`. SIRIUS 6.3.12 needs an active login for every subcommand, and concurrent processes break the shared `~/.sirius-6.3` login file (REPORTED, notes of an earlier project). Strictly one SIRIUS process at a time, project-wide.
- **NPClassifier** takes SMILES, not spectra (MEASURED from its code). The public API works (REVIEWED). A local copy needs Zenodo models and old TensorFlow (UNVERIFIED on HPCC).
- **The corpus stays on the EBC host.** Only counts and small tables cross `EBC-KT`. No spectra and no raw feature tables.
- **House rules:** `/usr/bin/python3.12`, SLURM constraint line, `$SCRATCH`, jobs of 1-1.5 h, compressed output, no `BASH_SOURCE`.
- **No personal information** in any message or file.

## 4. What is measured so far

- **UniProt (MEASURED):** 35 reviewed and 2,755 unreviewed fungal entries with keyword KW-0125. No reviewed entry for astaxanthin synthase.
- **Astaxanthin synthase (MEASURED and REVIEWED):** The 2006 paper (PMID 16416328) states the cDNA accession AX034665 and a 3,166 bp gene with 17 introns, AX034666. UniProt **Q3HR17** (gene Ast/crtS) is identical to NCBI CAC16261 over all 557 residues (REVIEWED); B5LEC8 differs at residue 60 (G in CAC16261 and Q3HR17, D in B5LEC8). The paper's *E. coli* system made canthaxanthin only, about 1% conversion; 3-hydroxylation is shown in the fungus only (REVIEWED, full text).
- **KEGG map00906 (MEASURED):** 128 compounds, 149 reactions, 61 KOs, 44 ECs. Fungal module M01042 (K15745, K17841, K17842, K17819). K23037 `CRTS, ASY`. K17739 `THNR` (EC 1.1.1.252) is on KEGG and missing from my version 1 (REVIEWED). K19580 `rppA` is a bacterial type III PKS (REVIEWED, recalled): not for fungi.
- **KEGG melanin compounds (MEASURED):** scytalone C00779, 1,3,6,8-naphthalenetetrol C04033, vermelone C01624, L-DOPA C00355, dopaquinone C00822, L-dopachrome C01693, melanin C05606, phaeomelanin C17936. KOs: K17740 `SCD1`, K00505 `TYR`, K00421 and K05909 laccase.
- **Compounds (MEASURED, KEGG):** formula and neutral exact mass for the carotenoid and melanin tables in `collab/to-ebc/` requests 001 and 003. ChEBI and PubChem agree on the InChIKey for astaxanthin.
- **SIRIUS pilot (MEASURED):** 20 spectra, module run: 4 min 49 s, 2.4 GB peak memory, exit 0. CANOPUS tables carry NPClassifier pathway, superclass and class columns plus ClassyFire classes. Two ILP solvers (CLP, CBC) failed to load on the node because of its older glibc; SIRIUS fell back to GLPK. A container run (job 32355241) tests whether this goes away.
- **Corpus (REPORTED):** 11,236 acquisition groups, 355.4 M features; 16.2% with MS2, 10.3% with an accepted predicted formula, 0.64% with a top-1 library hit; 21.5% of hits are NIST14 "Spectral Match to X" labels. FUNGAL_REF has 116 groups from 84 accessions, each with one dominant taxon: *Komagataella* 20, *S. cerevisiae* 12, genus-only "Fungi" 15, *Aspergillus* about 10, *Penicillium* about 6. 174 mirrored duplicate groups exist.

## 5. Approaches considered

| Option | Description | Trade-off |
|---|---|---|
| A | KEGG-centred | Fast; licence limits publication; gaps in downstream fungal steps |
| **B (chosen)** | Open backbone (Rhea, ChEBI, UniProt), fungal downstream steps curated from papers, KEGG IDs as cross-references | More curation; publishable; every link cited |
| C | Language-model extraction first | Fast; highest risk; use only as a helper inside B, with verification |

Storage: **flat tables** (decided by the user, 2026-11-21), with a more complex representation only if the tables prove too limiting.

## 6. Design

### 6.1 Flat tables (PROPOSED)
- `compound`: `chebi_id`, `inchikey`, `inchikey_block1`, `smiles`, `formula`, `neutral_exact_mass`, `kegg_id` (cross-reference only), `np_class`.
- `reaction`: `reaction_id` (Rhea, or a curated id), `ec` (nullable), `kegg_rn` (cross-reference only), `direction`, `curated` (bool), `note`.
- `reaction_participant`: `reaction_id`, `compound` (neutral parent), `role` (substrate or product), `stoich`.
- `step`: the census step id, route variant, `ko` (cross-reference), HMM id, `family`, `optional`.
- `step_reaction`: `step_id`, `reaction_id`, `organism`, `gene_name`, `protein_accession`, `tier`, `assay`, `pmid`, `quote`.
- `taxon_capability`: `taxon`, `step_id`, `copies`, `best_identity_to_t1`, `quality_flags`, `genomes_n`.
- Derived views only, never asserted by hand: `path_possible(taxon, compound)` with `min_tier`, `min_identity` and the list of steps found.

### 6.2 Evidence tiers (PROPOSED, from the pathway review)
- **T1a:** the fungal enzyme assayed, or its product identified after heterologous expression.
- **T1b:** gene deletion or complementation with intermediate analysis in the fungus.
- **T1c:** deletion with colour phenotype only.
- **T2:** reaction shown in another organism for the same enzyme family.
- **T3:** inferred from family or orthology.
- **T4:** hypothesis, no direct evidence.
A tier is bound to the named organism and protein. Orthologues in other species are T3 until a paper shows the reaction. A path with any T4 step is not reported as possible. The capability score reports the lowest tier and the lowest identity to the T1 enzyme along the path. No identity threshold is fixed in this version; the builder proposes one from the decoy data.

### 6.3 Pathway scope (REVIEWED corrections applied; step sets still need the user's review)

**Carotenoid, core and astaxanthin.**
- GGPP synthase; fused phytoene synthase and lycopene cyclase; phytoene desaturase recorded as **two reactions** (four-step, lycopene-forming; five-step, 3,4-didehydrolycopene-forming) with the tier attached per enzyme, T1a for *Phycomyces* carB and *Neurospora* AL-1 only; the step carries `ec: null` or both ECs.
- Astaxanthin synthase CrtS with accession **Q3HR17** (identical to CAC16261). Keep `ec: null`; KEGG ECs are cross-references only. CrtR (the P450 reductase partner, PMID 18837978) is a separate step. The astaxanthin route is a **grid**, not a line: after echinenone, 4'-ketolation or 3-hydroxylation can come first, and both reach phoenicoxanthin (REVIEWED, full text of PMID 16416328). The same enzyme converts torulene to 4-ketotorulene to 3-hydroxy-4-ketotorulene.
- Tiers: 4-ketolation T1a; 3-hydroxylation T1b (fungus only, not shown in *E. coli*).

**Carotenoid cleavage and apocarotenoids.**
- `carX` (*Fusarium*, retinal; PMID 17293483) and *Ustilago* Cco1 (PMID 19584000) as references of the carX step; `carT` and `cao-2` (torulene; PMID 17493127, 17610084); `ylo-1` and `carD` (PMID 18627463, 21749649). The order of cyclisation, cleavage and oxidation differs between *Neurospora* and *Fusarium*, so the KEGG module's straight line is not the whole story (REVIEWED, PMID 21749649).
- Decoy: *Neurospora* CAO-1 (resveratrol cleavage, PMID 23893079).

**Mucorales trisporoid steps (keyed by organism).** *Phycomyces* carS and acaA (PMID 21854466, 22445959); tsp3 (PMID 17822929); TSP2 (PMID 18931040); tdh (PMID 8828234). Later oxidations: T4. All REVIEWED (abstract read), restricted to Mucorales and to mating-type context.

**Marked T4 or excluded:** the torulene to torularhodin enzyme (no paper found by the reviewer; a few PubMed queries only), native fungal carotenoid glycosylation or esterification.

**DHN-melanin, as route variants (replaces my linear list).**
- **PKS routes (three reference clades):** (a) a true THN synthase, *C. lagenarium* PKS1 (PMID 10501004, 23261597); (b) a hexaketide PKS plus the chain-shortening hydrolase Yg1, *Wangiella* WdPKS1 and WdYg1 (PMID 11179356, 18676950); (c) a heptaketide PKS plus Ayg1, *A. fumigatus* Alb1 and Ayg1 (PMID 11350964). Ayg1/Yg1 is required with the hexaketide and heptaketide routes, not with the THN synthase route. Decoy: *A. nidulans* wA.
- **Two reductase types (separate reference clades, real substrate overlap):** 4HNR/Arp2 type (T4HN to scytalone, EC 1.1.1.252; PMID 10956664); 3HNR/THR1/BUF1/BRM2 type (T3HN to vermelone; PMID 8672814, 10956664). The census must not call one step absent when only one paralogue is present.
- **Scytalone dehydratase** (SCD1, Arp1, SDH1). The vermelone dehydration is a **curated reaction with its own evidence record**, not a copy of the scytalone edge; which enzyme acts is unsettled (REVIEWED, O4).
- **Oxidases (two families):** laccase (Abr2 type) and multicopper oxidase (Abr1 type); evidence is phenotype only, T1c (PMID 10515939, 16988817, 21277986).
- Lineage restrictions: *A. terreus* lacks DHN-melanin and makes Asp-melanin (MelA, TyrP; PMID 27133313). Botrytis and other Leotiomycetes carry two PKS (PMID 26514268).
- Add shunt products to the compound table: flaviolin (C18012, ChEBI 42646).

**L-DOPA melanin.**
- Tyrosinase (EC 1.14.18.1) and laccase (EC 1.10.3.2) are **separate steps**. A laccase does not hydroxylate tyrosine. *Cryptococcus* makes melanin only from exogenous diphenols (PMID 17098915, 8300520), so the route is substrate-dependent and the compound that enters the table depends on the supplied substrate.
- Other melanin routes as small branches or explicit exclusions: pyomelanin (HppD, HmgA; PMID 19028908), Asp-melanin.
- Dopaquinone to dopachrome onwards: non-enzymatic in the table; T4 for any fungal enzyme.

### 6.4 Build steps (PROPOSED)
1. Candidate enzymes from UniProt reviewed fungal entries by keyword and EC. **Candidates only**; each is tiered by its own paper.
2. Reactions from Rhea by EC and ChEBI; map every participant to its neutral parent.
3. Compounds from ChEBI for InChIKey, formula and mass; confirm with PubChem by InChIKey; match on the first block. Never join through KEGG's PubChem field.
4. KEGG identifiers as cross-references only.
5. Curate downstream fungal steps from papers with PMID, assay type and a quoted sentence.
6. Optional language-model extraction, each link verified a second time.
7. NPClassifier on each compound's SMILES.

### 6.5 Exports for the comparison
- **Compound table:** `chebi_id, inchikey_block1, formula, neutral_exact_mass, adducts`. Adducts for positive mode: [M+H]+, [M+Na]+, [M+NH4]+, [M+K]+, [M-H2O+H]+, [M+2H]2+, **and the radical cation M+.**, which carotenoids give in ESI (REVIEWED, PubMed). Negative mode: [M-H]-.
- **Capability table:** per taxon and compound, the best path, its lowest tier and lowest identity, and which steps the census found.

### 6.6 Comparison, reduced (PROPOSED, thresholds accepted by the user 2026-11-21)

**D0, the detectability gate (first).** A compound goes into M2 only if:
- (a) at least 6 capable accessions exist in the Orbitrap-positive filamentous stratum;
- (b) at least 3 of them show the compound at library-name level or at formula plus MS2;
- (c) the positive-control cohort shows it at library-name level in at least 3 accessions (carotenoids: plant, algae and cyanobacteria groups; also the known *Penicillium*/*Aspergillus* products as a check on the whole pipeline).
If (c) fails, carotenoids are undetectable in this corpus and a fungal null means nothing. The EBC agent is producing the counts (requests 001, 003, 005).

**M2, capability against observation (the only model in the first round).**
- Stratum: Orbitrap positive only. Orbitrap negative and qTOF positive are reported separately when a cell has at least 6 capable and 30 non-capable accessions (PROPOSED; the power table assumes 100 non-capable and must be recomputed with the real number). "Unknown" instrument groups stay out of the primary test.
- Unit: the accession, not the group. Hold out whole accessions.
- Filamentous taxa form the primary analysis. Yeasts (*Saccharomyces*, *Komagataella*) are a separate stratum, because 32 of 116 fungal groups are yeast and yeast capability is confounded with growth form and extraction. Genus-only "Fungi" groups (15) are dropped. Mirrored duplicate groups are collapsed first.
- Two definitions of "observed", reported separately: library name, and accepted formula plus MS2 class. Formula-only counts are MSI level 4 at best.
- Controls: positive taxa and the known-product check; negative taxa are those whose genomes lack the pathway (I confirm from the census, not from memory). A null for melanin or carotenoids says "not detected under the culture condition"; these pathways are often silent in standard media.
- Null: permute capability across taxa within instrument and polarity (single-taxon studies degenerate; the permutation is across accessions). BH within the cell. The family of tests is fixed before looking.
- **Minimum effect (what we can detect):** the observation rate among capable accessions at least 45%, at least 40 percentage points above the non-capable accessions; permutation p at most 0.01, BH q at most 0.05, same direction when any one accession is dropped. Anything weaker is reported as "not testable", never "no effect". Power table (MEASURED, exact one-sided Fisher, 100 non-capable accessions, 80% power): with 6, 10 and 16 capable accessions and a 5% background rate, the capable rate must be at least 50%, 39% and 32% at alpha 0.05, and 63%, 48% and 39% at alpha 0.01.

**M3, class anchors (second round, only after the control set).** CANOPUS classes for MS2 spectra against NPClassifier classes of the table's compounds. Valid only if a matched negative-control set is run through SIRIUS in the identical way, because the 734 spectra already sent were chosen from fungal marker tables (EBC review). Coarse classes only (pathway, superclass). Melanin classes are too broad to be useful on their own (REVIEWED). Record class-system versions. Mark features with default adducts (59) and mirror duplicates (22).

**M1, mass shifts: dropped.** The shifts (H2, O, H2O) occur everywhere in LC-MS data, cyclisation has zero shift, and isotope, adduct and in-source-loss pairs produce the same differences (REVIEWED and EBC). It may return later as an exploratory check only.

**The sterol branch is not used as a control.** It is not a valid M2 control (REVIEWED).

### 6.7 Melanin-specific points
- DHN-melanin is a polymer. Only monomers and intermediates are testable. Intermediates accumulate in blocked or inhibited strains (REVIEWED, PubMed).
- L-DOPA is not specific: plants and animals make it, media contain it, and in *Cryptococcus* the fungus does not make it (PMID 17098915).
- Dopaquinone and dopachrome may be unstable analytes (UNVERIFIED; asked in request 005).
- Detectability is measured (D0) before any melanin model runs.

### 6.8 Collaboration protocol
The rules in `EBC-KT/collab/README.md` apply. The BFD side signs as `BFD-HPCC-master agent`. A guarded sync script pulls with rebase and refuses to push anything but new files under `collab/to-ebc/`; it never forces a push (history audited 2026-11-21: every change to `collab/` is an addition).

### 6.9 Review status
Done: Fable (methods), Fable (pathways), EBC agent (comparison design, from the summary). Pending: the user's approval of this version. A second pass on the revised step sets by the same reviewers is optional.

## 7. Order of work (PROPOSED)
1. D0: the EBC counts, including the plant and algae positive control and the melanin detectability check. Nothing else in the comparison runs before this answers.
2. Stage 1: verified step sets and references for all three pathways, with decoys (the corrections in 6.3), and the flat tables.
3. Stage 2: the census on real data for the verified step sets (the real-data plan, revised).
4. Stage 3: M2 on the passing compounds.
5. Stage 4: the matched control set through SIRIUS in the container (serially), then M3.
6. Stage 5: Chapters A (carotenoid) and B (melanin). A negative result is reported as such.
In parallel now: the SIRIUS container pilot, and the full SIRIUS run on the 734 inputs, in one serial job or a serial chain.

## 8. Risks
- Carotenoids may be undetectable (D0 decides).
- Only about 13 usable taxa for M2 (Fable's estimate); a null is likely uninterpretable for weak effects.
- Culture extracts include medium and pellet; the corpus cannot remove that confound.
- Evidence tiers depend on curation. A wrong T1 corrupts every downstream claim. The reviewer read abstracts only, apart from the Ojima paper; nine yaml PMIDs and the RSY1 synonym were not checked.
- A shared login file forces SIRIUS to run serially, which makes the full run slow.
- Licence terms for a shared table with KEGG cross-references; ChEBI and PubChem licences unverified.

## 9. Open items for the user
1. Whether the tables will be released outside the group. The licence rule changes if so.
2. Whether to ask the reviewers for a second pass on the corrected step sets before stage 1.
3. Which fungal taxa to use as negative controls, once the census can confirm the pathway is absent (candidates only so far: *Saccharomyces*, *Komagataella*).

## 10. Left to the builder
Storage details of the flat tables (DuckDB or parquet), the capability score formula and the identity threshold, adduct tolerance, HMM and tree details, and the order within stages.
