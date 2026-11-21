# Design amendment: a pathway knowledge graph for fungal carotenoid and melanin biosynthesis

Written 2026-11-21 by the BFD-HPCC-master agent. Status: **DRAFT, awaiting three reviews and the user's approval.**
Amends and extends `2026-11-20-pathway-census-design.md` (approved 2026-11-20). The census stays as the survey engine. This amendment adds a compound layer, a corpus comparison, and melanin as a second pathway.

Labels: **MEASURED** = checked by me against a database or file on 2026-11-21. **REPORTED** = copied from another agent's document (the EBC agent's `quant_v5_handoff.md` and follow-up plan); not re-checked. **PROPOSED** = my choice, with no data behind it. **UNVERIFIED** = from memory or not checked.

---

## 1. Intent

### 1.1 The question
Which compounds can fungal enzymes make along carotenoid and melanin pathways, which fungal lineages encode the enzymes, and does a metabolite corpus show the compounds where the genomes say they are possible?

### 1.2 Why a graph
- The census says which enzymes a genome carries. It does not say what they can make. The graph links enzymes to reactions and compounds, so a count of enzymes becomes "this lineage could make X, by path P, with evidence tier T".
- Every link carries its source. A language model can then answer questions over the graph and cite each answer.
- The compound layer carries InChIKey, formula and exact mass. That makes comparison with a metabolite corpus possible.
- The corpus is the GNPS-processed `quant_v5` data held by the EBC agent on another host. The two agents communicate through the repository `stajichlab/EBC-KT` (`collab/` folder, rules in its README).

### 1.3 Why melanin is in scope
Melanin intermediates may be easier to detect in LC-MS than carotenoids. From memory (UNVERIFIED), scytalone, 1,3,6,8-tetrahydroxynaphthalene and L-DOPA are small polar molecules; carotenoids ionise poorly in electrospray. The EBC agent expects few or no carotenoid hits (its stated expectation, PROPOSED). If melanin intermediates are detectable, melanin is the better test of the whole comparison. The census also needs a melanin definition anyway (spec 5.8).

## 2. Goals, non-goals, anti-goals

### Goals
- **G1.** Verified step sets for three pathways: carotenoid (upstream and downstream), DHN-melanin, and L-DOPA melanin.
- **G2.** A graph of compounds, reactions, enzyme families, evidence and taxa (section 6.1), with an evidence tier on every enzyme-to-reaction link and a source on every edge.
- **G3.** The census (real data) attached to the graph: which genomes and lineages carry which enzyme families, with copy number and quality flags.
- **G4.** An export for the comparison: a compound table (ChEBI, InChIKey, formula, exact mass, adduct m/z values) and a capability table per taxon.
- **G5.** Three comparison models (section 6.6), each with negative controls and a permutation null.
- **G6.** A written description of the diversity of potential biology across the kingdom, with every statement tied to a table or figure.
- **G7.** Everything reproducible from the repository. MEASURED, REPORTED, PROPOSED and UNVERIFIED are labelled.

### Non-goals
- Activity, affinity or kinetics prediction. The UPO result (spec section 2) still holds.
- A claim that a compound is produced because genes are present. The graph says "possible", never "produced".
- Wet-lab work.
- Publishing KEGG content (section 3).
- A genome-scale metabolic model.
- Spectral identification as a deliverable. Identification levels are reported, not improved.

### Anti-goals (the work fails even if it runs)
- **An edge with no source.** Every edge needs a citation or database identifier.
- **Family similarity read as function.** `carX` and `carT` are close relatives with different substrates. A step is T3 (inferred) unless a paper shows the reaction.
- **Mass alone read as identity.** The EBC revalidation showed m/z and retention-time bin matching was at chance level (REPORTED). Corpus claims need a library name, MS2 or a null.
- **Isomers merged or split wrongly.** beta-carotene, lycopene and gamma-carotene are all C40H56, exact mass 536.4382 (MEASURED). A formula match cannot choose between them.
- **KEGG's PubChem field read as a compound ID.** MEASURED: KEGG lists astaxanthin as "PubChem: 10773". That is a substance ID (SID). PubChem compound 10773 is sodium hippurate. SID 10773 maps to compound 5281224, which is astaxanthin. The join must go through ChEBI or the InChIKey.
- **Language-model extraction without verification.** An extracted edge needs a quoted sentence, its PMID, and a second check before it enters the graph.
- **Presence read as production.** The derived edge `can_produce` must state its minimum evidence tier and which steps are present.

## 3. Constraints

- **KEGG terms (read from kegg.jp on 2026-11-21, MEASURED).** Academic users may use the website freely; the REST API is for academic use by academic users; providing a service built on KEGG needs an academic service licence in the FTP subscription; non-academic use needs a commercial licence. Whether a shared graph that carries KEGG content counts as a service is a question for the institution. Design rule (PROPOSED): store KEGG identifiers as cross-references only and copy no KEGG pathway or reaction content into any published output.
- **Open sources (PROPOSED backbone):** UniProt, ChEBI, Rhea, PubChem. Their licences are UNVERIFIED: I could not read the licence statements by script. Verify before any release.
- **Rhea is reachable (MEASURED):** `https://www.rhea-db.org/rhea?query=...&format=tsv&columns=...` returns reactions. Example: RHEA:56772, all-trans-adonixanthin + 2 AH2 + 2 O2 = all-trans-(3S,3'S)-astaxanthin + 2 A + 3 H2O.
- **MetaCyc** needs a licence for bulk data (UNVERIFIED). Not used.
- **Compute and house rules:** the pilot constraints stand (`/usr/bin/python3.12`, SLURM constraint line, `$SCRATCH`, job size 1-1.5 h, compressed output, no `BASH_SOURCE`).
- **SIRIUS/CANOPUS:** installed as modules 5.8.1 and 6.3.12 (MEASURED). The login is not active at 09:30 on 2026-11-21 (MEASURED) and a person must renew it.
- **NPClassifier:** its input is a SMILES string, not a spectrum (MEASURED from its code). Models come from Zenodo record 5068687 and the repo pins old TensorFlow. Whether it runs here is UNVERIFIED.
- **The corpus is not on this cluster.** Only counts and small tables cross the `EBC-KT` repo. No spectra and no raw feature tables.
- **No personal information** in any message or file.

## 4. What is measured so far

- **UniProt (MEASURED, 2026-11-21):** fungal entries with keyword "Carotenoid biosynthesis" (KW-0125): 35 reviewed, 2,755 unreviewed. Reviewed entries cover GGPP synthase, the fused phytoene synthase/cyclase, phytoene desaturase, `carX`, `carT`/`cao-2`, and `ylo-1`/`carD`. There is no reviewed entry for astaxanthin synthase.
- **Astaxanthin synthase (MEASURED, from the 2006 paper and NCBI):** CAC16261 is the protein from cDNA AX034665; the paper (PMID 16416328, Ojima et al.) states the cDNA accession and a 3,166 bp genomic gene AX034666, and shows a single P450 doing 4-ketolation and 3-hydroxylation. B5LEC8 differs from CAC16261 at one of 557 residues.
- **KEGG carotenoid map00906 (MEASURED):** 128 compounds, 149 reactions, 61 KOs, 44 ECs. Fungal module M01042 (neurosporaxanthin: K15745, K17841, K17842, K17819). K23037 `CRTS, ASY` is astaxanthin synthase.
- **KEGG melanin coverage (MEASURED):** compounds exist for scytalone (C00779), 1,3,6,8-tetrahydroxynaphthalene (C04033), vermelone (C01624), L-DOPA (C00355), dopaquinone (C00822), L-dopachrome (C01693), melanin (C05606) and phaeomelanin (C17936). KOs: K17740 `SCD1` scytalone dehydratase (EC 4.2.1.94), K00505 `TYR` tyrosinase, K00421 and K05909 laccase. KEGG has no dedicated melanin pathway map; melanin lies under tyrosine metabolism (map00350). I did not find a polyketide synthase KO for the DHN pathway (UNVERIFIED; my search output was truncated).
- **Compounds, formula and exact mass (MEASURED, KEGG; ChEBI IDs in request 001):** beta-carotene 536.4382, torulene 534.4226, echinenone 550.4175, canthaxanthin 564.3967, astaxanthin 596.3866, zeaxanthin 568.4280, retinal 284.2140.
- **ChEBI and PubChem agree (MEASURED):** astaxanthin ChEBI 40968 and PubChem compound 5281224 give the same InChIKey, MQZIGYBFDRPAKN-UWFIBFSHSA-N.
- **Corpus (REPORTED):** 11,236 acquisition groups, 355.4 M features. 16.2% of features have MS2, 10.3% an accepted predicted formula, 0.64% a top-1 library hit. 21.5% of hits are NIST14 "Spectral Match to X" labels. Peak areas are group-local. Naive fungal-marker claims were retracted after a permutation null.

## 5. Approaches considered for building the graph

| Option | Description | Trade-off |
|---|---|---|
| A | KEGG-centred: use map00906 and melanin KOs as the graph | Fast. Licence limits publication. Gaps in downstream fungal steps (no Mucorales trisporic acid step found, no torularhodin step found; UNVERIFIED for others). |
| **B (recommended)** | Open backbone (Rhea reactions, ChEBI compounds, UniProt enzymes) plus fungal downstream steps curated from papers; KEGG IDs as cross-references only | More curation. Publishable. Every link cited. |
| C | Language-model extraction from the literature first | Fast to fill. Highest risk of wrong edges. Use only as a helper inside B, with verification. |

B is my recommendation. It follows the licence rule and fills the gaps KEGG leaves.

## 6. Design

### 6.1 Graph schema (PROPOSED)
- **Compound:** `chebi_id`, `inchikey`, `smiles`, `formula`, `exact_mass`, `kegg_id` (cross-reference only), `np_class` (from NPClassifier on the SMILES).
- **Reaction:** `rhea_id` (or a curated id when Rhea has none), `ec`, substrates and products (Compound), `kegg_rn` (cross-reference only).
- **EnzymeFamily:** the census step id (from the pathway YAML), `ko` (cross-reference), HMM id, reference accessions.
- **Evidence:** source type (paper or database), id (PMID or accession), organism, assay type (in vivo complementation, heterologous expression, in vitro, inferred), tier.
- **Taxon and GenomeCapability:** the census result per genome (copies, quality flags) rolled up to taxon.
- **Edges:** `substrate_of`, `product_of`, `catalyses` (EnzymeFamily to Reaction, carries Evidence), `present_in` (EnzymeFamily to Taxon, carries copy counts), and the derived `can_produce` (never asserted by hand).

### 6.2 Evidence tiers (PROPOSED)
- **T1:** reaction shown in a fungus (paper cited).
- **T2:** reaction shown in another organism for the same enzyme family.
- **T3:** inferred from family or orthology.
- **T4:** hypothesis, no direct evidence.
A `can_produce` path reports its lowest tier and the steps that are present. A path with any T4 step is not reported as possible.

### 6.3 Pathway scope and starting step sets (UNVERIFIED; the domain-expert review must check each)
- **Carotenoid:** the pilot's 8 steps plus downstream additions: astaxanthin synthase reductase partner (crtR, PMID 18837978), the Mucorales trisporic acid cleavage steps, *Rhodotorula* torulene to torularhodin (enzyme unknown to me), and fungal carotenoid glycosylation if any.
- **DHN-melanin:** polyketide synthase (PKS1 type) to 1,3,6,8-tetrahydroxynaphthalene, then a reductase (Arp2/THR1 type) to scytalone, scytalone dehydratase (Scd1) to 1,3,8-trihydroxynaphthalene, a second reductase to vermelone, Scd1 again to 1,8-dihydroxynaphthalene, then a laccase or multicopper oxidase for polymerisation. Corrections from the earlier review stand: Rsy1 is scytalone dehydratase, not a separate step; Ayg1 belongs only to the heptaketide route; *A. terreus* lacks DHN-melanin while relatives have it.
- **L-DOPA melanin:** tyrosine to L-DOPA to dopaquinone (tyrosinase or laccase), then non-enzymatic steps. It is a family census inside the graph, not a long pathway.
- **Control branch (decision for the user, section 9):** the sterol branch from farnesyl diphosphate, because the corpus already shows ergosterol-region hits (REPORTED).

### 6.4 Build steps for the graph (PROPOSED)
1. Seed enzymes: UniProt reviewed fungal entries by keyword and EC.
2. Reactions: Rhea by EC and ChEBI identifiers.
3. Compounds: ChEBI for structure, InChIKey, formula and mass; confirm with PubChem by InChIKey. Never join by KEGG's PubChem field.
4. KEGG: identifiers as cross-references; no content copied.
5. Downstream fungal steps: curated from papers with PMID, quoted evidence sentence and assay type.
6. Assisted extraction (optional): each edge needs the sentence, the PMID and a second independent check.
7. Classes: run NPClassifier on each compound's SMILES.

### 6.5 Exports for the comparison
- **Compound table:** `chebi_id, inchikey, formula, exact_mass, adducts` where adducts are computed for the modes the corpus uses (positive [M+H]+, [M+Na]+, [M+NH4]+, [M+2H]2+; negative modes as the EBC agent reports).
- **Capability table:** per taxon, per compound, the best path, its minimum tier and which steps the census found.
- Match on the first block of the InChIKey (connectivity layer) when stereochemistry is unknown (PROPOSED).

### 6.6 Comparison models (PROPOSED; each with controls and a null)
**M1: reaction mass shifts.** Each reaction edge has a formula difference. Test whether feature pairs within one acquisition group differ by those shifts more often than expected. Null: random shifts of matched size. Unit: group, then meta-analysis across groups. Group-local rule (corpus): never merge areas or row IDs across groups. Expected to be weak: networks in the corpus are sparse (REPORTED: 72% of features isolated in the case-study group).

**M2: capability against observation.** Response: compound observed at identification level 2 (library name) or by formula plus MS2 class. Predictor: the capability score from 6.5. Controls: negative taxa that cannot make the compound (candidate: *Saccharomyces*, *Komagataella*, which dominate the fungal groups; confirm absence of the pathway in their genomes) and positive taxa. Null: permute capability across taxa within study, instrument and polarity. Validation: hold out whole studies. Report AUC with intervals. Depends on the EBC counts and on detectability.

**M3: class-level anchors.** CANOPUS predicts classes for MS2 spectra; NPClassifier gives classes for graph compounds. Test enrichment of the matching classes in groups containing capable taxa against matched groups without. Null: permute class labels within instrument. Classes are coarse, but this covers features with MS2 and no library hit. Needs the SIRIUS login.

Order of dependency: M2 needs the EBC counts; M3 needs SIRIUS output; M1 needs only the graph and the feature tables.

### 6.7 Melanin-specific points
- DHN-melanin is a polymer. Metabolomics sees only monomers and intermediates.
- L-DOPA is not specific: plants and animals make it.
- Detectability must be measured before any model is run (request to the EBC agent).

### 6.8 Collaboration protocol
The rules in `EBC-KT/collab/README.md` apply: identify the sender, append-only files, counts and small tables only, label every number, no personal information, and no biology claimed from mass alone. Each agent pulls before it pushes. The BFD side signs as `BFD-HPCC-master agent`.

### 6.9 Review plan
1. An independent review (Fable) of this document.
2. A domain-expert review of the biochemistry (section 6.3 and the compounds), told to verify against PubMed and UniProt and cite each check.
3. The EBC agent reviews sections 6.5 to 6.7 through `collab/`.

## 7. Order of work (PROPOSED)
1. Stage 0: this amendment, three reviews, approval.
2. Stage 1: verified step sets and references for all three pathways, with decoys.
3. Stage 2: graph build (6.4).
4. Stage 3: census on real data for the verified step sets (the real-data plan, revised).
5. Stage 4: comparison models, using the EBC counts, SIRIUS output, and NPClassifier.
6. Stage 5: Chapters A (carotenoid) and B (melanin), with every number traced to a table.
SIRIUS login, the NPClassifier check and the EBC counts proceed in parallel with stages 0 and 1.

## 8. Risks
- Carotenoids may not be detectable in the corpus at all.
- Too few fungal groups with genomes on my side.
- Media and extraction protocol confound the fungal signal (REPORTED).
- Evidence tiers depend on curation quality. A wrong T1 corrupts every downstream claim.
- Language-model extraction can invent edges.
- NPClassifier may not run on current TensorFlow.
- Licence terms for a shared graph (KEGG; others UNVERIFIED).

## 9. Open items for the user
1. Include the sterol branch as a control? (My recommendation: yes, as a small branch.)
2. Which reviewers beyond the three named?
3. Whether the graph is to be released outside the group. The licence rule changes if so.

## 10. Left to the builder
Graph storage format and query layer, the exact capability score, adduct lists, tolerances, tree and HMM details from the real-data plan, and the order within stages.
