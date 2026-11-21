# Review: pathway census design (carotenoid and melanin)

Reviewed 2026-11-20. Reviewer: Claude (Fable), independent read.
Design reviewed: `2026-11-20-pathway-census-design.md` (not edited).
Also read: `analysis/enzymes_UPO/DESIGN.md`, `RESULT_T2_negative.md`, `STRUCTURE_VS_RESIDUES.md`,
`upo_place/README.md`, and the BFD layout note.

Labels in this review:
- **VERIFIED**: I ran a query or read a file or a PubMed abstract today. The command or source is given.
- **RECALLED**: from my memory of the literature. Not checked today. Treat as a lead, not a fact.
- **OPINION**: my judgement.

Literature checks used PubMed. Each cited abstract has its DOI link.

## Verdict

The design has sound anti-goals and a sound non-goal (no activity prediction).
It is not ready to build. Four problems can each produce confident wrong results:

1. The BFD-effect test does not measure what question 4 asks.
2. The melanin step lists contain errors, and the L-DOPA "pathway" is not a pathway that a completeness score can describe.
3. No species tree was found. Goal G4 and pipeline steps 6 and 7 depend on one.
4. Step assignment by tree placement is not shown to work for the large families. The UPO precedent does not transfer without new validation.

OPINION: build the carotenoid chapter first as a pilot. Decide on the melanin chapter after that.

---

## Findings, ranked

### 1. The BFD-effect comparison is confounded and tests the wrong quantity (highest importance)

**VERIFIED facts** (DuckDB on `Fungi_BFD_runs/tables/samples.parquet`, `busco_genome.parquet`, joined to `prior_annotation_status.tsv` on `folder = ASMID`):

| Group | Genomes | Median BUSCO complete | Fraction below 80 % BUSCO | Species |
|---|---|---|---|---|
| `baseline_ok` | 5,625 | 98.8 | 0.13 | 2,791 |
| `no_prior_annotation_present` | 18,055 | 97.7 | 0.21 | 5,502 |
| no status row | 673 | no BUSCO row | - | 204 |

- The two groups differ in assembly quality before any biology is measured (13 % versus 21 % of genomes below 80 % BUSCO).
- Only 835 of 7,568 species have genomes in both groups. Only 58 of 196 orders have at least 5 genomes in each group. 25 orders have no baseline genome. 56 orders have no BFD-only genome.
- Outside Ascomycota, Basidiomycota, Mucoromycota and Microsporidia, one of the two groups is very small. Examples: Kickxellomycota 187 baseline and 3 BFD-only; Chytridiomycota 82 and 8; Mortierellomycota 85 and 15.
- The baseline is not "the NCBI pipeline". The `annotation_provider` column has 669 distinct providers for `baseline_ok` rows. The largest is DOE JGI (1,287 rows). The design's risk statement ("The funannotate and NCBI pipelines predict different gene models") describes a two-pipeline contrast that the data do not have.
- The largest genera are population-scale sets: *Fusarium* 2,174, *Aspergillus* 1,632, *Saccharomyces* 1,566, *Candidozyma* 1,240 genomes.

**OPINION, the main point.** Question 4 mixes two different effects. The design must separate them.

- **Annotation effect.** Does BFD's gene calling change the step calls for a genome? The correct test is paired: the same assembly, deposited proteins versus BFD proteins. The 5,625 baseline genomes allow this test directly. Deposited proteins are at `1KFG/2026/NCBI_fungi/source/NCBI_ASM/<acc>/protein.faa` (path from the layout note; I did not open these files today). This test has no lineage, quality or sampling confound. The design does not include it.
- **Sampling effect.** Do the 18,055 extra genomes show pathway states (gains, losses, expansions) that no baseline genome shows? This is a question about new lineages and new states. A difference in means between the two groups does not answer it.

The design's step 8 and section 5.6 compare the BFD matrix between the two groups. If the matrix uses BFD gene models for both groups (G2 reads this way; the design does not say so explicitly), then a group difference reflects which genomes submitters chose to annotate. That choice is not random. It correlates with assembly quality (verified above), with taxon, and probably with year and sequencing project (not checked). A "finding" under section 5.6 can therefore be a submitter-selection effect.

Further gaps:
- Section 5.6 needs "survival of the genome-quality and downsampling controls". No threshold, model or test is named. A reader cannot tell in advance what result fails.
- "Effect size and confidence interval" needs a stated resampling unit. Genomes are not independent. A confidence interval from genome-level resampling is too narrow. Use species or genus as the resampling unit.
- The number of tests is steps x lineages x (presence, copy number) x pathways. No multiple-testing rule is given.

**Recommended change.** State two estimands. (a) Paired annotation effect on baseline genomes: per step, count calls gained, lost and unchanged between deposited and BFD proteins. (b) Sampling effect: per step, list pathway states that occur only in clades without a baseline genome, with tblastn or miniprot confirmation. Drop the unpaired group contrast, or keep it as a descriptive table only.

### 2. Melanin step lists: errors, and one pathway that a completeness score cannot describe

The design marks these lists UNVERIFIED. I checked some items. The rest is recall.

**VERIFIED from PubMed abstracts:**

- **Rsy1 is scytalone dehydratase. It is not a separate step.** The design lists "scytalone dehydratase (Scd1/Arp1)" and then "Rsy1" as another item. An abstract on *Pyricularia oryzae* names the enzyme "scytalone dehydratase (SDH1/RSY1)" ([DOI](https://doi.org/10.1021/acs.jafc.1c04984)). Listing both makes one enzyme count as two steps.
- **Ayg1 is not a general DHN-melanin step.** In *Aspergillus fumigatus*, Ayg1 shortens the heptaketide YWA1 to 1,3,6,8-tetrahydroxynaphthalene (T4HN). The abstract states that no enzyme with similar activity was reported in another organism at that time, 2004 ([DOI](https://doi.org/10.1074/jbc.M406758200)). The *Colletotrichum lagenarium* PKS1 makes T4HN directly ([DOI](https://doi.org/10.1271/bbb.63.1445)). So a genome with a PKS1-type synthase and no Ayg1 can have a complete pathway. A fixed step list will call it incomplete.
- **Not all Ascomycota, and not all *Aspergillus*, use DHN-melanin.** *Aspergillus terreus* lacks the DHN-melanin genes. It makes Asp-melanin with an NRPS-like enzyme (MelA) and a tyrosinase (TyrP) ([DOI](https://doi.org/10.1016/j.chembiol.2016.03.014)). The "expected pathway per lineage" idea fails inside one genus. See finding 3.
- The Ayg1 catalytic triad is Ser257, Asp352, His380 ([DOI](https://doi.org/10.1074/jbc.M406758200)). This is a usable motif check for that step.

**RECALLED, not checked today:**

- The "WA type" PKS label is risky. *Aspergillus nidulans* wA makes the naphthopyrone YWA1, and I recall that the *A. nidulans* conidial pigment is not DHN-melanin. The same non-reducing PKS clade also makes other products in other fungi. Clade membership of the PKS does not prove melanin synthesis.
- Abr1 is a ferroxidase-like multicopper oxidase and Abr2 is a laccase. They are two different proteins. The design groups them as one step.
- Some fungi have a second naphthol reductase (4HNR) in addition to the THN reductase (Buf1/Arp2/Thr1 type).
- Tyrosinase catalyses both tyrosine hydroxylation and L-DOPA oxidation. This answers the design's open question: a separate hydroxylase is not needed where a tyrosinase is present.
- *Cryptococcus neoformans* melanin synthesis uses a laccase (Lac1) on exogenous diphenol substrates. I recall that it has no tyrosinase route. My PubMed query for this returned no records, so this stays unverified.
- Tyrosinase-based melanin is not limited to Basidiomycota. *Neurospora crassa* tyrosinase and *Aspergillus oryzae* MelO are classic Ascomycota examples.
- Other fungal melanin routes exist and are absent from the design: pyomelanin (homogentisate), GHB-melanin in *Agaricus*, catechol melanin.

**OPINION.** The L-DOPA "pathway" has one or two enzymes and no pathway-specific steps. "Completeness" has no meaning for it. The melanin role of a laccase or tyrosinase cannot be inferred by orthology in families this large (VERIFIED sizes: Tyrosinase PF00264 has 165,751 proteins in 15,682 genomes; the three Cu-oxidase Pfams have about 218,000 to 229,000 proteins each in about 23,100 genomes). A tree clade that contains Lac1 is a clade of laccases. It is not a clade of melanin laccases. Report this sub-chapter as a family census of tyrosinases and laccases. Do not report it as a pathway census.

### 3. "Expected pathway per lineage" is circular

OPINION. The anti-goal says completeness must be judged against the lineage's expected pathway. The expected set is unknown for most of the 196 orders. The census is the measurement that would define it. If the builder sets expectations at phylum level from a few model species, the output repeats the input. The *A. terreus* case (VERIFIED, finding 2) shows that the expectation fails below genus level.

Recommended change: remove "expected lineages" from the definition file as a scoring input. Report observed step sets per clade. Keep the model-organism pathway only as a labelled reference row.

### 4. No species tree was found

**VERIFIED searches:**

- `DeltaGain_Fungi` (depth 6): the only tree files are the UPO gene trees under `analysis/enzymes_UPO/data/refset_zamocky2022/`.
- `Fungi_BFD_runs/analysis`: only transcription-factor ortholog gene trees (`T20_ortho_*_tree.nwk`).
- `Fungi_BFD/README_phylogeny.md` documents a PHYling pipeline. Its output goes to `phylogeny/<taxon>/.../tree/concat.treefile`. Neither `Fungi_BFD/phylogeny` nor `Fungi_BFD_runs/phylogeny` exists. No `phylo.csv` exists in either folder.
- `Fungi_BFD/docs/BFD_Introduction.md` has no mention of a species tree (grep for "tree" and "phylog").
- The UPO project also lists the species tree as an input. I found no species tree used in its files.

**Limit of this check.** Some `find` commands timed out in `Fungi_BFD_runs` (deep folders), `DeltaGain/results`, `DeltaGain/data` and `DeltaGain_Fungi/runs`. I did not search `1KFG/2026`. The correct statement is "not found in the places listed", not "does not exist".

**Consequence.** The non-goal "A new species tree. The tree is an input" leaves G4 without its input. Step 6 also needs "close relatives", which needs a tree or a substitute.

**OPINION, options:**
- Build a backbone tree with the existing PHYling pipeline on one high-BUSCO genome per genus or per species. This is new work. Its cost is not measured.
- Use the taxonomy ranks in `samples.parquet` as nested random effects (genus in family in order). This needs no tree. It does not support ancestral-state reconstruction.
- Use sister-clade or within-genus contrasts.
- `Fungi_BFD_runs/results/mash_sketch` and `results/ANI` exist (VERIFIED as folders; contents not inspected). They could support species-level dereplication and a "nearest relative" definition.

In all options, reduce to one genome per species (7,568 species, VERIFIED) before any comparative test.

### 5. Step assignment at scale is not shown to work

**VERIFIED family sizes** (`pfam.parquet`, distinct proteins and genomes; I did not check the score threshold of the table rows):

| Pfam | Proteins | Genomes | Step that uses it |
|---|---|---|---|
| PF00106 adh_short | 2,601,444 | 23,335 | THN reductase |
| PF00067 p450 | 2,021,480 | 23,184 | astaxanthin synthase |
| PF00171 Aldedh | 512,554 | 23,203 | ylo-1 type dehydrogenase |
| PF00109 ketoacyl-synt | 373,940 | 23,149 | PKS |
| PF01593 Amino_oxidase | 320,805 | 23,159 | phytoene desaturase |
| PF00394 Cu-oxidase | 217,678 | 23,100 | laccase, Abr1/Abr2 |
| PF00264 Tyrosinase | 165,751 | 15,682 | tyrosinase |
| PF00348 polyprenyl_synt | 113,160 | 23,210 | GGPP synthase |
| PF03055 RPE65 | 61,220 | 15,315 | carotenoid cleavage oxygenase |
| PF00494 SQS_PSY | 56,923 | 23,064 | phytoene synthase |
| PF02982 Scytalone_dh | 19,540 | 11,676 | scytalone dehydratase |
| PF18916 Lycopene_cyc | 1,032 | 1,026 | lycopene cyclase |

- The UPO placement tool worked on one family of about 97,520 proteins. It had a 187-leaf curated reference tree and a leave-one-out validation (184/187 correct, from `upo_place/README.md`). This design has 10 or more families. Several are 3 to 27 times larger. No reference tree exists for any of them.
- PF00494 is present in 23,064 genomes at about 2.5 proteins per genome. The Pfam name is "SQS_PSY". RECALLED: squalene synthase (ERG9) is in this family and is present in almost all fungi. PF00348 is present in 23,210 genomes at about 4.9 proteins per genome. RECALLED: GGPP synthase is needed for protein prenylation and for diterpene clusters, so it is not a carotenoid-specific gene. Both steps will be "present" in almost every genome if the call is loose. The carotenoid completeness score will be inflated.
- PF18916 has only 1,032 proteins in 1,026 genomes. PF00494 is in 23,064 genomes. I do not know if this Pfam model covers the fungal cyclase domain of the fused al-2/carRA/crtYB proteins. If it does not, the planned "start from the Pfam hits" prefilter will miss this step. Test this on the known reference proteins before the pipeline is built.

**OPINION, what the design needs:**
- A reference tree for each step must include non-step paralogs as decoy clades. If the tree has only step references, placement puts every candidate inside a step clade. The `unassigned` class then never triggers.
- Repeat the leave-one-out and matched-identity validation of `upo_place` for each step. A step that fails validation gets a family-level label only.
- For the THN reductase (SDR family) and the polymerising multicopper oxidase, sequence orthology is unlikely to be enough. Use physical linkage to the PKS (same scaffold, set window) as the evidence, and label the call "candidate by context".
- `swissprot_best.parquet` exists (VERIFIED: 157,393,345 rows, best SwissProt hit per protein). The design does not list it. It is a cheap prefilter and a cross-check for step calls.
- The step-HMM score cutoff needs calibration against known paralogs. An HMM built from a few references has no gathering threshold.

### 6. Carotenoid step list: smaller errors

**VERIFIED from PubMed abstracts:**
- "Carotenoid cleavage oxygenase (cao-2/carX family)" mixes two enzymes. In *Fusarium fujikuroi*, CarX makes retinal from beta-carotene, and deletion of carX did not stop neurosporaxanthin synthesis ([DOI](https://doi.org/10.1128/EC.00392-06)). CarT is the torulene-cleaving enzyme of the neurosporaxanthin pathway ([DOI](https://doi.org/10.1111/j.1365-2958.2007.05665.x)). These must be two steps with separate reference clades. Both are in the RPE65 family.
- The *Fusarium* bifunctional gene carRA encodes phytoene synthase and carotene cyclase, and carB encodes the desaturase (same two abstracts). The design's statement on the fused gene is correct. The design assigns "carRA" only loosely to "several lineages"; RECALLED: the Mucorales names are carRA (*Phycomyces*) and carRP (*Mucor*).

**RECALLED, not checked today:**
- al-3 = GGPP synthase, al-2 = fused cyclase and phytoene synthase, al-1 = desaturase, ylo-1 = aldehyde dehydrogenase, cao-2 = the *Neurospora* CarT ortholog. These agree with the design.
- "crtY" is the bacterial lycopene cyclase name. I recall no characterised standalone fungal cyclase of that type. A standalone crtY-type hit in a fungal assembly is more probably contamination than novelty. See finding 8.
- Astaxanthin synthase (crtS) needs a cytochrome P450 reductase partner (crtR). The design does not list it. It is a general P450 reductase, so it is a weak pathway marker.
- Mucorales cleave beta-carotene for trisporic acid synthesis. The design does not list this branch.

### 7. Fused genes and copy counting

OPINION. "A fused gene must count toward both steps" is correct. It is not enough.
- Count and place **domains**, not proteins. Place the cyclase domain and the synthase domain on separate trees.
- A gene-model error can split a fused gene into two models, or merge it with a neighbour. Carotenoid genes are often adjacent. The length and truncation filter (step 4) will then remove true genes or keep merged ones. Record split and merged cases as their own class.
- Count genes, not transcripts. Protein IDs have a `-T1` suffix. I did not check if any genome has more than one transcript per gene.
- Assemblies of dikaryons or diploids can carry both haplotypes. This doubles the copy number. `busco_genome.parquet` has a `duplicated_pct` column (VERIFIED). Use it as a covariate or a flag. The design uses only completeness and N50.

### 8. Contamination control is missing

VERIFIED: the UPO design requires contamination screening before any lichen or Glomeromycota result. This design has no such control. The design itself names Lecanoromycetes as heavily sampled.

OPINION. This matters more here than for UPOs. Carotenoid genes (crtE, crtB, crtI, crtY) are common in bacteria and in algal photobionts. A contaminant contig in a lichen assembly can produce a false "gain", a false "standalone cyclase" and a high novelty rank in section 5.5 ("nearest characterised reference being far"). Add the same checks as the UPO design: neighbour-gene taxonomy on the contig, GC content, and a flag for contigs that carry only the candidate.

### 9. The tblastn absence check

OPINION.
- Scope is not bounded. "Genomes that lack a step that their close relatives have" can be thousands of genomes times many steps. Measure the count from the matrix first. Size the jobs after that.
- For large families, tblastn always finds a hit (a paralog). A hit does not show that the step gene is present. The recovered sequence must be placed on the step tree, the same as an annotated protein.
- tblastn fragments hits at introns. A spliced aligner (for example miniprot) gives a model that can be placed. I did not check if miniprot is installed on HPCC.
- The check is one-sided. It tests absences only. It never tests false presences (contamination, pseudogenes, merged models). The matrix is then biased toward presence.
- Three outcomes ("present but unannotated, fragmented, or truly absent") need written rules. A fourth outcome is needed: "assembly gap, cannot decide".
- The trigger needs relatives. A BFD-only lineage with no close relative is never checked.

### 10. The prioritisation scheme

OPINION. The scheme is defensible as a **triage list**. It is not defensible as a **ranking**.

- Good: no fitted weights, terms set before results, one column per term, hard filters apart from score terms. This matches the UPO lesson (VERIFIED: clade carried most of the signal; pocket partial r = +0.129 given clade; structure did not transfer out of clade).
- Problem 1: the terms pull in opposite directions. "Nearest reference is far" and "altered active-site residue" mark the most novel copies. They also mark the copies least likely to perform the step, and the copies most likely to be paralogs with another function, pseudogenes or contaminants. One ordered list hides this conflict.
- Problem 2: no outcome data exist, so no order can be validated. A lexicographic order or an additive score implies a precision that is not there. Rank 3 versus rank 30 has no meaning.
- Problem 3: "expressibility" needs a host and a goal. The design names neither. Transmembrane helices are a defect for a secreted enzyme and normal for a membrane-bound desaturase or P450.
- Problem 4: "altered active-site residue" needs known catalytic residues for each step. These are not in the design. Today I verified only the Ayg1 triad.
- Problem 5: near-identical copies from 2,174 *Fusarium* genomes will fill the list. Dereplicate candidates by sequence cluster before listing.

Recommended change: output named bins, not ranks. For example: (A) canonical ortholog in an untested lineage, (B) lineage-specific expansion member, (C) divergent copy with intact motifs, (D) copy in a conserved gene cluster. Give each bin its stated purpose and its main risk. Let the testing partner choose bins. State in the output that the order inside a bin has no meaning.

### 11. Scope and housekeeping

- **Stale numbers marked as measured.** VERIFIED today: `prior_annotation_status.tsv` has 5,703 `baseline_ok` rows, not 5,511. It also has 18,145 `no_prior_annotation_present`, 2,359 `no_prior_annotation_not_downloaded`, 1,184 `annotated_NOT_downloaded` and 10 `annotated_folder_missing_protein`. `samples.parquet` has 24,363 genomes, not 23,680. 673 of them have no status row. The design asks for a recount; the recount changes the numbers.
- **The table list is stale.** VERIFIED: `tables/` now has `tcdb.parquet`, `deeploc.parquet`, `tf_inventory.parquet`, `swissprot_best.parquet` and `gene_intergenic_distances.parquet`. The design copies the gaps from a 13-day-old note and labels them MEASURED.
- **Over-scope (OPINION).** One generic tool, three pathways, a BFD-effect test and a candidate list in one version is too much. The claim "about 70 % duplicated code" for option A is not measured. Build carotenoid first. Its families are the smallest and its steps are the best defined. Make the tool generic after one pathway works.
- **Section 7 gives the builder the decisions that decide validity**: the comparative method, the orthology rule and the score form. The UPO design required the decision rule to be written before coding. Apply the same rule here.
- **Section 8 is empty.** The open items for the user should include at least: the estimand for question 4, the species-tree source, and if L-DOPA melanin stays in scope as a pathway.
- **Copy-number tests need a count model** and a statement of the unit (species, not genome). The design names neither.

---

## What I checked and what I did not

Checked today (VERIFIED):
- Status counts, provider counts, group sizes, BUSCO by group, phylum, order and species overlap.
- Pfam family sizes for 12 step-related families.
- Presence or absence of tree files in the folders listed in finding 4.
- Six PubMed abstracts: RSY1 = scytalone dehydratase; Ayg1 function and triad; PKS1 product; CarX; CarT and carRA/carB; Asp-melanin in *A. terreus*.

Not checked (stated as RECALLED or left open):
- All *Neurospora*, *Phycomyces*, *Mucor*, *Xanthophyllomyces* and *Rhodotorula* gene assignments.
- Abr1/Abr2 identity, the wA product in *A. nidulans*, the *Cryptococcus* laccase route, tyrosinase bifunctionality.
- Whether PF18916 covers the fungal cyclase domain.
- Whether deposited protein files are usable for the paired test.
- Whether a species tree exists outside the searched folders.
- Any runtime or job-size estimate. I have no data for these.
