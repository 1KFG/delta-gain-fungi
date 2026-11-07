# Prioritizing Divergent Fungal Proteins for Industrial Applications

Nov 7, 2026 · @Jason Stajich

## Framing: novelty is one axis, usefulness is another

Delta-gain measures how much new protein space a genome adds. The industrially actionable question is where that novelty sits, which needs a prioritization layer on top of the novelty metric.

Divergent-but-domain-identifiable family members are valuable in three distinct situations, and the screening logic differs for each.

1. **Known function, new operating envelope.** The chemistry is already commercial; what is missing is a variant that works at the pH, temperature, solvent or substrate load the process demands. Prioritize by ecological provenance of the source genome, not sequence distance alone.
2. **Known function, altered selectivity.** Regioselectivity in P450s, enantioselectivity in hydrolases, chain-length specificity in lipases. The signal here is divergence at active-site-proximal positions while the catalytic scaffold stays conserved — the opposite of case 1.
3. **Unknown function worth discovering.** Orphan domain architectures and families whose fungal members have no characterized representative. Lower hit rate, higher upside, and the place where re-annotation adds knowledge rather than variants.

Score proteins separately for each case rather than producing one ranked list.

## A prioritization scheme: orthogonal scores, then tiered filters

Compute a small set of orthogonal scores per protein and filter on them rather than summing them. Summing collapses the distinction between "divergent everywhere" and "divergent in the right place", which is the information worth keeping.

| Score | What it captures | Why it reorders a candidate list |
| --- | --- | --- |
| Distance to nearest **characterized** homolog | Novelty relative to experimental knowledge, not to TrEMBL | A protein 40% identical to its nearest characterized relative outranks one 40% from its nearest database entry |
| Catalytic machinery occupancy | Whether known catalytic/binding residues survive in the alignment | Separates live candidates from pseudogenes, mis-predictions and neofunctionalized non-enzymes |
| Divergence partitioned by position class | Identity computed separately over catalytic columns, second-shell (\~8 Å) columns, and the rest | Case 2 wants divergence in the second shell with the site conserved; case 1 wants it everywhere else |
| Structural plausibility | pLM pseudo-perplexity as triage, then predicted structure and TM-score to the family fold | Catches bad gene models — the main risk when gene calls are your own |
| Expression plausibility | Secretion signal, absence of TM helices, length versus family distribution, intron-count outliers | Filters proteins that cannot be produced as secreted industrial enzymes |
| Genomic context | Cluster membership and tailoring complement for biosynthesis | A divergent core enzyme in a novel cluster beats the same enzyme in isolation |

Build the characterized-reference set from evidence-coded resources rather than from bulk sequence databases: Swiss-Prot `ACT_SITE`/`BINDING` features, [CAZy](http://www.cazy.org/) characterized entries, [BRENDA](https://www.brenda-enzymes.org/), [MIBiG](https://mibig.secondarymetabolites.org/) for biosynthetic gene clusters, and the PDB for site geometry. Mapping catalytic residues onto HMM alignment columns and checking occupancy is the cheapest high-yield filter in the whole scheme.

Then filter in tiers instead of ranking by a composite:

- **Tier 1** — site intact, good predicted fold, secreted where relevant, far from anything characterized. The drop-in-variant pool.
- **Tier 2** — same but with a partially divergent site. The selectivity-engineering pool.
- **Tier 3** — orphan architectures with no characterized family member. The discovery pool, kept separate because validation differs.

## Where divergent fungal enzymes are wanted, by how crowded the space is

Knowing which spaces are saturated matters before committing screening effort. The three rows marked "open" are the best match between divergent fungal family members and an unmet industrial need.

| Application space | Families to mine | How crowded | What actually differentiates a new variant |
| --- | --- | --- | --- |
| Lignocellulose degradation | GH5/6/7/10/11, AA9, AA14, AA16, accessory oxidoreductases | Heavily mined | Performance on real pretreatment hydrolysate, not model substrate; LPMO regioselectivity; electron-donor supply |
| Polyester and polyamide depolymerization | Cutinases, CE5, α/β-hydrolases, lipases | Open | Activity at or near polymer Tg; polyurethane and polyamide activity barely exists |
| Oxyfunctionalization | UPOs (heme-thiolate peroxygenases), P450s + cognate CPRs | Open for UPOs | Functional secreted expression; substrate scope beyond small aromatics; peroxide tolerance |
| Chitin and chitosan processing | GH18, GH20, GH46, CE4, AA11, AA15 | Moderately open | Defined product length and deacetylation pattern rather than bulk hydrolysis |
| Protein modification | Transglutaminases, prolyl endopeptidases, specific proteases | Open on the fungal side | Texturization of plant proteins; gluten processing; fungal space under-explored versus bacterial |
| Glycan remodeling | Fungal GTs, transglycosylating GHs | Open but hard | Defined oligosaccharide products; assay difficulty is why it stays open |
| Halogenation and tailoring | Flavin-dependent halogenases | Very thin in fungi | New families would be a real contribution for late-stage functionalization |
| Terpene synthases | Fungal TSs and their tailoring enzymes | Under-sampled | Product promiscuity blocks sequence-based assignment; expression screening is the only route |

Two cautions on the crowded row. Lignocellulose has absorbed decades of industrial effort, so a sequence variant alone rarely adds value; the remaining gaps are narrow, and [2024–2025 reviews](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2026.1743670/full) note that assay realism — testing in industrial matrices rather than on clean substrate — is what separates lab hits from process-ready catalysts. For secondary metabolism, prioritization at the cluster level is now the standard framing: the known/unknown quadrant scheme for clusters versus products, and tools such as [FunBGCeX](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11109268/), which deliberately targets biosynthetic enzymes lacking a detectable domain — directly relevant if your catalog contains domain-poor proteins you would otherwise discard.

## The known enzyme space for the three best-matched targets

These are the reference sets to anchor a distance-to-characterized score against. In all three the characterized space is small relative to the family size, which is exactly the condition that makes a divergence-based catalog useful.

### UPOs — small characterized set, large sequence space

The [UPObase](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8293615/) collection holds more than 1,900 putative UPO sequences, split into a short-type Family I of more than 1,100 members (mean \~30 kDa, mostly without predictable signal peptides) and a long-type Family II of more than 900 (mean 44.4 kDa, confined to Ascomycota and Basidiomycota, usually with a C-terminal disulfide bridge and a predicted signal peptide). Against that, roughly a dozen enzymes are actually characterized.

| Enzyme | Source fungus | Type | Why it matters as a reference |
| --- | --- | --- | --- |
| AaeUPO | *Agrocybe aegerita* | Long (Family II) | The model enzyme; evolved mutant libraries are the only engineered UPOs commercially available |
| CglUPO | *Chaetomium globosum* | Short (Family I) | Oxidizes testosterone, which AaeUPO and MroUPO cannot — the clearest published case of substrate scope tracking divergence |
| MroUPO | *Marasmius rotula* | Short (Family I) | Second most-studied; benchmark for the short-type scope |
| MweUPO | *Marasmius wettsteinii* | Short (Family I) | Close comparator to MroUPO |
| CciUPO | *Coprinopsis cinerea* | Subfamily I | Expressed in an industrial host, so an expression reference |
| MthUPO, TteUPO, HinUPO | *Myceliophthora thermophila*, *Thielavia terrestris*, *Humicola insolens* | Thermophile-derived | Provenance-based leads for the operating-envelope case |
| LfuCPO | *Leptoxyphium fumago* | Chloroperoxidase, short-type adjacent | One of only two commercially available members |

The gap to exploit: a [2019 genome-mining study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8293615/) searched more than 800 fungal genomes and recovered only 113 putative UPO sequences across 35 species, and a [2025 phylogeny of the peroxidase-peroxygenase superfamily](https://repository.tudelft.nl/record/uuid:4e78a3a7-94cd-4c19-afec-36c38368ae28) notes that most known sequences come from Dikarya while Mucoromycota, Glomeromycota and Chytridiomycota remain insufficiently investigated. Re-annotated genomes from early-diverging lineages are therefore high-value input.

Budget for the real bottleneck, not the sequence: published expression yields sit at mg-per-liter, purification recovery averages about 10%, unoptimized wild-type AaeUPO in yeast gave almost undetectable protein, and H₂O₂ both drives and inactivates the enzyme. Prioritize candidates with a predicted signal peptide and long-type architecture if expression is the limiting step.

### Polyester and polyamide hydrolases — one workhorse, thin coverage elsewhere

| Enzyme | Source fungus | Substrates reported | Structure |
| --- | --- | --- | --- |
| HiC (Novozym 51032) | *Humicola insolens* | PET, PBAT, Akestra, polyurethane adhesives; active on both ester and urethane bonds | [PDB 4OYY](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9603852/) |
| FsC | *Fusarium solani pisi* | PET | [PDB 1CEX](https://www.mdpi.com/2076-2607/10/6/1180) |
| FoCut5a | *Fusarium oxysporum* | PET | — |
| Cut3 | *Aspergillus niger* | PET and polycaprolactone surface modification | — |
| McCut | *Malbranchea cinnamomea* | Cutin/polyester | PDB 5X88 |
| Gc cutinase | *Glomerella cingulata* | Cutin/polyester | PDB 3DCN |
| TtcutA | *Thielavia terrestris* CAU709 | Polyurethane | — |
| CpCut1 | *Cladosporium* sp. P7 | Polyester-polyurethane; little PET activity | — |
| Assorted esterases | *Curvularia senegalensis*, *Aspergillus fumigatus* S45, *Pestalotiopsis microspora* E2712A, *Cladosporium pseudocladosporioides* T1.PL.1, *Alternaria* sp. PURDK2 | Polyurethane, including polyether PUR | — |

Two features of this table are the opportunity. First, HiC carries most of the applied literature essentially alone, and CpCut1 shows that a close homolog can be excellent on polyurethane while nearly inactive on PET — so sequence divergence within cutinases maps onto polymer specificity in a way nobody has systematically charted. Second, the bacterial side (LCC and its engineered ICCG variant, IsPETase) has absorbed most engineering effort, which leaves fungal diversity comparatively unsampled.

The quantitative bottleneck is crystallinity and glass-transition temperature, not turnover on soluble esters: HiC activity drops up to 25-fold as PET crystallinity rises from 7% to 35%, and its \~70 °C optimum matters because it approaches PET's Tg of about 72 °C. That makes thermophile-derived cutinases the single best-defined provenance filter you can apply. Polyamide is barely addressed in the fungal literature at all, and polyether polyurethane remains the hardest case — both are places where a modest measured activity would be novel.

### Chitin-active enzymes — a notable fungal gap in one family

| Family | Reference enzymes | Role | Fungal coverage |
| --- | --- | --- | --- |
| GH18 | Serratia marcescens ChiA/ChiB/ChiC (the mechanistic benchmark); Chit42 (*Trichoderma*) | Endo- and exo-chitinase; some members also transglycosylate | Large fungal family, few characterized members |
| GH19, GH20, GH46 | Plant/bacterial GH19; GH20 β-*N*-acetylhexosaminidases; GH46 chitosanases | Chain shortening, monomer release, chitosan hydrolysis | Sparse in fungi |
| AA10 | CBP21 (*Serratia marcescens*), the first chitin-active LPMO | Oxidative cleavage of crystalline chitin | Predominantly bacterial |
| AA11 | *Aspergillus oryzae* AA11 | Chitin-active LPMO | The fungal chitin LPMO family — and the clearest open slot in the table |
| AA15 | Insect-derived AA15 | Chitin-active LPMO | Non-fungal origin |
| CE4 | CDA from *Colletotrichum lindemuthianum*; *Absidia corymbifera*; *Pestalotiopsis* sp. | Deacetylation of chitin, chitosan and oligomers | Structural information remains limited |

The exploitable asymmetry: only AA10, AA11 and AA15 contain chitin-active LPMOs, AA10 is largely bacterial and AA15 is insect-derived, so AA11 is where fungal genomes are the only source and the characterized set is smallest. Pair that with CE4, where most characterized deacetylases show broad specificity and deacetylate several positions — which is precisely the problem, since the industrial need is for defined deacetylation patterns and defined oligomer lengths, not bulk conversion. A divergent CE4 with narrow positional specificity is a genuinely valuable object.

This target also has the best-defined combination play: GH18 plus an AA11 LPMO plus a CE4 deacetylase, drawn from one lineage so their pH and temperature optima already match, converting shell or insect-meal waste into defined chitooligosaccharides.

## Where the need is for combinations, not single enzymes

This is where a catalog beats a hit, because the unit of value is a matched set rather than a best-in-class enzyme.

- **Cocktail completion.** Industrial hydrolysis is limited by accessory activities, not the headline cellulase. Recovering the full accessory complement from one lineage — LPMO plus its electron donor, feruloyl and glucuronoyl esterases, swollenins — can beat mixing best-in-class enzymes from unrelated organisms.
- **Cascades with matched conditions.** Multi-step enzymatic routes fail on optimum mismatch. Organizing the catalog by source-organism physiology lets you assemble cascades whose members already share an envelope. Almost nobody uses this axis, and whole genomes rather than isolated sequences are what make it possible.
- **Cofactor and co-substrate partners.** Divergent NADPH-regenerating dehydrogenases, and catalases or oxidases tuned for controlled H₂O₂ supply, are perpetually needed and rarely the focus of discovery. For UPOs specifically, the peroxide-supply partner is as limiting as the peroxygenase.
- **Redox partners for P450s.** Prioritizing a divergent P450 without its cognate cytochrome P450 reductase wastes the hit; self-sufficient fusion architectures are worth flagging in the catalog as a separate class.
- **Whole-cluster transplantation.** For secondary metabolism the unit is the cluster. Shortlist clusters that are complete, novel in tailoring complement, and come from a host that is heterologously tractable or close to one.

## Making it defensible

Two things decide whether this is credible outside your group.

**Annotation-artifact control.** Because the input is genomes NCBI left unannotated, the first question any reviewer asks is whether the novel proteins are real. Build the controls in from the start: fold prediction on the shortlist, RNA-seq support where SRA data exist for the same or a close species, synteny and orthology consistency with annotated relatives, and explicit reporting of the fraction of "novel" proteins that fail each check. A reported failure rate is more persuasive than a clean-looking catalog.

**A funnel you can actually run.** The sequence-to-function gap is where discovery projects stall. A realistic shape:

1. 10⁵–10⁶ proteins in the catalog
2. computational tiering to \~10³
3. cell-free expression or small-scale yeast, *Pichia* or *Aspergillus* expression on 10²–10³ with a generic activity panel
4. full characterization on \~10

Decide the assay panel *before* the prioritization, because the assay determines which families are screenable at all. This inverts the usual order and saves a great deal of wasted ranking.

The highest-value publishable artifact is probably not the catalog but the prioritization method with its hit rate measured — showing that the scoring enriches for active, well-folded, non-artifactual enzymes relative to random draws from the same novel space. That turns delta-gain's novelty metric from descriptive into actionable, and it is the part others will cite and reuse.

## Three executable strategies on the delta-gain data

### 1. A FunBGCeX-style extension, aimed at the genomes it could not see

The decisive detail is in [FunBGCeX's](https://pmc.ncbi.nlm.nih.gov/articles/PMC11109268/) inputs: it takes a fungal genome **as a GenBank file with CDS features that include translations**, and the published run used 1,990 *annotated* NCBI fungal reference genomes. Delta-gain's input is precisely the set that requirement excludes. That is a clean, defensible niche rather than a reimplementation — and it is reinforced by the [observation](https://www.osti.gov/pages/servlets/purl/3025697) that many fungal assemblies reached NCBI without annotations, diminishing the usability of those datasets.

The method transfers almost directly, with your gene calls substituting for NCBI's:

| FunBGCeX step | What it does | Your adaptation |
| --- | --- | --- |
| Target detection | `hmmscan` against a target HMM, E-value 1e-5 | Same, but over your re-annotated proteomes; run the same HMM over annotated genomes too as a positive control |
| Neighborhood scan | Up to 20 flanking genes checked for core biosynthetic proteins against an in-house HMM library + DIAMOND | Same; your gene-model error rate is the new variable, so measure it on genomes where an NCBI annotation also exists |
| Clustering judgment | Biosynthetic character, duplication (self-resistance), relatedness to neighbors, size (<200 aa ignored unless biosynthetic), intergenic distance | Keep as-is; the <200 aa rule also protects you from fragmented gene calls |
| Known-BGC comparison | DIAMOND BLASTp vs FunBGCs proteins at ≥25% identity / ≥50% coverage; clustered if hits share a parent BGC origin within 15,000 bp; stop after three uncorrelated genes | Use FunBGCs directly ([Zenodo 10.5281/zenodo.8126803](https://pmc.ncbi.nlm.nih.gov/articles/PMC11109268/)) — no need to rebuild it |
| Novelty score | Average identity × proportion of homologous genes relative to a known BGC (homologous = ≥45% identity, ≥50% coverage) | This is the component that plugs straight into delta-gain: it is a per-cluster novelty metric with the same shape as yours |

The part most worth copying is how they handled domainless proteins, because it is the direct answer to your divergent-members problem. Of 5,070 FunBGCs proteins, Pfam found 520 domains; the remaining 572 domainless proteins were grouped by all-vs-all BLASTp (E-value 0.01), aligned with Clustal Omega, and converted into **59 new HMM profiles** with `hmmbuild`. Do the same with your divergent clusters: a group of family members too far from Pfam to annotate confidently becomes a new HMM, which then becomes a detection tool. That inverts the usual framing — the proteins that look like annotation noise become the instrument.

Scale expectation from their run: \~1,990 genomes gave 1,050 Pyr4-containing BGCs, of which 690 lacked a second core gene, and manual inspection narrowed to 16 SHC/OSC-encoding and 7 UbiA-like candidates, with four taken into *Aspergillus oryzae* NSAR1. Note also that antiSMASH recovered only 325 Pyr4 BGCs from the same genomes — the gain came from the custom HMM library, not from more genomes. So the honest claim for a delta-gain version is "new clusters from genomes no pipeline has looked at, plus new HMMs from proteins no pipeline could classify."

### 2. Mining additional UPOs

This is the most tractable of the three, because the search model is already published and the target is a single compact family.

1. Build the query set from UPObase's long- and short-type profiles rather than one seed sequence: Family I (>1,100 members, \~30 kDa, mostly no predictable signal peptide) and Family II (>900 members, 44.4 kDa, Ascomycota and Basidiomycota only, usually a C-terminal disulfide bridge plus a signal peptide). Two HMMs, not one.
2. Search your re-annotated proteomes with both, plus profile-profile or structure-based search (the heme-thiolate fold is well conserved even where sequence is not) to catch members the HMMs miss.
3. Score each hit on the catalytic machinery: the heme-thiolate cysteine ligand and the conserved acid-base glutamate are the pass/fail test. A "UPO" lacking the proximal cysteine is not one.
4. Stratify by expression feasibility — predicted signal peptide, long-type architecture, disulfide pattern — because this, not sequence novelty, is what killed most previous candidates.
5. Benchmark distance against the roughly twelve genuinely characterized enzymes (AaeUPO, CglUPO, MroUPO, MweUPO, CciUPO and the thermophile-derived set), not against the 1,900 putative sequences.

The quantitative headroom is good: the published genome-mining effort recovered just 113 putative UPOs from more than 800 genomes across 35 species, so the family is sparsely represented per genome and the catalog is nowhere near saturation.

### 3. A targeted early-diverging deliverable

This is the strongest-positioned product of the three, because the gap is explicitly stated in the literature: most known UPO sequences come from Dikarya, and Mucoromycota, Glomeromycota and Chytridiomycota remain insufficiently investigated. The same under-sampling holds for BGCs.

Be realistic about scale before promising a deliverable. Published phylogenomic efforts worked with tens of genomes per phylum — one study assembled 131 zygomycete genomes (97 Mucoromycota, 34 Zoopagomycota), and for most species only a single representative genome exists; Chytridiomycota has long been worse, with very few annotated genomes released. So step one is a census, not a search: query NCBI for assemblies in these three phyla, split by whether an annotation exists, and report that split. That number is itself a result, and it sizes everything downstream.

What makes the deliverable credible:

- **Phylum-stratified reporting.** Report novelty per phylum rather than pooled, so the early-diverging contribution is visible instead of swamped by Dikarya.
- **Two families, not a catalog dump.** UPOs and one chitin-active family (AA11 or CE4) are the right first pair: both have small characterized reference sets, and chitin is mechanistically plausible in Chytridiomycota and Mucoromycota given their cell-wall biology.
- **Honest caveats on Glomeromycotina.** Large genomes (>100 Mb) with 20–30k genes, obligate biotrophy, and multinucleate spores carrying different genotypes mean gene-model confidence will be lowest exactly where novelty looks highest. Report per-phylum artifact rates, or a reviewer will raise it first.
- **A concrete output.** Per phylum: the census, the HMM hits with catalytic-residue status, predicted structures for the shortlist, and the new HMMs built from divergent groups — deposited, so the instrument is reusable even before any enzyme is assayed.

## Open question

The answer changes the plan substantially: is there access to an expression-and-assay capability, even a modest one?

- **Yes** — go at UPOs (AA11-adjacent expression caveats included), fungal cutinases against polyamide and polyether polyurethane, and the GH18 + AA11 + CE4 chitin set. All three have small characterized reference sets and a defined assay.
- **No, computational for now** — the stronger play is the method paper plus cluster-level and cocktail-level prioritization, where the claims rest on genomic context rather than measured activity.

## Sources

- [Identification and Expression of New Unspecific Peroxygenases — Recent Advances, Challenges and Opportunities](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8293615/) (Kinner, Rosenthal & Lütz, *Front. Bioeng. Biotechnol.* 2021) — UPO families, UPObase counts, genome-mining numbers, expression limitations
- [Molecular diversity of unspecific heme peroxygenases with promising applications in sustainable oxyfunctionalizations](https://repository.tudelft.nl/record/uuid:4e78a3a7-94cd-4c19-afec-36c38368ae28) (Zámocký et al. 2025) — superfamily phylogeny; under-investigated early-diverging lineages
- [Surfing the wave of oxyfunctionalization chemistry by engineering fungal unspecific peroxygenases](https://research.tudelft.nl/en/publications/surfing-the-wave-of-oxyfunctionalization-chemistry-by-engineering) (Beltrán-Nogal et al., *Curr. Opin. Struct. Biol.* 2022) — long and short UPO family engineering
- [Challenges and perspectives in using unspecific peroxygenases for organic synthesis](https://www.frontiersin.org/journals/catalysis/articles/10.3389/fctls.2024.1470616/pdf) (Huang et al. 2024)
- [Fungal Enzymes Involved in Plastics Biodegradation](https://www.mdpi.com/2076-2607/10/6/1180) (*Microorganisms* 2022) — fungal cutinases, lipases and esterases on PET and PUR; crystallinity effects
- [An Overview into PET Hydrolases and Efforts in Tailoring Enzymes for Improved Plastic Degradation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9603852/) — HiC, PDB codes, degradation conditions
- [Identification and characterization of a fungal cutinase-like enzyme CpCut1 from *Cladosporium* sp. P7 for polyurethane degradation](https://www.researchgate.net/publication/378773374_Identification_and_characterization_of_a_fungal_cutinase-like_enzyme_CpCut1_from_Cladosporium_sp_P7_for_polyurethane_degradation) — the PET-inactive, PUR-active divergence case
- [Comparison of the enzymatic depolymerization of PET and Akestra using *Humicola insolens* cutinase](https://www.frontiersin.org/journals/chemical-engineering/articles/10.3389/fceng.2022.1048744/full) — HiC temperature optimum versus polymer Tg
- [Enzymatic Modifications of Chitin, Chitosan, and Chitooligosaccharides](https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2019.00243/full) (*Front. Bioeng. Biotechnol.* 2019) — AA10/AA11/AA15, CE4 deacetylases, CBP21
- [Current Perspectives on Chitinolytic Enzymes and Their Agro-Industrial Applications](https://www.mdpi.com/2079-7737/10/12/1319) — GH18/19/23/48 classification and mechanism
- [Structure and function of a CE4 deacetylase isolated from a marine environment](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5673215/) — CE4 substrate promiscuity and the limited structural coverage
- [Discovery of fungal onoceroid triterpenoids through domainless enzyme-targeted global genome mining](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11109268/) (*Nat. Commun.* 2024) — FunBGCeX; mining clusters whose enzymes lack detectable domains
- [Deciphering Chemical Logic of Fungal Natural Product Biosynthesis through Heterologous Expression and Genome Mining](https://pmc.ncbi.nlm.nih.gov/articles/PMC9906657/) — the known/unknown 2×2 cluster prioritization framing
- [Precision enzyme discovery through targeted mining of metagenomic data](https://link.springer.com/article/10.1007/s13659-023-00426-8) (*Nat. Prod. Bioprospect.* 2024) — structure-informed targeted mining for process-relevant properties
- [Fungal biotechnology review](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2026.1743670/full) (*Front. Microbiol.* 2026) — ML-guided BGC prioritization and assay realism in industrial matrices
- [Divergent Evolution of Early Terrestrial Fungi Reveals the Evolution of Mucormycosis Pathogenicity Factors](https://academic.oup.com/gbe/article/15/4/evad046/7079964) (*Genome Biol. Evol.* 2023) — 131 zygomycete genomes; Glomeromycotina genome sizes and gene counts
- [Sequencing the Genomes of the First Terrestrial Fungal Lineages: What Have We Learned?](https://pmc.ncbi.nlm.nih.gov/articles/PMC10386755/) — uneven sampling across Mucoromycota and Zoopagomycota; AMF sequencing obstacles
- [Fungal genome annotation availability analysis](https://www.osti.gov/pages/servlets/purl/3025697) — assemblies deposited to NCBI without annotations

