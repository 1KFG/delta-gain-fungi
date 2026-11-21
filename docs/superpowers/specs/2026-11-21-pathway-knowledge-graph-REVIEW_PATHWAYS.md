# Domain review: pathway step lists (carotenoid, DHN-melanin, L-DOPA melanin)

Review of `2026-11-21-pathway-knowledge-graph-design.md` sections 4, 6.2 and 6.3, of `pathways/carotenoid.yaml`, and of `refs/carotenoid/UNVERIFIED.md`.
Date: 2026-11-21. Reviewer: domain-review agent (fungal secondary metabolism). No file other than this one was changed.
LC-MS detectability is out of scope (another agent covers it).

## 0. How to read the labels

- **V-abs**: VERIFIED. I read the PubMed abstract (or, for PMID 16416328, the full text file `ojima_flow.txt`) on 2026-11-21.
- **V-db**: VERIFIED. I read the database record on 2026-11-21 (UniProt REST, Rhea REST, KEGG REST, NCBI protein).
- **V-title**: Partly verified. PubMed returned the record and I read the title only. The content is not checked.
- **RECALLED**: From memory. Not checked. Do not enter into the graph without a check.

Literature statements are based on articles retrieved from PubMed. DOI links are in section 8.
I did not read any full text except Ojima et al. 2006. A V-abs statement is only as strong as the abstract.

## 1. Main findings

1. **DHN-melanin, first reductase (6.3) is wrong as written.** THR1 of *Colletotrichum lagenarium* is the 1,3,8-trihydroxynaphthalene (T3HN) reductase, not the T4HN-to-scytalone reductase. There are two reductase paralogues with different substrate preference. See E1.
2. **"Ayg1 belongs only to the heptaketide route" is wrong.** The Ayg1 orthologue WdYg1 acts in the hexaketide route of *Wangiella dermatitidis*. See E2.
3. **"PKS1 type to 1,3,6,8-THN" hides three different routes** (direct THN synthase; hexaketide plus Yg1; heptaketide plus Ayg1). The census needs three PKS reference clades and one Ayg1/Yg1 step. See E3.
4. **L-DOPA melanin: a laccase does not convert tyrosine to L-DOPA.** *Cryptococcus neoformans* makes melanin only from exogenous diphenols. Tyrosinase and laccase must be separate steps. See E6, E7.
5. **carotenoid.yaml, phytoene desaturase EC.** One EC (1.3.99.30) is used for all references. The *Phycomyces* enzyme does four desaturations (lycopene-forming), the *Neurospora* enzyme does five. See E8.
6. **Astaxanthin synthase accession gap is solved.** NCBI CAC16261 is 100% identical (557/557) to UniProt Q3HR17. B5LEC8 differs at residue 60 (D in B5LEC8, G in CAC16261 and Q3HR17). See section 2.3.
7. **Torulene to torularhodin: I found no paper that identifies the enzyme.** Two PubMed queries gave no primary paper. The step must be T4.
8. **The Rhea example in section 3 (RHEA:56772) is not on the fungal route.** Adonixanthin is not an intermediate reported by Ojima et al. See E10.

## 2. Task 1: check of each stated fact

### 2.1 Carotenoid: yaml steps

| Step (yaml) | Statement | Result | Source |
|---|---|---|---|
| ggpp_synthase | P24322 = al-3, *N. crassa*; Q9P885 = carG, *Mucor circinelloides*; EC 2.5.1.29 | Correct. Both entries are reviewed. UniProt lists EC 2.5.1.29 with 2.5.1.1 and others. | V-db (UniProt) |
| phytoene_synthase / lycopene_cyclase | P37295 al-2, Q8X0Z1 carRA, Q9P854 carRA, Q9UUQ6 CARRP, Q7Z859 crtYB; EC 2.5.1.32 and 5.5.1.19; fused protein | Correct. All five are reviewed and bifunctional. | V-db (UniProt); Rhea RHEA:34475, RHEA:32219, RHEA:32239 |
| phytoene_desaturase | EC 1.3.99.30 for all four refs | **Error E8** (below). | V-db, V-abs |
| phytoene_desaturase | O13506 crtI, *X. dendrorhous* | Entry exists, unreviewed, no EC. PMID 10397832 not checked by me. | V-db |
| carx | S0ENQ2, EC 1.13.11.63, retinal-forming, not needed for neurosporaxanthin | Correct. CarX also cleaves gamma-carotene, torulene and beta-apo-8'-carotenal to retinal in vitro. | V-db; V-abs PMID 17293483 |
| cart | A1KQY3 carT, A7UXI1 cao-2, EC 1.13.11.59; torulene to beta-apo-4'-carotenal | Correct. Both enzymes were assayed in vitro. CAO-2 did not cleave gamma-carotene. Rhea: torulene + O2 = 4'-apo-beta-carotenal + 3-methyl-2-butenal (RHEA:31519). | V-db; V-abs PMID 17493127, 17610084 |
| ylo1 | Q1K615 ylo-1, F6IBC7 carD, EC 1.2.1.82 | Correct. Both proteins made neurosporaxanthin from beta-apo-4'-carotenal in vitro. See order flag O2. | V-db; V-abs PMID 18627463, 21749649; RHEA:31515 |
| astaxanthin_synthase | B5LEC8, PMID 18769767; `ec: null`; crtR as electron donor (PMID 18837978) | Correct. B5LEC8 is unreviewed and has no EC. PMID 18769767 reports cloning and deletion of crtS. crtR deletion stops astaxanthin and gives beta-carotene. KEGG K23037 lists EC 1.14.99.63, 1.14.15.24 and 1.14.99.-. I recommend that `ec` stays null and the KEGG ECs are cross-references only (RECALLED: those EC entries were defined for non-P450 enzymes). | V-db; V-abs PMID 18769767, 18837978 |

I did not re-check the PMIDs 1826006, 12695851, 8163509, 12172798, 11172012, 10951210, 10589832, 11862485, 10397832 in the yaml.

### 2.2 Carotenoid: errors and corrections

**E8. Phytoene desaturase EC (yaml).** The yaml gives one EC, 1.3.99.30, "from UniProt P54982".
- EC 1.3.99.30 = 15-cis-phytoene + 5 A = 3,4-didehydrolycopene (RHEA:30975). EC 1.3.99.31 = 15-cis-phytoene + 4 A = all-trans-lycopene (RHEA:15585). V-db.
- *Neurospora* AL-1 does up to five desaturations. Products are 3,4-didehydrolycopene and lycopene. V-abs PMID 11017770. UniProt P21334 and Q8X0Z0 carry EC 1.3.99.- (V-db).
- *Phycomyces* carB does "the four dehydrogenations that convert phytoene to lycopene". V-abs PMID 9079885. UniProt P54982 still carries EC 1.3.99.30 (V-db). The UniProt EC and the paper do not agree.
- Correction: do not give the step one EC. Record two reactions (four-step, lycopene-forming; five-step, 3,4-didehydrolycopene-forming). Attach the five-step reaction as T1 only to AL-1. The number of steps cannot be read from the sequence family. This decides whether a lineage can make torulene and neurosporaxanthin.

**E9. UniProt keyword as seed (section 4 and 6.4 step 1).** The count is correct: 35 reviewed fungal entries with KW-0125 (V-db, same query, today). The list contains *Saccharomyces cerevisiae* BTS1 (Q12051), *Schizosaccharomyces pombe* spo9 (O59703) and *Eremothecium gossypii* BTS1 (Q758K0). These are GGPP synthases. RECALLED: these yeasts make no carotenoids. Many other entries (for example *Aspergillus oryzae* Q2U4X9, *Trichophyton verrucosum* D4D802) are named by similarity; I did not see a paper for them. The keyword is not evidence that a species has the pathway. Use it to find candidates, then tier each entry by its own paper.

**E10. Rhea example (section 3).** RHEA:56772 exists and the equation is quoted correctly (V-db). Rhea assigns it EC 1.14.99.64. Ojima et al. report this route in *X. dendrorhous*: beta-carotene to echinenone, then canthaxanthin or 3-hydroxyechinenone, then phoenicoxanthin, then astaxanthin (V-abs/full text, PMID 16416328). Adonixanthin is not in that route. A Rhea search for "phoenicoxanthin" returned no reaction (V-db). Rhea has beta-carotene to echinenone (RHEA:55660) and echinenone to canthaxanthin (RHEA:55664), both under EC 1.14.99.63. The fungal 3-hydroxylation reactions need curated reaction ids.

**E11. Gene-name clash (anti-goal "family similarity read as function").** `carS` means two different things. In *Phycomyces*, carS is the beta-carotene cleavage oxygenase (V-abs PMID 22445959, 21854466). In *Fusarium fujikuroi*, carS is a regulator that represses carotenoid synthesis (V-title PMID 33383912, 29371556). The UNVERIFIED.md item 2 writes "carS and others" for the Mucorales step; that is correct for *Phycomyces* only. Key every gene name by organism.

**E12. cao-1 is a decoy.** *Neurospora* CAO-1 is in the same family as CAO-2. It did not convert carotenoids. It cleaves resveratrol and piceatannol (V-abs PMID 23893079; UniProt Q7S860 reviewed). Add it to the decoy set of the PF03055 tree.

### 2.3 Astaxanthin synthase (section 4 and UNVERIFIED.md item 8)

| Statement in the design | Result | Source |
|---|---|---|
| Gene is 3,166 bp, accession AX034666 | Correct. The paper also states 17 introns. | V (full text, lines 17, 352-354) |
| cDNA accession AX034665 | Correct. | V (full text, lines 404-405) |
| Protein is 557 residues | Correct. | V (full text, line 363) |
| CAC16261 is the protein of AX034665 | The paper text does not contain "CAC16261". I did not open the AX034665 nucleotide record. NCBI CAC16261.1 is "unnamed protein product [Phaffia rhodozyma]", 557 aa. | V-db for the protein record; the link to AX034665 is not checked by me |
| B5LEC8 differs from CAC16261 at one of 557 residues | Correct. Position 60: G in CAC16261, D in B5LEC8. | V-db (my alignment of the two FASTA records) |
| A single P450 does 4-ketolation and 3-hydroxylation | Supported by genetics in *X. dendrorhous*. In *E. coli* with a yeast P450 reductase the enzyme made canthaxanthin only (about 1% conversion); no hydroxylated product was found. The 3-hydroxylase activity is therefore shown in the fungus only. | V (full text, lines 470-500, 673-690) |
| No reviewed UniProt entry for astaxanthin synthase | Correct. B5LEC8, Q3HR17 and A0A0A8KXB4 are all unreviewed. | V-db |

New fact: **UniProt Q3HR17 (gene Ast / crtS) is identical to CAC16261 over all 557 residues** (V-db, my alignment). Q3HR17 can fill the `accession: null` reference. UniProt links Q3HR17 to PMID 25504221, not to PMID 16416328; the link to the Ojima paper rests on sequence identity with the patent sequence.

### 2.4 Section 4: other facts

| Statement | Result | Source |
|---|---|---|
| KEGG M01042 = K15745, K17841, K17842, K17819 | Correct. The module is "lycopene => neurosporaxanthin". KEGG lists it complete in *Neurospora*, *Sordaria*, *Podospora*. | V-db (KEGG) |
| K23037 CRTS, ASY | Correct. | V-db |
| KEGG compounds C00779, C04033, C01624, C00355, C00822, C01693, C05606, C17936 | All eight identifiers and names are correct. | V-db |
| K17740 SCD1 (EC 4.2.1.94), K00505 TYR, K00421 and K05909 laccase | Correct. | V-db |
| "I did not find a polyketide synthase KO for the DHN pathway" | KEGG has K19580 `rppA`, 1,3,6,8-tetrahydroxynaphthalene synthase, EC 2.3.1.233. RECALLED: rppA is a bacterial type III PKS, not the fungal type I enzyme. Do not use it for fungi. KEGG also has **K17739 THNR, tetrahydroxynaphthalene reductase, EC 1.1.1.252**, which the design does not list. | V-db for the KO records |
| KEGG has no trisporic acid or torularhodin step (section 5) | KEGG compound search for "trisporic" and for "torularhodin" returned nothing. KO search for "trisporic" returned nothing. | V-db |
| Exact masses (7 compounds) | All agree with my calculation from the formulas (C40H56 536.4382; C40H54 534.4226; C40H54O 550.4175; C40H52O2 564.3967; C40H52O4 596.3866; C40H56O2 568.4280; C20H28O 284.2140). I did not check the formula of each compound in ChEBI. | Calculated |
| map00906 counts (128, 149, 61, 44); ChEBI/PubChem InChIKey | Not checked. | - |

### 2.5 DHN-melanin (6.3): errors and corrections

**E1. First reductase.** The design says "a reductase (Arp2/THR1 type) to scytalone".
- *C. lagenarium* THR1 converts 1,3,8-trihydroxynaphthalene to vermelone. V-abs PMID 8672814. UniProt P87025 names it "Trihydroxynaphthalene reductase (T3HN reductase)" (V-db).
- Thr1p can also use T4HN in vivo, and a second, T4HN-specific reductase is deduced in *C. lagenarium*. V-abs PMID 12967213.
- *Magnaporthe grisea* has two enzymes. Tetrahydroxynaphthalene reductase (4HNR; UniProt Q12634, EC 1.1.1.252) prefers T4HN by a factor of 310. Trihydroxynaphthalene reductase (3HNR) prefers T3HN by a factor of 4.2. The two are 46% identical. V-abs PMID 10956664. The 3HNR cloning paper is PMID 8112349 (V-title).
- *Aspergillus fumigatus* Arp2 is called 1,3,6,8-THN reductase; arp2 deletants accumulate flaviolin only. V-abs PMID 11350964, 10515939. UniProt E9QUT3 (V-db).
- Correction: T4HN to scytalone is the 4HNR/Arp2 type (EC 1.1.1.252; RHEA:21908). T3HN to vermelone is the 3HNR/THR1/BUF1/BRM2 type. The two types are paralogues in one short-chain dehydrogenase family. They need separate reference clades. Substrate overlap is real, so the census must not call one step absent when only one paralogue is present.

**E2. Ayg1.** The design says "Ayg1 belongs only to the heptaketide route".
- *A. fumigatus*: Alb1 makes a heptaketide; Ayg1 shortens it to 1,3,6,8-THN. V-abs PMID 11350964.
- *W. dermatitidis*: WdYg1, an Ayg1 orthologue, is needed to convert the **hexaketide** 2-acetyl-1,3,6,8-THN to T4HN. V-abs PMID 18676950.
- *Botrytis cinerea*: BcYGH1 converts the BcPKS13 product for conidial melanin. V-abs PMID 26514268.
- *A. niger*: OlvA is the Ayg1 homologue. V-abs PMID 21277986.
- Correction: Ayg1/Yg1 is a chain-shortening hydrolase. It is needed where the PKS releases a hexaketide or heptaketide. It is not needed where the PKS is a true THN synthase.

**E3. PKS type.** "PKS1 type" covers three routes:
- True THN synthase: *C. lagenarium* PKS1 gives T4HN in *A. oryzae* (V-abs PMID 10501004). Its thioesterase domain cyclises and deacetylates an enzyme-bound hexaketide (V-abs PMID 23261597). UniProt P79068 (unreviewed).
- Hexaketide plus Yg1: *W. dermatitidis* WdPKS1 (V-abs PMID 11179356, 18676950). UniProt Q9Y7A7 (unreviewed).
- Heptaketide (naphthopyrone YWA1) plus Ayg1: *A. fumigatus* Alb1/PksP, UniProt Q4WZA8 "Naphthopyrone synthase" (V-db; V-abs PMID 11350964).
- Decoy: *A. nidulans* wA is a naphthopyrone synthase in the same PKS clade (UniProt Q03149; V-title PMID 11251292). RECALLED: *A. nidulans* does not make DHN-melanin from it; the green pigment needs the laccase yA (UniProt P17489, V-db).
- All these are non-reducing type I iterative PKS. The product chain length cannot be read from a domain hit.

**E4. "Rsy1 is scytalone dehydratase".** RECALLED as correct (RSY1, "rosy", is the *Magnaporthe* scytalone dehydratase locus). UniProt P56221 gives only the gene name SDH1 (V-db). I could not confirm the RSY1 synonym in a record today.

**E5. Second dehydration.** The design says "Scd1 again to 1,8-dihydroxynaphthalene". Rhea has only the scytalone reaction for EC 4.2.1.94 (RHEA:24396); a Rhea search for "vermelone" returned no reaction (V-db). RECALLED: scytalone dehydratase is reported to dehydrate vermelone in *Magnaporthe*. RECALLED: for *A. fumigatus*, some pathway schemes put Abr1 at the vermelone step. I did not check either. Treat vermelone to 1,8-DHN as a curated reaction with its own evidence record, not as a copy of the scytalone edge.

**Correct statements (DHN).**
- Scytalone dehydratase: SCD1 *C. lagenarium* (UniProt Q00455, V-db; V-abs PMID 8953707); Arp1 *A. fumigatus* (O14434; V-abs PMID 9383199, 10515939: the arp1 deletant accumulates scytalone and flaviolin); SDH1 *M. oryzae* (P56221).
- *A. terreus* lacks the genes for DHN-melanin. It makes Asp-melanin: MelA (NRPS-like) makes aspulvinone E and the tyrosinase TyrP activates it for polymerisation. V-abs PMID 27133313. UniProt Q0CRX1, Q0CRX0 (V-db). Rhea has the aspulvinone E synthesis reaction RHEA:63824 (V-db).
- Laccase or multicopper oxidase: *A. fumigatus* has both abr1 (multicopper oxidase) and abr2 (laccase) in the cluster (V-abs PMID 10515939). Deletion of abr2 changes conidia from grey-green to brown (V-abs PMID 16988817). The abstracts show phenotype, not the polymerisation reaction.

### 2.6 L-DOPA melanin (6.3): errors and corrections

**E6. "tyrosine to L-DOPA to dopaquinone (tyrosinase or laccase)".**
- Tyrosinase, EC 1.14.18.1: Rhea gives L-tyrosine + O2 = L-dopaquinone + H2O (RHEA:18117) and 2 L-dopa + O2 = 2 L-dopaquinone + 2 H2O (RHEA:34287). V-db. Free L-DOPA is not a required intermediate in the Rhea reaction.
- Laccase, EC 1.10.3.2, oxidises diphenols (RHEA:11276, V-db). RECALLED: laccase does not hydroxylate tyrosine. The "tyrosine to L-DOPA" edge must not carry a laccase.

**E7. *Cryptococcus*.** The *C. neoformans* enzyme is a laccase (CNLAC1; UniProt Q55P57, J9VY90; EC 1.10.3.2). It oxidised dopamine to decarboxy-dopachrome. V-abs PMID 8300520. *C. neoformans* "produces melanin only from exogenous substrates" and needs diphenolic precursors. V-abs PMID 17098915. A second laccase, LAC2, exists (UniProt J9VQZ4, V-db). So in *Cryptococcus* the pathway is one enzyme plus a substrate from the environment. The compound that enters the graph depends on the substrate supplied (L-DOPA, dopamine, homogentisate and others).

**Section 6.7 "L-DOPA is not specific".** Agreed. Add: in *Cryptococcus* the L-DOPA is not made by the fungus (V-abs PMID 17098915).

### 2.7 Section 6.2 (evidence tiers)

Section 6.2 states no facts. Two problems for use:
- T1 says "reaction shown in a fungus". Many fungal steps are shown by gene deletion and phenotype only (abr1, abr2, arp2, Sordaria, Podospora). The schema field "assay type" (6.1) must be used with the tier. I use T1a = enzyme assay or heterologous product identification with the fungal protein; T1b = gene deletion or complementation with intermediate analysis in the fungus; T1c = deletion with colour phenotype only.
- T1 must be bound to the organism and the protein. A T1 for AL-1 (five steps) is not a T1 for carB of *Phycomyces* (four steps).

## 3. Task 2: missing steps and enzymes

### 3.1 Carotenoid: upstream and core

| Step | Gene / enzyme | Organism | Source | Label |
|---|---|---|---|---|
| IPP isomerase, GGPP synthase, structural genes of the astaxanthin route | idi, crtE, crtYB, crtI, crtS (cloned; crtYB, crtI, crtS deleted) | *X. dendrorhous* | PMID 18769767 | V-abs |
| Core genes | carRA, carB | *Blakeslea trispora* | PMID 15345447; UniProt Q67GH9 (carRA, reviewed) | V-title; V-db |
| Torulene formation order | AL-1 does not accept gamma-carotene. Torulene forms by cyclisation of 3,4-didehydrolycopene. | *N. crassa* | PMID 11017770 | V-abs |
| Regulator, not a step | carS (repressor) | *F. fujikuroi* | PMID 33383912, 29371556 | V-title |
| Opsin, not a step | carO in the car cluster | *F. fujikuroi* | PMID 17293483 | V-abs |

### 3.2 Carotenoid cleavage and apocarotenoids

| Step | Gene / enzyme | Organism | Source | Label |
|---|---|---|---|---|
| beta-carotene to 2 retinal | CarX (S0ENQ2) | *F. fujikuroi* | PMID 17293483 | V-abs (in vitro) |
| beta-carotene to 2 retinal | Cco1. Deletion raises beta-carotene. Three opsin genes in the genome. | *Ustilago maydis* | PMID 19584000 | V-abs (in vitro and in vivo) |
| torulene to beta-apo-4'-carotenal | CarT; CAO-2 | *F. fujikuroi*; *N. crassa* | PMID 17493127; 17610084 | V-abs |
| beta-apo-4'-carotenal to neurosporaxanthin | YLO-1; CarD | *N. crassa*; *F. fujikuroi* | PMID 18627463; 21749649 | V-abs |
| Decoy: resveratrol cleavage | CAO-1 (Q7S860) | *N. crassa* | PMID 23893079 | V-abs |

### 3.3 Mucorales trisporic acid pathway

| Step | Gene / enzyme | Organism | Source | Label |
|---|---|---|---|---|
| beta-carotene to beta-apo-12'-carotenal (C25), cleavage at C11'-C12' | CarS (UniProt G8Z9M3, unreviewed) | *Phycomyces blakesleeanus* | PMID 21854466; 22445959 | V-abs (product in *E. coli*) |
| beta-apo-12'-carotenal to beta-apo-13-carotenone (C18), cleavage at C13-C14 | AcaA (UniProt A0ABR3AZC2, unreviewed) | *P. blakesleeanus* | PMID 21854466 | V-abs (in vitro) |
| Carotene cleavage, induced by trisporic acid | tsp3 (and tsp4 in *Rhizopus oryzae*); UniProt A5HXP5 "Putative carotene oxygenase" | *B. trispora*, *R. oryzae* | PMID 17822929 | V-abs. Function shown as loss of carotene in *E. coli*; the abstract names no product. |
| 4-dihydrotrisporin to trisporin (NADP-dependent, short-chain dehydrogenase) | TSP2 (UniProt B9DR41, unreviewed) | *Mucor mucedo* | PMID 18931040 | V-abs (enzyme purified) |
| 4-dihydromethyltrisporate to methyltrisporate (NADP-dependent, aldo/keto reductase) | tdh (UniProt Q01213, reviewed, EC 1.1.1.-) | *M. mucedo* | PMID 8828234 | V-abs (enzyme purified) |
| Review of the pathway | - | Mucorales | PMID 19665150 | V-title |

- The *Phycomyces* paper states that the genome has five carotenoid-oxygenase-like sequences and that the function of two is known. Three are open. V-abs PMID 21854466.
- RECALLED: the steps between the C18 ketone and trisporic acid (oxidations at C4 and of the C1 methyl group, and the methyl ester step) have no identified gene beyond TSP2 and tdh. I did not check this.
- I did not find the name "tsp1" in PubMed. RECALLED: TSP1 is the gene name of the tdh product. Not checked.
- Rhea search for "trisporate" and "12'-apo-beta-carotenal" returned no reaction (V-db). All these reactions need curated ids.

### 3.4 *Rhodotorula* and other red yeasts

- Torulene to torularhodin (oxidation of a methyl group to a carboxyl group): **no enzyme found.** PubMed queries "torularhodin[Title] AND (gene OR enzyme OR synthase OR oxidase OR pathway)" and a second query on torulene oxidation returned no primary paper on the enzyme. A review exists (PMID 35430937, V-title; not read). Rhea and KEGG have no torularhodin record (V-db). RECALLED: the enzyme is unknown and a mixed-function oxidase is proposed.
- RECALLED: red yeasts have crtYB-type and crtI-type genes. My PubMed query for a functional paper returned nothing useful. Not checked.
- UNVERIFIED.md item 7 is right that there is no Pucciniomycotina reference. Without one, the fused synthase/cyclase and desaturase trees have no anchor in that subphylum.

### 3.5 Carotenoid glycosylation and esterification in fungi

- One PubMed query ("carotenoid glycoside fungus glycosyltransferase") returned bacterial and plant enzymes and one engineered yeast paper (PMID 42915850, V-title: engineered *Yarrowia lipolytica*). It returned no native fungal carotenoid glycosyltransferase.
- A query for astaxanthin esters or glucosides in *Xanthophyllomyces* returned nothing relevant.
- RECALLED: no native fungal carotenoid glycosyltransferase or acyltransferase gene is characterised. RECALLED: *X. dendrorhous* astaxanthin is mainly not esterified. Neither is checked. Recommendation: do not add a glycosylation step. If one is added it is T4.

### 3.6 Astaxanthin route in *Xanthophyllomyces*

- crtS (asy): beta-carotene to astaxanthin through echinenone, canthaxanthin or 3-hydroxyechinenone, and phoenicoxanthin. Side route: torulene to 4-ketotorulene to 3-hydroxy-4-ketotorulene. V (full text PMID 16416328).
- crtR: P450 reductase, UniProt B6D970 (unreviewed, EC 1.6.2.4). Deletion gives beta-carotene. V-abs PMID 18837978. The yaml leaves crtR out as a weak marker. I agree for the census. In the graph it must be a required partner of the reaction, or the edge is incomplete.

### 3.7 DHN-melanin: genes by organism

| Organism | Genes | Source | Label |
|---|---|---|---|
| *Aspergillus fumigatus* (conidia) | alb1/pksP (Q4WZA8), ayg1 (Q4WZB3), arp2 (E9QUT3), arp1 (O14434), abr1 (Q4WZB4), abr2 (E9RBR0); six-gene cluster, 19 kb | PMID 10515939, 11350964, 9383199, 16988817 | V-abs; V-db |
| *Aspergillus niger* (conidia) | fwnA (PKS), olvA (ayg1 homologue), brnA (abr1 homologue), pptA. Naphtho-gamma-pyrones depend on fwnA. | PMID 21277986 | V-abs |
| *Colletotrichum lagenarium* | PKS1 (P79068), SCD1 (Q00455), THR1 (P87025), plus a deduced T4HN-specific reductase | PMID 10501004, 23261597, 8953707, 8672814, 12967213 | V-abs |
| *Magnaporthe grisea / oryzae* | scytalone dehydratase SDH1 (P56221); 4HNR (Q12634); 3HNR | PMID 10956664; 8112349 | V-abs; V-title |
| *Wangiella (Exophiala) dermatitidis* | WdPKS1 (Q9Y7A7), WdYG1 (Q6BCB5) | PMID 11179356, 18676950 | V-abs |
| *Botrytis cinerea* | bcpks12 (sclerotia), bcpks13 (conidia), bcygh1, bcbrn1, bcbrn2, bcscd1; three pathway transcription factors | PMID 26514268 | V-abs |
| *Sordaria macrospora* | pks, teh, sdh, tih | PMID 17681008 | V-abs |
| *Alternaria alternata* | ALM, BRM1, BRM2 cluster | PMID 8392512 | V-title |
| *Podospora anserina* | PaPKS1 | PMID 17707627, 33974074 | V-title |
| *Bipolaris oryzae* | SCD1 | PMID 14716498 | V-title |
| *Sclerotinia sclerotiorum* | SCD1, THR1 | PMID 29024255 | V-title |
| *Verticillium dahliae* | VdPKS1; transcription factor VdCmr1 | PMID 28755294, 35637443, 29485393 | V-title |
| *Cochliobolus heterostrophus* | melanin depends on the transcription factor Cmr1 | PMID 17237364 | V-title |
| *Pestalotiopsis fici* | PfmaG, tetrahydroxynaphthalene reductase (UniProt W3XAA0, reviewed, EC 1.1.1.252) | UniProt | V-db |
| *Fonsecaea monophora* | transcriptome study of melanin pathways | PMID 27504027 | V-title (abstract cut off) |
| *Sporothrix* spp. | DHN-melanin in conidia and yeast cells; L-DOPA raises melanin; pigment from L-tyrosine in 72 of 73 strains | PMID 23042177 | V-abs |

- Laccase for DHN polymerisation in *Colletotrichum* and *Magnaporthe*: my PubMed query for a *C. lagenarium* LAC1 paper returned nothing. RECALLED only. A *W. dermatitidis* laccase query also returned nothing.
- Transcription factors (Cmr1, Pig1 type) are regulators. Keep them out of the step list.

### 3.8 Other fungal melanin routes

| Route | Enzymes | Organism | Source | Label |
|---|---|---|---|---|
| Laccase melanin from exogenous diphenols | LAC1, LAC2 | *Cryptococcus neoformans* | PMID 8300520, 11500433, 17098915 | V-abs |
| Pyomelanin from L-tyrosine through homogentisate | hppD (Q4WHU1, EC 1.13.11.27), hmgA (Q4WHT9, EC 1.13.11.5). hmgA deletion raises pigment; hppD deletion removes it. | *A. fumigatus* | PMID 19028908 | V-abs |
| Pyomelanin-like pigment from homogentisate | laccase-dependent use of exogenous HGA | *C. neoformans* | PMID 17098915 | V-abs |
| Asp-melanin | melA (NRPS-like, aspulvinone E), tyrP (tyrosinase) | *A. terreus* | PMID 27133313 | V-abs |
| Tyrosinase (DOPA route) | T (P00440, EC 1.14.18.1) | *N. crassa* | UniProt; PMID 2529259 | V-db; V-title |
| Tyrosinase | melB | *A. oryzae* | PMID 16233650, 22213164 | V-title |
| GHB-melanin (gamma-L-glutaminyl-4-hydroxybenzene, oxidised by tyrosinase) | tyrosinase (PPO) | *Agaricus bisporus* | PMID 10552675 | V-title. RECALLED: GHB is the spore melanin precursor. |
| Teliospore pigment | laccase and polyketide synthase genes | *Ustilago maydis* | PMID 25226432 | V-title |
| Catechol melanin in Basidiomycota | - | - | - | RECALLED. Not checked. |

- Pyomelanin: in the graph, HppD makes homogentisate; the pigment forms when HmgA does not remove it. Presence of hppD and hmgA is the tyrosine degradation pathway. It is present in most fungi (RECALLED). The census cannot tell pigment production from catabolism.
- The *A. fumigatus* conidial extract also shows "laccase-type activity" on L-DOPA (V-abs PMID 14970241). This is an example of the laccase substrate range. It is not evidence of a DOPA pathway.

## 4. Task 3: lineage-restricted steps (absence is not a defect)

| Step | Restriction | Evidence | Label |
|---|---|---|---|
| Whole carotenoid pathway | Absent in *Saccharomyces* and related yeasts although UniProt tags their GGPP synthase with the carotenoid keyword. | UniProt Q12051, O59703, Q758K0 | V-db for the entries; absence of the pathway is RECALLED |
| carT/cao-2 and ylo-1/carD (neurosporaxanthin) | KEGG lists module M01042 as complete in *Neurospora*, *Sordaria*, *Podospora*. Papers show the genes in *Neurospora* and *Fusarium*. Mucorales, *Ustilago* and *Xanthophyllomyces* end at beta-carotene, retinal or astaxanthin. | KEGG M01042; PMID 17610084, 17493127, 19584000 | V-db, V-abs. The full taxon range is not measured. |
| Fifth desaturation (3,4-didehydrolycopene) | Shown for *Neurospora* AL-1. *Phycomyces* carB does four. | PMID 11017770, 9079885 | V-abs |
| crtS astaxanthin synthase | Shown only in *X. dendrorhous*. | PMID 16416328, 18769767 | V-abs. "Only" is RECALLED; I did not search other taxa. |
| carS, acaA, tsp3, TSP2, tdh (trisporoids) | Mucorales sexual signalling. | PMID 21854466, 22445959, 17822929, 18931040, 8828234 | V-abs |
| CarX / Cco1 (retinal) | Linked to opsins. *U. maydis* makes beta-carotene mainly for retinal. | PMID 17293483, 19584000 | V-abs |
| Torularhodin step | Red yeasts (Pucciniomycotina). Enzyme unknown. | PMID 35430937 (title) | RECALLED |
| Ayg1 / Yg1 | Needed only with hexaketide or heptaketide PKS. Not needed with the *C. lagenarium* THN synthase. | PMID 23261597, 11350964, 18676950 | V-abs |
| Second PKS for melanin | *B. cinerea* and other Leotiomycetes have two PKS (sclerotia, conidia). | PMID 26514268 | V-abs |
| DHN-melanin as a whole | *A. terreus* lacks the genes; relatives have them. | PMID 27133313 | V-abs |
| Two reductases (4HNR and 3HNR) | Shown in *M. grisea* and *C. lagenarium*. The *A. fumigatus* cluster has one (arp2). | PMID 10956664, 12967213, 10515939 | V-abs. Whether *A. fumigatus* has a second reductase outside the cluster is not checked. |
| abr1 (multicopper oxidase) next to abr2 (laccase) | Described for the *A. fumigatus* cluster; *A. niger* brnA is the abr1 homologue. | PMID 10515939, 21277986 | V-abs. Range in other lineages not checked. |
| Tyrosine to L-DOPA in *Cryptococcus* | Not present; substrate is exogenous. | PMID 17098915 | V-abs |
| Asp-melanin (melA, tyrP) | *A. terreus*. | PMID 27133313 | V-abs. Range not checked. |

Consequence for the capability score: a "missing step" must be judged against the route that the lineage uses. A single linear DHN step list will mark *C. lagenarium* as lacking Ayg1 and *A. terreus* as defective. Both calls would be wrong.

## 5. Task 4: evidence tiers

Tier codes: T1a = fungal enzyme assayed or product identified by heterologous expression. T1b = fungal gene deletion or complementation with intermediate analysis. T1c = fungal deletion, colour phenotype only. T2, T3, T4 as in section 6.2. The tier holds for the named organism and protein only. Orthologues in other species are T3 until a paper shows the reaction.

### 5.1 Carotenoid

| Reaction | Enzyme (organism) | Tier | Source | Label |
|---|---|---|---|---|
| FPP + IPP to GGPP | al-3 (*N. crassa*), carG (*M. circinelloides*) | T1 (assay type not checked) | UniProt P24322, Q9P885, reviewed | V-db; papers not read |
| 2 GGPP to phytoene | al-2, carRA, carRP, crtYB | T1 (assay type not checked by me) | UniProt reviewed entries; PMID 18769767 for crtYB deletion and *E. coli* complementation of a truncated crtYB | V-db; V-abs for crtYB |
| phytoene to lycopene (4 steps) | carB (*Phycomyces*) | T1a | PMID 9079885 | V-abs |
| phytoene to 3,4-didehydrolycopene (5 steps) | AL-1 (*N. crassa*) | T1a | PMID 11017770 | V-abs |
| 5 steps in *Fusarium* carB, *X. dendrorhous* crtI | - | T3 from my checks | UniProt names only | V-db |
| lycopene to gamma-carotene to beta-carotene; 3,4-didehydrolycopene to torulene | al-2 / carRA / crtYB cyclase domain | T1 (PMID 11862485 in the yaml; not read by me). Loss of cyclase activity in truncated crtYB: PMID 18769767 | V-abs for crtYB only |
| beta-carotene to retinal | CarX (*F. fujikuroi*); Cco1 (*U. maydis*) | T1a | PMID 17293483; 19584000 | V-abs |
| torulene to beta-apo-4'-carotenal | CarT; CAO-2 | T1a | PMID 17493127; 17610084 | V-abs |
| beta-apo-4'-carotenal to neurosporaxanthin | YLO-1; CarD | T1a | PMID 18627463; 21749649 | V-abs |
| beta-carotene to echinenone to canthaxanthin (4-ketolation) | CrtS + P450 reductase | T1a | PMID 16416328 (canthaxanthin in *E. coli*) | V full text |
| 3-hydroxylation steps to phoenicoxanthin and astaxanthin | CrtS (*X. dendrorhous*) | T1b | PMID 16416328, 18769767 (genetics in the fungus; no hydroxylated product in *E. coli*) | V full text; V-abs |
| Electron supply to CrtS | CrtR | T1b | PMID 18837978 | V-abs |
| torulene to 4-ketotorulene to 3-hydroxy-4-ketotorulene | CrtS | T1b | PMID 16416328 | V full text |
| beta-carotene to beta-apo-12'-carotenal | CarS (*Phycomyces*) | T1a | PMID 21854466, 22445959 | V-abs |
| beta-apo-12'-carotenal to beta-apo-13-carotenone | AcaA (*Phycomyces*) | T1a | PMID 21854466 | V-abs |
| carotene cleavage (product not named in abstract) | tsp3 (*B. trispora*, *R. oryzae*) | T1a for substrate loss; product T3 | PMID 17822929 | V-abs |
| 4-dihydrotrisporin oxidation | TSP2 (*M. mucedo*) | T1a | PMID 18931040 | V-abs |
| 4-dihydromethyltrisporate oxidation | tdh (*M. mucedo*) | T1a | PMID 8828234 | V-abs |
| Other steps from the C18 ketone to trisporic acid | no gene known to me | T4 | - | RECALLED |
| torulene to torularhodin | unknown | T4 | no paper found | search result |
| Carotenoid glycosylation or esterification | none known to me | T4 | no paper found | search result |

### 5.2 DHN-melanin

| Reaction | Enzyme (organism) | Tier | Source | Label |
|---|---|---|---|---|
| acetyl/malonyl-CoA to 1,3,6,8-THN (direct) | PKS1 (*C. lagenarium*) | T1a | PMID 10501004, 23261597 | V-abs |
| to heptaketide (YWA1) | Alb1/PksP (*A. fumigatus*) | T1a/T1b | PMID 11350964 | V-abs |
| heptaketide to 1,3,6,8-THN | Ayg1 (*A. fumigatus*) | T1a (cell-free extract) | PMID 11350964 | V-abs |
| to hexaketide 2-acetyl-1,3,6,8-THN, then to T4HN | WdPKS1, WdYg1 (*W. dermatitidis*) | T1b | PMID 18676950, 11179356 | V-abs |
| 1,3,6,8-THN to scytalone | 4HNR (*M. grisea*) | T1a | PMID 10956664 | V-abs |
| 1,3,6,8-THN to scytalone | Arp2 (*A. fumigatus*) | T1b | PMID 10515939, 11350964 | V-abs |
| 1,3,6,8-THN to scytalone | THR1 (*C. lagenarium*), side activity | T1b | PMID 12967213 | V-abs |
| scytalone to 1,3,8-THN | SCD1, Arp1, SDH1 | T1b (Arp1: scytalone accumulates in the deletant; SCD1: complementation) | PMID 10515939, 8953707 | V-abs. Enzyme assays for the *Magnaporthe* enzyme are RECALLED (many papers listed on UniProt P56221, not read). |
| 1,3,8-THN to vermelone | THR1 (*C. lagenarium*) | T1b | PMID 8672814 | V-abs |
| 1,3,8-THN to vermelone | 3HNR (*M. grisea*) | T1a | PMID 10956664 | V-abs |
| 1,3,8-THN to vermelone | Arp2 (*A. fumigatus*) | T3 from my checks | - | The abstracts do not show this reaction for Arp2. |
| vermelone to 1,8-DHN | scytalone dehydratase | T3 from my checks (RECALLED T1a for *Magnaporthe*) | no Rhea reaction | Not checked. See O4. |
| 1,8-DHN to polymer | Abr2 laccase; Abr1 multicopper oxidase (*A. fumigatus*); BrnA (*A. niger*) | T1c | PMID 10515939, 16988817, 21277986 | V-abs. The substrate of each oxidase is not shown in the abstracts. |
| Botrytis, Sordaria, Podospora, Verticillium, Alternaria genes | see 3.7 | T1b or T1c | see 3.7 | V-abs or V-title |

### 5.3 L-DOPA and other melanins

| Reaction | Enzyme (organism) | Tier | Source | Label |
|---|---|---|---|---|
| L-tyrosine to L-dopaquinone; L-DOPA to L-dopaquinone | tyrosinase T (*N. crassa*); melB (*A. oryzae*); *A. bisporus* tyrosinase | T1a (RECALLED for the assay papers) | UniProt P00440; RHEA:18117, 34287; PMID 2529259, 22213164, 10552675 (titles) | V-db; V-title |
| Diphenol or catecholamine to quinone, then melanin | LAC1 (*C. neoformans*) | T1a | PMID 8300520 | V-abs |
| Tyrosine to L-DOPA by laccase | - | Not a valid edge | - | RECALLED |
| Dopaquinone to dopachrome and later steps | non-enzymatic in the design | T4 for any fungal enzyme | - | RECALLED. I did not search for fungal dopachrome tautomerase. |
| 4-hydroxyphenylpyruvate to homogentisate; homogentisate removal | HppD; HmgA (*A. fumigatus*) | T1b | PMID 19028908 | V-abs |
| 4-hydroxyphenylpyruvate to aspulvinone E; activation for polymerisation | MelA; TyrP (*A. terreus*) | T1a/T1b | PMID 27133313; RHEA:63824 | V-abs; V-db |

### 5.4 Steps where the biochemical order is uncertain

- **O1. Astaxanthin.** After echinenone the route branches: 4'-ketolation first (canthaxanthin) or 3-hydroxylation first (3-hydroxyechinenone). Both lead to phoenicoxanthin. The paper presents both branches. V full text PMID 16416328. Model it as a grid of reactions, not a line.
- **O2. Neurosporaxanthin.** The order of cyclisation, cleavage and oxidation differs between *Neurospora* and *Fusarium*. The carD paper states that the reactions "are different in their order" and that YLO-1 mediates the formation of apo-4'-lycopenoic acid. V-abs PMID 21749649. The KEGG module (lycopene to neurosporaxanthin, four KOs in a line) does not show this.
- **O3. Torulene.** In *Neurospora* the fifth desaturation comes before cyclisation; gamma-carotene is not a substrate of AL-1. V-abs PMID 11017770. Do not draw gamma-carotene to torulene as an AL-1 reaction.
- **O4. DHN late steps.** Which enzyme dehydrates vermelone in *A. fumigatus*, and which of Abr1 and Abr2 acts first, are not settled by the abstracts I read. RECALLED: published schemes differ.
- **O5. DHN first reduction.** In *C. lagenarium* two reductases act on T4HN. V-abs PMID 12967213.
- **O6. Trisporoids.** Only the two cleavage steps and two dehydrogenase steps have enzymes. The order of the later oxidations, and which steps occur in which mating type, are RECALLED as cooperative between (+) and (-) partners. TSP2 activity is found only in the (-) type although the protein is in both (V-abs PMID 18931040).
- **O7. Tyrosinase.** Rhea writes tyrosine to dopaquinone as one reaction. "tyrosine to L-DOPA to dopaquinone" as two enzyme edges is a simplification.

## 6. Recommended changes to section 6.3 and the yaml

1. DHN-melanin: replace the linear list with route variants. (a) THN synthase. (b) hexaketide PKS + Yg1. (c) heptaketide PKS + Ayg1. Then: T4HN reductase (4HNR/Arp2 type), scytalone dehydratase, T3HN reductase (3HNR/THR1/BUF1/BRM2 type), dehydration of vermelone (curated), oxidase(s) (laccase Abr2 type and multicopper oxidase Abr1 type as two families).
2. Remove "Ayg1 belongs only to the heptaketide route". Write "Ayg1/Yg1 is required with hexaketide and heptaketide PKS".
3. Remove "Arp2/THR1 type" for the first reduction. Name the two reductase types.
4. L-DOPA: make tyrosinase (EC 1.14.18.1) and laccase (EC 1.10.3.2) separate steps. Laccase carries only diphenol oxidation. Mark the *Cryptococcus* route as substrate-dependent.
5. Add melanin routes as small branches or as explicit exclusions: pyomelanin (hppD, hmgA), Asp-melanin (melA, tyrP), GHB-melanin.
6. Carotenoid yaml: split phytoene desaturase into two reactions; set the step EC to null or list both; fill the astaxanthin synthase reference with Q3HR17 (identical to CAC16261); add cao-1 as a decoy.
7. Carotenoid additions with evidence: *Phycomyces* carS and acaA; tsp3; TSP2; tdh; *Ustilago* cco1 (as a reference of the carX step); crtR as a partner in the graph.
8. Keep out, or mark T4: torularhodin step, carotenoid glycosylation, late trisporoid oxidations.
9. Key gene names by organism (carS clash).
10. Do not seed from the UniProt keyword alone (E9).

## 7. Limits of this review

- I read abstracts, not full texts (one exception). Assay details can differ from the abstract.
- "No paper found" means my PubMed queries returned none. It does not prove that none exists. No data supports a stronger statement.
- I did not check: map00906 counts, ChEBI identifiers, InChIKeys, the yaml PMIDs listed in 2.1, the AX034665 nucleotide record, the RSY1 synonym, the range of each gene across taxa.
- Statements marked RECALLED need a check before use.

## 8. References (PubMed; DOI links)

Abstract read (V-abs) unless marked.

| PMID | First author, year | DOI |
|---|---|---|
| 16416328 | Ojima 2006 (full text read) | [10.1007/s00438-005-0072-x](https://doi.org/10.1007/s00438-005-0072-x) |
| 18837978 | Alcaíno 2008 | [10.1186/1471-2180-8-169](https://doi.org/10.1186/1471-2180-8-169) |
| 18769767 | Niklitschek 2008 | [10.4067/S0716-97602008000100011](https://doi.org/10.4067/S0716-97602008000100011) |
| 17293483 | Prado-Cabrero 2007a | [10.1128/EC.00392-06](https://doi.org/10.1128/EC.00392-06) |
| 17493127 | Prado-Cabrero 2007b | [10.1111/j.1365-2958.2007.05665.x](https://doi.org/10.1111/j.1365-2958.2007.05665.x) |
| 17610084 | Saelices 2007 | [10.1007/s00438-007-0269-2](https://doi.org/10.1007/s00438-007-0269-2) |
| 18627463 | Estrada 2008 | [10.1111/j.1365-2958.2008.06349.x](https://doi.org/10.1111/j.1365-2958.2008.06349.x) |
| 21749649 | Díaz-Sánchez 2011 | [10.1111/j.1742-4658.2011.08242.x](https://doi.org/10.1111/j.1742-4658.2011.08242.x) |
| 23893079 | Díaz-Sánchez 2013 | [10.1128/EC.00084-13](https://doi.org/10.1128/EC.00084-13) |
| 19584000 | Estrada 2009 | [10.1016/j.fgb.2009.06.011](https://doi.org/10.1016/j.fgb.2009.06.011) |
| 11017770 | Hausmann 2000 | [10.1006/fgbi.2000.1212](https://doi.org/10.1006/fgbi.2000.1212) |
| 9079885 | Ruiz-Hidalgo 1997 | [10.1007/s004380050378](https://doi.org/10.1007/s004380050378) |
| 21854466 | Medina 2011 | [10.1111/j.1365-2958.2011.07805.x](https://doi.org/10.1111/j.1365-2958.2011.07805.x) |
| 22445959 | Tagua 2012 | [10.1016/j.fgb.2012.03.002](https://doi.org/10.1016/j.fgb.2012.03.002) |
| 17822929 | Burmester 2007 | [10.1016/j.fgb.2007.07.008](https://doi.org/10.1016/j.fgb.2007.07.008) |
| 18931040 | Wetzel 2009 | [10.1128/EC.00225-08](https://doi.org/10.1128/EC.00225-08) |
| 8828234 | Czempinski 1996 | [10.1099/00221287-142-9-2647](https://doi.org/10.1099/00221287-142-9-2647) |
| 10515939 | Tsai 1999 | [10.1128/JB.181.20.6469-6477.1999](https://doi.org/10.1128/JB.181.20.6469-6477.1999) |
| 11350964 | Tsai 2001 | [10.1074/jbc.M101998200](https://doi.org/10.1074/jbc.M101998200) |
| 9383199 | Tsai 1997 | [10.1046/j.1365-2958.1997.5681921.x](https://doi.org/10.1046/j.1365-2958.1997.5681921.x) |
| 16988817 | Sugareva 2006 | [10.1007/s00203-006-0144-2](https://doi.org/10.1007/s00203-006-0144-2) |
| 14970241 | Youngchim 2004 | [10.1099/jmm.0.05421-0](https://doi.org/10.1099/jmm.0.05421-0) |
| 21277986 | Jørgensen 2011 | [10.1016/j.fgb.2011.01.005](https://doi.org/10.1016/j.fgb.2011.01.005) |
| 10501004 | Fujii 1999 | [10.1271/bbb.63.1445](https://doi.org/10.1271/bbb.63.1445) |
| 23261597 | Vagstad 2012 | [10.1016/j.chembiol.2012.10.002](https://doi.org/10.1016/j.chembiol.2012.10.002) |
| 8953707 | Kubo 1996 | [10.1128/aem.62.12.4340-4344.1996](https://doi.org/10.1128/aem.62.12.4340-4344.1996) |
| 8672814 | Perpetua 1996 | [10.1094/mpmi-9-0323](https://doi.org/10.1094/mpmi-9-0323) |
| 12967213 | Tsuji 2003 | [10.1017/s0953756203008001](https://doi.org/10.1017/s0953756203008001) |
| 10956664 | Thompson 2000 | [10.1074/jbc.M006659200](https://doi.org/10.1074/jbc.M006659200) |
| 18676950 | Wheeler 2008 | [10.1128/EC.00179-08](https://doi.org/10.1128/EC.00179-08) |
| 11179356 | Feng 2001 | [10.1128/IAI.69.3.1781-1794.2001](https://doi.org/10.1128/IAI.69.3.1781-1794.2001) |
| 10541862 | Ye 1999 (not used for a claim) | [10.1007/s002940050496](https://doi.org/10.1007/s002940050496) |
| 26514268 | Schumacher 2016 | [10.1111/mmi.13262](https://doi.org/10.1111/mmi.13262) |
| 17681008 | Engh 2007 | [10.1111/j.1574-6968.2007.00867.x](https://doi.org/10.1111/j.1574-6968.2007.00867.x) |
| 27133313 | Geib 2016 | [10.1016/j.chembiol.2016.03.014](https://doi.org/10.1016/j.chembiol.2016.03.014) |
| 19028908 | Schmaler-Ripcke 2009 | [10.1128/AEM.02077-08](https://doi.org/10.1128/AEM.02077-08) |
| 8300520 | Williamson 1994 | [10.1128/jb.176.3.656-664.1994](https://doi.org/10.1128/jb.176.3.656-664.1994) |
| 11500433 | Zhu 2001 | [10.1128/IAI.69.9.5589-5596.2001](https://doi.org/10.1128/IAI.69.9.5589-5596.2001) |
| 17098915 | Frases 2007 | [10.1128/AEM.01947-06](https://doi.org/10.1128/AEM.01947-06) |
| 23042177 | Almeida-Paes 2012 | [10.1128/AEM.02414-12](https://doi.org/10.1128/AEM.02414-12) |
| 27504027 | Li 2016 (title only; abstract cut off) | [10.1016/j.simyco.2016.02.001](https://doi.org/10.1016/j.simyco.2016.02.001) |

Title only (V-title; PubMed record seen, DOI not retrieved): 8112349 (Vidal-Cros 1994), 8392512 (Kimura 1993), 11251292 (Fujii 2001), 15324811 (Watanabe 2004), 15345447 (Rodríguez-Sáiz 2004), 12839738 (Verdoes 2003), 19665150 (Schimek 2009), 29371556 (Avalos 2017), 33383912 (Marente 2020), 35430937 (Mussagy 2023), 17707627 (Coppin 2007), 33974074 (Gautier 2021), 14716498 (Kihara 2004), 29024255 (Liang 2018), 28755294 (Zhang 2017), 35637443 (Li 2022), 29485393 (Wang 2018), 17237364 (Eliahu 2007), 25226432 (Islamovic 2015), 2529259 (Kupper 1989), 16233650 (Obata 2004), 22213164 (Fujieda 2012), 10552675 (Espín 1999), 42915850 (Taratynova 2026).

Database records read on 2026-11-21: UniProt accessions named in the text; Rhea RHEA:56772, 30975, 15585, 31519, 32887, 31515, 24396, 21908, 55660, 55664, 18117, 34287, 11276, 63824, 34475, 32219, 32239; KEGG M01042, K15745, K17841, K17842, K17819, K23037, K17740, K17739, K19580, K00505, K00421, K05909 and the eight compound entries; NCBI protein CAC16261.1.
