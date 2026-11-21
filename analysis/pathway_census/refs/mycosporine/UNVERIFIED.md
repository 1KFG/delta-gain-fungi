# Unverified items: pathway `mycosporine`

- Source: `pathways/mycosporine.yaml`
- Pathway-level `verified`: **false**
- `validate()` errors: 0
- Refs: 8 total, 2 verified, 6 unverified.
- Fix round 1 (2026-11-21): B2J6X8 link T2 to T3; PMID 24258137 recorded; USER DECISION 3 widened; USER DECISION 7 added.

## Refs with verified: false

| step | accession | species | pmid |
|---|---|---|---|
| ddgs | WEV88679.1 | Phaffia rhodozyma (strain UCD 67-385) | 36983003 |
| omt | WEV88680.1 | Phaffia rhodozyma (strain UCD 67-385) | 36983003 |
| omt | B2J6X8 | Nostoc punctiforme (strain ATCC 29133 / PCC 73102) | 21890703 |
| atp_grasp_ligase | WEV88681.1 | Phaffia rhodozyma (strain UCD 67-385) | 36983003 |
| atp_grasp_ligase | A0A9Y1GGR7 | Phaffia rhodozyma | 36983003 |
| atp_grasp_ligase | B2J6X7 | Nostoc punctiforme (strain ATCC 29133 / PCC 73102) | 21890703 |

Reasons (curator's note, 2026-11-21):

- **WEV88679.1, WEV88680.1, WEV88681.1** (DDGS, OMT, ATPG of Phaffia rhodozyma UCD 67-385). There is no UniProtKB entry for these sequences: UniParc lists each one (477, 284 and 500 residues) with no UniProtKB accession (MEASURED). So UniProt gives no gene name, and the rule for `verified: true` cannot be met. The records are GenBank proteins of the nucleotide record OQ547787 ("Phaffia rhodozyma mycosporine biosynthetic gene cluster"). The GenBank record gives the gene names and cites the paper by title, journal and DOI 10.3390/ijms24065930, which is the DOI of PMID 36983003 (MEASURED). The full text of the paper names OQ547787 as the cluster. **USER DECISION 1:** accept a GenBank record that cites the paper as the tie. If yes, the three references become `verified: true`.
- **A0A9Y1GGR7** (ATPG, Phaffia rhodozyma strain VOH, mRNA OP021064). It is the only UniProtKB entry of a Phaffia cluster gene. It is not the assayed protein: 490 of 500 residues equal WEV88681.1 (MEASURED). The entry cites no paper. Its link is T3.
- **B2J6X8** (Npun_R5599) and **B2J6X7** (Npun_R5598) of Nostoc punctiforme. Unreviewed entries. Their only citation is PMID 23463784, "A Nostoc punctiforme sugar transporter necessary to establish a Cyanobacterium-plant symbiosis"; UniProt uses it for the genome sequence, and it is not a paper about these steps. The first version gave it as the `pmid` of the two references; fix round 1 replaced it with PMID 21890703, which the entries do not cite. The abstract of PMID 21890703 names the genes NpR5599 and NpR5598; UniProt writes Npun_R5599 and Npun_R5598. B2J6X8 carries EC 2.1.1.104 from automatic annotation and no Rhea reaction; B2J6X7 carries no EC and no reaction.
- Decoy **B2J6X6** (Npun_F5597, MysD) is `verified: false` for the same reason: the entry does not cite PMID 21890703, and the abstract writes NpF5597.

## Refs missing an accession or a PMID

None.

## Steps with no verified ref

- omt (optional=true)
- atp_grasp_ligase (optional=true)
- mg_glucosyltransferase (optional=true): no reference at all, tier T4
- glutamicol_formation (optional=true): no reference at all, tier T4

The step `ddgs` has verified references, and both are cyanobacterial (B2J6X9, Q3M6C3). No step has a verified fungal reference.

## Items for a person to check (hand-written, 2026-11-21, stage 1 task S6)

Labels: MEASURED = read from a database record, an abstract or a full text, or computed, on 2026-11-21. PROPOSED = curator's choice. REPORTED = stated by another agent or document, not checked here. UNVERIFIED = not checked.

What was read. PMID 36983003: abstract and full text (Europe PMC, PMC10057719; read by word search for the genes, the accessions, the methods and the products, not from start to end; the supplementary tables were not read). PMID 20813918, 22741921 and 28182402: abstract, and passages of the full text found by a word search (NCBI PMC); the full texts were not read from start to end. PMID 40036870: abstract only (no full text in Europe PMC). PMID 24907338: abstract only (the full-text fetch failed). Every other paper: abstract only.

### A. The fungal genes: what is known and what is not

1. **The fungal genes are known for one species.** One paper (PMID 36983003) shows three genes of one 7873 bp cluster in the yeast Phaffia rhodozyma UCD 67-385: DDGS, OMT and ATPG. Deletion of each gene stops mycosporine production; the atpg deletion accumulates 4-deoxygadusol; the three cDNAs in Saccharomyces cerevisiae give mycosporine-glutaminol (MEASURED, abstract).
2. **The fungal genes of these steps are unknown:** the glucosylation of mycosporine-glutaminol (`mg_glucosyltransferase`) and the formation of mycosporine-glutamicol (`glutamicol_formation`). Both steps have `refs: []` and tier T4. For the glucosylation the full text of PMID 36983003 states that the gene is not part of the cluster and that its function is absent from Saccharomyces cerevisiae. **An enzyme activity is published for the glucosylation, without a sequence** (PMID 24258137, Planta 1983, abstract read; missed in the first version and found by the review): a soluble UDP-glucose:mycosporin-2 glucosyltransferase from spores of Ascochyta fabae; the amide form is the best acceptor; pH optimum about 8.5; divalent cations are needed (Mn2+ more efficient than Ca2+ or Mg2+); the reaction was not reversible. The protein link stays T4, because no sequence is known and a link needs a protein. It is UNVERIFIED whether "mycosporin-2, amide form" is mycosporine-glutaminol, the aglycone of mycosporine-glutaminol-glucoside: the abstract does not say it. UDP-glucose as the donor is a lead for the enzyme family. For mycosporine-glutamicol no paper that was read gives a gene, an enzyme or a precursor.
3. **No other fungus has an assayed gene.** PMID 40036870 reports a gene cluster in Naganishia friedmannii that is induced by UV radiation (RT-PCR). The abstract names no gene and no accession, so the paper gives no reference and no link. The full text of PMID 36983003 states that the Phaffia genes were found with homologues from Pyrenophora tritici-repentis, Parastagonospora nodorum, Aspergillus clavatus and Magnaporthe oryzae (EDU45310, EAT88778, EAW13536, EDK03437, EDU45308, EAT88781); these are predicted proteins without an assay and are not references. UniProt has unreviewed fungal entries with the name "Demethyl-4-deoxygadusol synthase" (for example A0A484F947 of Colletotrichum orbiculare); they are automatic names and are not references.
4. **The step set is not a copy of the cyanobacterial gene list.** Cyanobacteria: MysA (DDGS), MysB (O-methyltransferase), MysC (ATP-grasp ligase, glycine), then MysD (D-Ala-D-Ala ligase-like) or an NRPS-like enzyme for the second amino acid (PMID 20813918, 21890703). The fungal route differs in three points (MEASURED):
   - The fungal ATP-grasp ligase adds glutaminol, not glycine. The product is mycosporine-glutaminol, not mycosporine-glycine.
   - There is no second-amine step. No paper that was read reports an iminomycosporine (shinorine, porphyra-334, palythine) from a fungus. MysD and the NRPS-like enzyme are not steps.
   - The fungal O-methyltransferase is not shown to be an ortholog of MysB: 17% identity to Npun_R5599 and 19% to Ava_3857, against 29 to 30% to catechol O-methyltransferases (global alignment). The Phaffia OMT has Pfam PF13578 at the gathering cutoff; MysB has PF01596.
   - The glucoside is the main fungal compound, and its gene is unknown.
5. **USER DECISION 2: bacterial references.** Four cyanobacterial proteins are references for tree placement (B2J6X9, Q3M6C3, B2J6X8, B2J6X7). They are never fungal evidence. Decide whether to keep them as references or to move them to a separate file.

Bacterial genes by organism (MEASURED from the abstracts and the UniProt records; "full text" where marked):

| Organism | Gene (paper name) | UniProt | Role | In the tables |
|---|---|---|---|---|
| Nostoc punctiforme ATCC 29133 | NpR5600, mysA (Npun_R5600) | B2J6X9 (reviewed) | DDGS, EC 4.2.3.154 | reference of `ddgs`, T2 |
| Nostoc punctiforme ATCC 29133 | NpR5599, mysB (Npun_R5599) | B2J6X8 | O-methyltransferase | reference of `omt`, T2 |
| Nostoc punctiforme ATCC 29133 | NpR5598, mysC (Npun_R5598) | B2J6X7 | ATP-grasp ligase, glycine | reference of `atp_grasp_ligase`, T2 on the glycine reaction |
| Nostoc punctiforme ATCC 29133 | NpF5597, mysD (Npun_F5597) | B2J6X6 | D-Ala-D-Ala ligase-like, second amino acid | decoy of `atp_grasp_ligase` |
| Anabaena variabilis ATCC 29413 | Ava_3858 | Q3M6C3 (reviewed) | DDGS, EC 4.2.3.154 | reference of `ddgs`, T3 |
| Anabaena variabilis ATCC 29413 | Ava_3857 | Q3M6C4 | O-methyltransferase (full text of PMID 20813918) | not in the tables |
| Anabaena variabilis ATCC 29413 | Ava_3856 | Q3M6C5 | ATP-grasp homolog (full text of PMID 20813918) | not in the tables |
| Anabaena variabilis ATCC 29413 | Ava_3855 | Q3M6C6 | NRPS-like enzyme (full text of PMID 20813918) | not in the tables: another fold, no fungal step |
| Actinosynnema mirum DSM 43827, Pseudonocardia sp. P1 | cluster (no gene named in the abstract) | none | shinorine in Streptomyces avermitilis (PMID 24907338) | not in the tables |

### B. Tiers (binding rule 1, applied as written)

A T1 or T2 link needs a quoted abstract sentence that names the gene of the protein. T1a also needs the assay or the product of that reaction in the sentence. A bacterial protein is T2 at most.

1. **DDGS of Phaffia, T1c, and OMT of Phaffia, T1a. USER DECISION 3 (decide the pair together).** The two tiers rest on the same sentence and the same experiment: the co-expression of DDGS and OMT in Saccharomyces cerevisiae gives 4-deoxygadusol. The full text has no strain with DDGS alone and none with OMT alone, so no control shows which enzyme does which step; the order rests on homology. The omt deletion has no reported intermediate. OMT is T1a because the rule as written is met (the sentence names OMT and the product of its reaction). DDGS is held at T1c because its own product is not named. The two consistent choices are: both T1a (the pair makes 4-deoxygadusol), or OMT at T1c like DDGS. No route result changes with either choice. Detail for DDGS: The sentence with the heterologous expression names DDGS, and the product in it is 4-deoxygadusol, which is the product of the next reaction. The product of DDGS (demethyl-4-deoxygadusol) is not named in the abstract. So the link uses the deletion sentence and is T1c (no mycosporine; no intermediate for this gene). T1a is the other reading: the co-expression of DDGS and OMT gives 4-deoxygadusol, which needs both enzymes.
2. **OMT and ATPG of Phaffia, T1a.** The sentence names each gene and the product of its reaction. Limits (MEASURED, full text): the assays are co-expressions, not single enzymes; the products were identified by HPLC retention time and absorption spectrum with co-elution, and the full text reports no mass spectrometry.
3. **B2J6X8 (mysB), held at T3 since fix round 1. USER DECISION 4.** The sentence names NpR5599 and states the product of the three genes together, not the product of the methylation. No fetched database record supports the activity of this protein: the UniProt entry is unreviewed, has no catalytic activity record with evidence, and exists "by homology". A sentence of the full text of PMID 20813918 (word search) reports 4-deoxygadusol from NpR5600 and NpR5599 together in vitro; it is not an abstract sentence. The same enzyme family as the fungal OMT is not shown (A4). T2 is the other reading. **B2J6X9 (mysA) stays T2:** the same sentence names NpR5600, and the reviewed UniProt entry gives RHEA:49560 with experimental evidence (ECO:0000269) from PMID 20813918 and 22741921. The mysC link (glycine) is T2: the sentence states the product of its reaction.
4. **Q3M6C3 (Ava_3858), held at T3.** The abstract sentence of PMID 22741921 names the enzyme class and does not name the protein. The reviewed UniProt entry cites the paper for the catalytic activity, and the full text reports the activity of Ava_3858. Raise to T2 if the UniProt record is accepted as the tie.
5. **Short gene labels.** The abstract of PMID 21890703 writes NpR5600, NpR5599, NpR5598 and NpF5597. UniProt writes Npun_R5600, Npun_R5599, Npun_R5598 and Npun_F5597. The curator reads them as the same genes (PROPOSED; same organism, same strain, same number).
6. **Ordered locus names.** UniProt gives the four cyanobacterial entries an ordered locus name and no other gene name. `kbsources.uniprot_entry` does not read ordered locus names, so it returns no gene name for them. The curation scripts read the ordered locus name from the full UniProt record. `kbsources.py` is not changed. **USER DECISION 5:** accept an ordered locus name as a gene name (as an ORF name was accepted in Task 5). If no, three entries become `verified: false`: the references B2J6X9 and Q3M6C3 and the decoy C6WFL3 (Amir_2000). Six entries of this pathway have an ordered locus name only (B2J6X9, Q3M6C3, B2J6X8, B2J6X7, B2J6X6, C6WFL3). If yes, the Task 7 verifier must read `orderedLocusNames`, or its gene-name check fails for these six.

### C. Routes

- No declared route gives a product in the tables (D1). The route prefix `ddgs` gives demethyl-4-deoxygadusol at T1c for Phaffia rhodozyma. The step `atp_grasp_ligase` alone gives mycosporine-glutaminol at T1a for Phaffia rhodozyma.
- `mycosporine_glutaminol_glucoside` and `mycosporine_glutamicol` end in a T4 step. A path with a T4 step is never reported as possible.
- `mycosporine_glutamicol` lists `ddgs` and `omt` as its first steps. That is an assumption from the shared 4-deoxygadusol core (PROPOSED). It is not shown in a glutamicol-producing fungus.
- The T4 steps are in the YAML with `refs: []`. The carotenoid and melanin tables keep T4 steps in `step.tsv` only. Here the route variants need the step ids (PROPOSED).

### D. Compounds

1. **USER DECISION 6: compounds with no ChEBI id.** A ChEBI search by name, by InChIKey and by the first InChIKey block found no ChEBI entry for these compounds (MEASURED). No id is invented, and the compounds are not in `compound.tsv`:

| Compound | Formula (PubChem) | PubChem CID | Effect |
|---|---|---|---|
| 4-deoxygadusol | C8H12O5 | 85926716 | the product of `omt` and the substrate of `atp_grasp_ligase` are missing: no route gives a product |
| mycosporine-glycine | C10H15NO6 | 14444485 | the bacterial reaction has no participant |
| mycosporine-glutaminol-glucoside | C19H32N2O11 | 171990409 | the main fungal compound is not in the tables |
| mycosporine-glutamicol-glucoside | C19H31NO12 | 171990350 | not in the tables |
| asterina-330 | not found in PubChem by this name | none | not in the tables |
| mycosporine-alanine | C11H17NO6 | 101450154 | not in the tables |
| mycosporine-glutamine | not found in PubChem by this name | none | not in the tables (reported in Colletotrichum graminicola) |
| glutaminol | not resolved (PubChem returns L-glutamic acid for the name) | none | not listed as a substrate |

2. ChEBI has gadusol(1-) (CHEBI:228248) and demethylgadusol(1-) (CHEBI:228249) as anions only under those names. Neutral gadusol is CHEBI:172445 under a systematic name (same InChIKey as PubChem "gadusol", MEASURED). Neutral demethylgadusol was not found, so KEGG C22843 has no compound row.
3. **PubChem check.** For shinorine (CHEBI:131712) and porphyra-334 (CHEBI:35671) PubChem has no compound with the ChEBI InChIKey, and the PubChem compounds of those names have another first block (the imine is drawn on the other nitrogen). Palythine (CHEBI:182637) has the same first block as PubChem "palythine". The three are in `compound.tsv` with the ChEBI structure (PROPOSED). PubChem "shinorine" has the first block of CHEBI:175245 (Mytilin A).
4. **KEGG.** C05382 and C10607 are KEGG cross-references of the ChEBI records. C21281 (demethyl-4-deoxygadusol) and C22846 (gadusol) are matched by KEGG name and formula, because neither database links the other (PROPOSED). KEGG C21281 has no stereo centre in its name; the table uses the (R) form of RHEA:49560. KEGG ids are cross-reference columns only.
5. CHEBI:7039 has the name "Mycosporine". Its InChIKey is that of PubChem "mycosporine serinol" (MEASURED). Gadusol, mycosporine-serinol and palythine are context compounds: no reaction of the table makes them.
6. Rhea: RHEA:49560 is the only reaction of the fungal route in Rhea. Searches for gadusol, deoxygadusol, mycosporine, shinorine and valiolone found no other reaction of the route. The other reactions have curated ids.

### E. Decoys

- The check is a check of the UniProt annotation: family id on the entry, no EC and no Rhea reaction of the step. For `omt` and `atp_grasp_ligase` the step has no EC and no Rhea reaction, so the check is only that the decoy carries another annotation.
- Assay read (abstract): C6WFL3 (2-epi-valiolone synthase, PMID 22741921) and B2J6X6 (MysD, PMID 21890703). All other decoys: annotation only. No decoy was assayed with the substrate of the fungal step, so no decoy is shown to lack the reference reaction.
- `omt`: the Phaffia OMT is 29 to 30% identical to the decoys O42898 and P21964 and 17% to the bacterial reference (MEASURED). The decoys carry PF01596; the Phaffia OMT has PF01596 below the gathering cutoff only.
- **USER DECISION 7: decoy O42898 is not excluded as a true member of the `omt` step.** It is "Probable catechol O-methyltransferase 1" of Schizosaccharomyces pombe. Its catalytic activity in UniProt has no evidence code (MEASURED), and no assay of it was read. It is closer to the Phaffia OMT (30%) than MysB is (17%). MysB, the bacterial enzyme of this reaction, has the Pfam family of the decoys (PF01596). So a census must not count a hit to O42898 as a negative; the `omt` step needs tree placement. Keep O42898 as a decoy, or drop it. Whether Schizosaccharomyces pombe makes a mycosporine was not checked (UNVERIFIED). P21964 (human COMT) has an annotated function on catechols and was not assayed with demethyl-4-deoxygadusol.
- `atp_grasp_ligase`: the family id of the decoys is the PROSITE profile PS50975. It is on the UniProt entry of the strain VOH protein (A0A9Y1GGR7) and of MysC (B2J6X7). It was not computed for WEV88681.1.
- `mg_glucosyltransferase`: no decoy (no reference). `glutamicol_formation`: no decoy (no reference).
- Not decoys and not references: the NRPS-like enzyme Q3M6C6 (another fold); the fish 2-epi-5-epi-valiolone synthase of gadusol biosynthesis (PMID 25965179; no accession was fetched).

### F. Lineages (`fungal_lineages.tsv`)

- 29 rows: organism, compound, PMID, and whether the abstract or the full text was read. The table is a list of reports. It is not generalised to the kingdom.
- The reports are yeasts of the Basidiomycota (Phaffia, Naganishia, Rhodotorula, Cryptococcus, Cystobasidium, Pucciniomycotina), a few ascomycetes (Hortaea werneckii, rock-inhabiting microcolonial fungi, Colletotrichum graminicola, Ascochyta fabae) and lichen thalli. For the lichens the abstracts do not say which partner makes the compound.
- Negative reports are kept: Phaffia rhodozyma CBS 6938 (no cluster, no mycosporine) and Metarhizium anisopliae conidia.
- Organism names are the names in the abstracts. They were not checked against NCBI Taxonomy, except Phaffia rhodozyma and Naganishia friedmannii. Current names may differ (UNVERIFIED).
- The list is not complete. The older literature on filamentous fungi was not searched in full.

### G. Context from the EBC agent (context only, not evidence for the pathway)

REPORTED by the EBC agent (`IDEA_environment_host_signatures.md`), not checked here: names of mycosporine-like amino acids occur in marine accessions (36/87, odds ratio 16.9) and in 0 of 82 fungal-reference accessions. The task brief gives other numbers for the same report (176 corpus accessions; 32 of 88 marine accessions, odds ratio 14.5). Both sets are REPORTED. This is context only. It is not evidence for or against any step, tier or route of this pathway. The compounds that the fungal literature reports (mycosporine-glutaminol-glucoside and its relatives) are not the iminomycosporines of marine algae and cyanobacteria (shinorine, porphyra-334, palythine). Whether the EBC names include the fungal compounds is not known here.
