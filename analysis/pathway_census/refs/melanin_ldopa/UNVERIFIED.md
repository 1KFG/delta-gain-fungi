# Unverified items: pathway `melanin_ldopa`

- Source: `pathways/melanin_ldopa.yaml`
- Pathway-level `verified`: **false**
- `validate()` errors: 0
- Refs: 15 total, 14 verified, 1 unverified.

## Refs with verified: false

| step | accession | species | pmid |
|---|---|---|---|
| tyrosinase | Q2TL94 | Trametes sanguinea | 16151802 |

**Q2TL94, curator's note (2026-11-21).** The UniProt entry (unreviewed, 618 residues, EC 1.14.18.1, RHEA:18117 and RHEA:34287) has no gene name, so the gene-name rule cannot be met. `gene_name: tyrosinase` in the YAML is the protein name of the entry; it is not a UniProt gene name. The entry cites PMID 16151802, and the abstract of that paper reports monophenolase and diphenolase activity of the recombinant enzyme made in Aspergillus niger (MEASURED). The paper says Pycnoporus sanguineus; UniProt says Trametes sanguinea (taxid 158606). NCBI Taxonomy returned no taxid for "Pycnoporus sanguineus" (esearch, All Names, 2026-11-21), so there is no alias. USER DECISION at the review gate: keep the reference and its two T1a links, or remove them. The rule for entries without a UniProt gene name is the same open rule as for Q9HFV6 in `melanin_dhn`.

## Refs missing an accession or a PMID

None.

## Steps with no verified ref

None.

---

## Items for a person to check (hand-written, 2026-11-21, stage 1 task S5)

The section above is the output of `bin/list_unverified_refs.py pathways/melanin_ldopa.yaml`, plus the curator's note on Q2TL94.
The items below were written by the curator on 2026-11-21.
Labels: MEASURED = read from a database record or abstract, or computed, on 2026-11-21; REVIEWED = stated by the pathway review (`docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-REVIEW_PATHWAYS.md`), not re-checked; PROPOSED = curator's choice; UNVERIFIED = not checked.
The curator read PubMed abstracts and database records only. No full text was read.

### A. References with a gap (the user decides)

1. **`tyrosinase` Q2TL94: `verified: false`.** See the note above. No gene name in UniProt.
2. **`tyrosinase` Q2UP46 (melB of Aspergillus oryzae): ORF name only, and the name melB is not in the entry.** The entry is unreviewed, has the ORF name `AO090001000117` and 638 residues. The tie to melB is a chain (MEASURED): (i) the entry cites PMID 23749993 for the crystal structure of residues 2 to 580; (ii) the abstract of PMID 23749993 is about "Aspergillus oryzae full-length pro-tyrosinase" and cites, by journal, volume and pages, the two papers on the activation of this pro-tyrosinase: ChemBioChem 13, 193-201 and J. Biol. Inorg. Chem. 18, 19-26; (iii) PubMed gives these citations to PMID 22213164 and PMID 23053534, and both abstracts name melB. Gap: PMID 16233650 gives 616 residues for melB; Q2UP46 has 638. A0A1S9DK56 (Aspergillus oryzae, another strain, 616 residues, cites the structure paper PMID 32356371) is 91% identical to Q2UP46 in a global alignment (MEASURED). The difference was not located. USER DECISION: accept Q2UP46 as melB, or use A0A1S9DK56, or both. The T1a link (RHEA:18117) depends on this tie.
3. **`tyrosinase` P00440 (Neurospora crassa T).** The entry cites no assay paper. The T1c link rests on PMID 2147580, which the entry does not cite. The tie is the organism, the word protyrosinase, and the shared first and last author with the gene paper PMID 2529259 that the entry cites (Kupper, Lerch; MEASURED from PubMed). PMID 3037257 ("Monophenol monooxygenase from Neurospora crassa") has no abstract in PubMed.
4. **`tyrosinase` C7FF05 (Agaricus bisporus PPO4).** The entry does not cite PMID 28500345. The tie is the name: the abstract says "tyrosinase 4 from Agaricus bisporus" and the title says AbPPO4; UniProt names C7FF05 PPO4 with the alternative name "Tyrosinase 4" (MEASURED). The paper the entry cites (PMID 20676921) reports that the Escherichia coli product had no PPO activity.
5. **`laccase_diphenol` J9VY90 (LAC1) and J9VQZ4 (LAC2) of Cryptococcus neoformans H99.** UniProt gives both functions "by similarity". Neither entry cites PMID 15643074, the paper of their T1c links. The tie is the gene name and the organism; for LAC2 also the position: the abstract says LAC2 is adjacent to LAC1, and the ORF names are CNAG_03464 and CNAG_03465 (MEASURED). The abstract does not name the strain. If the strain of the paper is serotype D, the proteins are those of another taxon (Cryptococcus deneoformans); this was not checked.
6. **`laccase_diphenol` Q55P57 (LAC1).** UniProt organism: Cryptococcus deneoformans strain B-3501A (serotype D). The papers (PMID 8300520, 8760791) say Cryptococcus neoformans. NCBI gives "Cryptococcus neoformans var. neoformans" the taxid of Cryptococcus deneoformans (40410), and "Cryptococcus neoformans" and "Cryptococcus neoformans var. grubii" the taxid 5207 (MEASURED). The two species keys are not merged. So `species_key("Cryptococcus neoformans")`, the name in the papers, selects the H99 links and not the Q55P57 links.
7. **`tyrosinase` A7BHQ9 (Pholiota microspora tyr1), Q00234 (Aspergillus oryzae melO), C7FF04 (Agaricus bisporus PPO3).** Verified references; all their links are T3 (section B). UniProt organism of A7BHQ9: Pholiota microspora; the paper says Pholiota nameko (one NCBI taxid, 1538424; alias added).
8. **`asp_melanin_tyrosinase` Q0CRX0 and `aspulvinone_e_synthase` Q0CRX1.** Reviewed entries of Aspergillus terreus NIH 2624. UniProt also has A0A336U966 (TyrP of another Aspergillus terreus strain); not used.

### B. Tiers: where they differ from the pathway review, and why

Rule used (lessons of Tasks 3 and 4): a link is T1a only when the abstract that was read states the assay or the product for the named protein and organism. Deletion with intermediate analysis is T1b. A colour phenotype only is T1c. Anything else is T3; a T3 link may keep a sentence as context. One quoted sentence per link.

1. **USER DECISION. Tyrosinase: the review gave T1a (RECALLED) for Neurospora crassa T, Aspergillus oryzae melB and Agaricus bisporus. The tables give T1a to four link pairs only.** Rule of this table (PROPOSED): a tyrosinase link is T1a when the quoted sentence states an activity assay on the substrate class of the reaction (monophenol, phenolase or L-tyrosine for RHEA:18117; diphenol or L-DOPA for RHEA:34287). A sentence that states only "tyrosinase activity" or "activity determination" is held at T3. No abstract that was read names L-tyrosine or L-DOPA as the assay substrate of a reference protein. Result:
   - T1a: C7FF05 PPO4 (both reactions, "mono- and diphenols", PMID 28500345); Q2TL94 (both, "monophenolase and diphenolase", PMID 16151802); Q2UP46 melB (RHEA:18117 only, "phenolase activity", PMID 23053534).
   - T1c: P00440, RHEA:18117. An overexpression strain forms melanin on L-tyrosine (PMID 2147580). This is a colour phenotype of a gain-of-function strain, not a deletion; T1c is the nearest tier (PROPOSED).
   - T3 with a sentence as context: Q2UP46 on RHEA:34287 (PMID 22213164), Q00234 melO (PMID 7893753), C7FF04 PPO3 (PMID 21598903), A7BHQ9 (PMID 17617709). T3 without a paper: P00440 on RHEA:34287.
   - A stricter reading (the substrate itself must be named) leaves no T1a tyrosinase link. A looser reading ("tyrosinase activity" is enough) raises melO, melB on RHEA:34287 and P00440 to T1a.
   - Lead: PMID 16939623 states that TYR2 of Trichoderma reesei "was active on both l-tyrosine and l-dopa". No UniProt entry cites the paper and no entry has the gene name tyr2 (two unreviewed Trichoderma reesei tyrosinases, G0RFJ8 and G0RRF5, were found by name search), so the protein could not be tied to an accession. Not used.
2. **The quoted sentence of the PPO4 links says "the proteolytically activated mushroom tyrosinase".** In this abstract that is the recombinant tyrosinase 4. The same abstract uses "mushroom tyrosinase" for the commercial isoenzyme mixture in another sentence. The curator reads the sentence as clearly about the recombinant enzyme. If the user does not accept this, the two PPO4 links go to T3.
3. **USER DECISION. LAC1 (Q55P57), dopamine to decarboxy dopachrome: T1a (review: T1a).** The quoted sentence of PMID 8300520 states the assay and the product: "Dopamine was shown by mass spectroscopy to be oxidized to decarboxy dopachrome, an intermediate of melanin synthesis." The sentence does not name the enzyme. The abstract is about one purified enzyme and its gene (CNLAC1), and UniProt Q55P57 cites the paper for function and catalytic activity (MEASURED). By the strict naming rule the link goes to T3. The curator kept T1a because no other enzyme is in the abstract. No route result with L-DOPA depends on it.
4. **LAC1 and LAC2 on L-DOPA: T1c (review: T1a for LAC1, from PMID 8300520).** The abstract of PMID 8300520 does not mention L-DOPA. The links rest on colour phenotypes: disruption of CNLAC1 gives a Mel- mutant (PMID 8760791, Q55P57); lac1 mutants have melanin defects and a lac2 mutation delays melanin formation (PMID 15643074, J9VY90 and J9VQZ4). PMID 8760791 states the oxidation of L-dopa to quinones as background, not as a result. PMID 15643075 states "Both laccases contribute to melanin production using L-dopa as a substrate"; the sentence names neither gene, so it was not used.
5. **Laccase and homogentisic acid: T3 for both LAC1 references** (review: V-abs, no tier). PMID 17098915 states that melanisation from homogentisic acid is laccase dependent; it names no laccase gene and no strain.
6. **HppD and HmgA of Aspergillus fumigatus: T1b**, as in the review (PMID 19028908; one sentence for each gene). **HmgA of Aspergillus nidulans (Q00667): T1a**, not in the review (PMID 7673153: purified fusion protein with homogentisate dioxygenase activity).
7. **MelA (Q0CRX1): T1a** (review: T1a/T1b). PMID 27133313 states the product of MelA; PMID 28791090 states heterologous expression and the product for atmelA. Q0CRX1 does not cite PMID 28791090; the tie is the gene name.
8. **USER DECISION. TyrP (Q0CRX0): T3 for RHEA:74195 and RHEA:74199 (review: T1a/T1b).** The abstract of PMID 27133313 states that aspulvinone E "is activated for polymerization by TyrP". It states no assay and no product of TyrP. UniProt Q0CRX0 cites the paper for catalytic activity and biophysicochemical properties and gives the two Rhea reactions (MEASURED). The user may raise the links to T1a on the strength of the UniProt record. No route result changes (the route has no product, section D).
9. **Dopaquinone onward: T4**, as in the review. Two T4 links of the step `dopachrome_onward` (in `kb/melanin_ldopa/step.tsv` only): L-dopaquinone to L-dopachrome (curated) and L-dopachrome to 5,6-dihydroxyindole-2-carboxylic acid (RHEA:13041, EC 5.3.3.12). UniProt has no reviewed fungal entry with EC 5.3.3.12; the unreviewed fungal entries with this EC are automatic annotations (MEASURED). One PubMed query ("fungal dopachrome tautomerase") returned no fungal enzyme. This is not a full literature search.

### C. Leads that are not in the tables

1. **Other tyrosinases with reviewed UniProt entries:** Agaricus bisporus PPO1 (Q00024) and PPO2 (O42713), Aspergillus fumigatus tyr1 (Q8J130), Podospora anserina TYR (Q92396). Their cited papers report sequences only (PMID 12743763, abstract read) or have no PMID. Not used.
2. **Commercial mushroom tyrosinase.** PMID 10552675 reports kinetics with gamma-L-glutaminyl-4-hydroxybenzene, L-tyrosine and L-DOPA; it does not name the isoform, so no link uses it. GHB-melanin of Agaricus bisporus (review, section 3.8) is not in the tables: gamma-L-glutaminyl-4-hydroxybenzene was not looked up in ChEBI.
3. **Tyrosinases of secondary metabolism.** UniProt gives EC 1.14.18.1 to ustQ of ustiloxin biosynthesis (B8NM74, Aspergillus flavus, with RHEA:18117 and RHEA:34287; A0A8E5HUS5, Ustilaginoidea virens, EC read from a search result). The abstract of PMID 24841822 does not show tyrosinase activity of UstQ. These entries are neither references nor decoys. A census hit close to them is not a melanin tyrosinase.
4. **Zymoseptoria tritici HPPD (O42764, reviewed).** PMID 9570125: expression in Escherichia coli gives a brown pigment, prevented by an HPPD inhibitor. The product is inferred from the colour; no tier fits well, so it is not a reference.
5. **Tyrosine aminotransferase** (L-tyrosine to 4-hydroxyphenylpyruvate) and **HmgX / HmgR** of Aspergillus fumigatus (PMID 22046314) are not steps.
6. **Catechol melanin and the Ustilago maydis teliospore pigment** (review, section 3.8: RECALLED or title only). Not curated.
7. **Aspergillus fumigatus conidial extract oxidises L-DOPA** ("laccase-type activity", PMID 14970241, abstract read). The abstract names no protein. It shows the substrate range of laccases; it is not a link.
8. **KEGG.** K00505 (TYR; tyrosinase), K00421 (lccA; laccase), K05909 (laccase), K00457 (HPD, hppD) and K00451 (HGD, hmgA) are cross-references in `step.tsv` (names read from KEGG `list` on 2026-11-21; no KEGG reaction content copied). KEGG has two laccase KOs; which one holds the fungal laccases was not checked.

### D. Routes that the tables cannot follow, and the exogenous-substrate flag

1. **USER DECISION. `asp_melanin`: no product in `path_possible`.** The two TyrP products of Rhea have charged ChEBI entries only: CHEBI:193114 (C17H11O7, charge -1) and CHEBI:193115 (C17H9O6, charge -1). The ChEBI records have no conjugate-acid relation, and a search with the first block of each InChIKey found the anion only (MEASURED). The compound table holds neutral parents keyed by ChEBI id, so the TyrP reactions have no product participant, and the route stops after aspulvinone E (T1a, Aspergillus terreus, first step alone). Options: (i) a curated id scheme for compounds without a neutral ChEBI entry, the same open question as for `melanin_dhn` route (b); (ii) leave the route without a product. The products are unstable according to UniProt ("highly unstable"), so they are poor targets in any case.
2. **USER DECISION. `exogenous_substrate` is text, not a schema field.** The schema (`defs.py`) has no field for it and was not changed. Each route note starts with `exogenous_substrate: true` or `exogenous_substrate: false`; only `laccase_exogenous` is true. The `lineage_note` of step `laccase_diphenol` carries the same flag. Tests check the text. If later code must read the flag, a schema field is safer.
3. **`path_possible` does not know that a substrate is exogenous.** In the first step of a route every link is applied with no requirement on its substrates. So `laccase_exogenous` returns L-dopaquinone (T1c) and dopaminechrome (T1a) as "possible". This means: possible if L-DOPA or dopamine is supplied. A capability table must carry the flag next to this route.
4. **`pyomelanin` and `tyrosine_degradation` share the step `hppd`.** `tyrosine_degradation` is not a pigment route; it is declared so that `hmga` is counted against pyomelanin. Aspergillus fumigatus completes both with its own links (T1b). The census cannot tell pigment production from tyrosine catabolism (REVIEWED).
5. **Polymers.** L-DOPA melanin, pyomelanin, the homogentisic acid pigment of Cryptococcus and Asp-melanin have no formula. The reactions that end in a polymer have no product participant.
6. **Free L-DOPA is not a product of any reaction of the table** (Rhea writes L-tyrosine to L-dopaquinone as one reaction; review, O7). A model that needs L-DOPA as an analyte has no step that makes it.

### E. Compounds

1. **Dopaminechrome (CHEBI:189002).** The abstract of PMID 8300520 says "decarboxy dopachrome". The curator mapped this name to ChEBI "dopaminechrome (keto form)", C8H7NO2, which is the formula of L-dopachrome (C9H7NO4) minus CO2 (MEASURED for the formulas; the mapping of the name is PROPOSED). ChEBI also has an enol form (CHEBI:188999), not used.
2. **Rhea species and neutral parents.** Rhea lists zwitterions for L-tyrosine, L-DOPA and L-dopaquinone (CHEBI:58315, 57504, 57924; charge 0, charged atoms in the SMILES) and anions for the acids. The table holds CHEBI:17895, 15765, 16852, 15772, 2003, 15999, 44747, 47904 and 17704. For each pair the atoms other than H and the first InChIKey block are equal (MEASURED).
3. **PubChem cross-check.** For all 11 compounds PubChem has a compound with the ChEBI InChIKey, and the exact mass agrees with the formula mass within 0.0005 (MEASURED).
4. **Isomers.** L-dopachrome and 5,6-dihydroxyindole-2-carboxylic acid are both C9H7NO4 (193.0375). 4-hydroxyphenylpyruvic acid (C9H8O4) and homogentisic acid (C8H8O4) differ.
5. **KEGG ids** are the KEGG cross-references of the ChEBI records: C00082, C00355, C00822, C01693, C04185, C03758, C01179, C00544, C01036, C02006 (MEASURED). Dopaminechrome has none in ChEBI.
6. Dopaquinone and dopachrome may be unstable analytes (spec section 6.7, UNVERIFIED).
7. ChEBI and PubChem licences are UNVERIFIED (spec section 3). `np_class` is empty for all compounds.

### F. Decoys

1. **What separates tyrosinase and laccase.** They are not related in sequence (type-3 copper family, Pfam PF00264; multicopper oxidase family, Pfam PF07732, PF00394, PF07731), so neither is a sequence decoy of the other. They are near neighbours in function: both oxidise L-DOPA to the quinone. The annotation that separates them: tyrosinase is EC 1.14.18.1 with RHEA:18117 and RHEA:34287; laccase is EC 1.10.3.2 with the generic RHEA:11276. UniProt states for LAC1 (Q55P57) that it does not oxidise monophenols such as tyrosine (from PMID 8300520). The laccase reactions of the table are curated reactions with EC 1.10.3.2; no laccase is linked to RHEA:18117 or RHEA:34287.
2. **Rule for a decoy (checked by `build.py` for every decoy):** the entry carries the Pfam ids of the step; it carries neither the EC of the step nor a Rhea reaction of the step; its own EC and Rhea ids are those in the reason text. This is a check of the annotation. An assay that excludes the reference substrates exists for three decoys only: ivoB (Q5BGU9: "not tyrosine, 3,4-dihydroxyphenylalanine", PMID 2126551), AoCO4 (Q2UNF9: no activity on L-tyrosine and on L-dopa, PMID 19798497) and GriF (B1VTI5: no monophenolase activity, PMID 16282322). For orsC, phomQ1, gedJ, phlC, the ascorbate oxidase and abr1, activity with L-DOPA was not excluded by an assay that was read.
3. **Entries that are not decoys because they carry the EC of the step:** ustQ (B8NM74, EC 1.14.18.1) for `tyrosinase`; the laccases yA of Aspergillus nidulans (P17489), lacc of Neurospora crassa (P06811) and LAC1 of Melanocarpus albomyces (Q70KY3), all EC 1.10.3.2, for `laccase_diphenol`. Consequence: the laccase step is not specific to melanin. Any laccase may oxidise a supplied diphenol.
4. **Q2UNF9 (AoCO4) is an unreviewed entry** with the ORF name AO090001000383, no EC and the protein name "DNA, SC001". It is the Aspergillus oryzae catechol oxidase: the entry cites three crystal-structure papers of that enzyme (PMID 24043469 read; the abstract names AoCO4), and it has 408 residues, the length that PMID 19798497 gives for AoCO4 (MEASURED). PMID 19798497 is not cited by the entry.
5. **The tyrosinase decoys are far from the references.** MEASURED (global alignment, BLOSUM62, identical columns / alignment columns): the seven references are 20 to 55% identical to each other; the five decoys and TyrP are 12 to 18% identical to the references and 20 to 45% to each other (TyrP and ivoB 45%). No reviewed non-tyrosinase entry inside the long-tyrosinase group was found. So the decoys test the family border; they do not test the border between a melanin tyrosinase and another long tyrosinase.
6. **TyrP (Q0CRX0) is not a decoy of `tyrosinase`.** It is a reference of `asp_melanin_tyrosinase`. Whether it oxidises L-tyrosine or L-DOPA is not stated in the records that were read.
7. **The laccase decoys** are 22 to 29% identical to the three references (MEASURED). abr1 (Q4WZB4) is a reference of `mco_abr1` in pathway `melanin_dhn`. The bilirubin oxidase Q12737 was not used: UniProt gives it no gene name.
8. `hppd`: no decoy. The brief names no decoy family. No member of TIGR01263 with another function was looked for.
9. `hmga`: no decoy. The brief names no decoy family. No member of TIGR01015 with another function was looked for.
10. `aspulvinone_e_synthase`: no decoy. apvA (Q0CWD0, 65% identical to melA) does the same reaction, so it cannot be a decoy. Other NRPS-like enzymes of Aspergillus terreus with other products exist (PMID 29305695, abstract read: atromentin from ATEG_03090); no UniProt entry was tied to them here.
11. Identity threshold: none proposed.
