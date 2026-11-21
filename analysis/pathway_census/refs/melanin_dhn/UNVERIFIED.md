# Unverified items: pathway `melanin_dhn`

- Source: `pathways/melanin_dhn.yaml`
- Pathway-level `verified`: **false**
- `validate()` errors: 0
- Refs: 22 total, 22 verified, 0 unverified.

## Refs with verified: false

None.

## Refs missing an accession or a PMID

None.

## Steps with no verified ref

None.

---

## Items for a person to check (hand-written, 2026-11-21, stage 1 task S4)

The section above is the output of `bin/list_unverified_refs.py pathways/melanin_dhn.yaml`.
The items below were written by the curator on 2026-11-21.
Labels: MEASURED = read from a database record or abstract, or computed, on 2026-11-21; REVIEWED = stated by the pathway review (`docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-REVIEW_PATHWAYS.md`), not re-checked; PROPOSED = curator's choice; UNVERIFIED = not checked.
The curator read PubMed abstracts only. No full text was read.

### A. References that are verified with a gap (the user decides)

1. **`t4hn_reductase` Q9HFV6 (4HNR): the UniProt entry has no gene name.** `verified: true` is set although the rule asks for a gene name match; no gene name exists to compare (MEASURED). `gene_name` in the YAML is `4HNR`, the protein abbreviation of PMID 11467929. The tie to the paper is by the record: the entry cites PMID 10956664 and PMID 11467929 (PDB 1JA9), and its organism is Pyricularia grisea (NCBI synonym Magnaporthe grisea, taxid 148305) (MEASURED). If the user wants the rule applied strictly, set `verified: false`.
2. **`t3hn_reductase` Q12634: ORF name only (`MGG_02252`), and the UniProt name conflicts with the papers.** UniProt names Q12634 "Tetrahydroxynaphthalene reductase", EC 1.1.1.252. The pathway review used it as 4HNR. It is the trihydroxynaphthalene reductase (3HNR): it cites PMID 8112349 (the 1994 reductase cDNA) and two crystal-structure papers of 1,3,8-trihydroxynaphthalene reductase (PMID 8939741, 11342131; PDB 1YBV), and it does not cite PMID 10956664 (MEASURED). Q9HFV6 and Q12634 are 45% identical in a global alignment (MEASURED); the abstract of PMID 10956664 gives 46% for 4HNR against 3HNR. The abstract of PMID 10956664 names no accession. The UniProt organism of Q12634 is Pyricularia oryzae strain 70-15; that of Q9HFV6 is Pyricularia grisea. NCBI holds them as two species (taxid 318829 and 148305), so the tables do not merge them, although the paper calls both enzymes Magnaporthe grisea.
3. **Aspergillus niger A2QUI2, A2R3S8, A2R3T0 (PMID 21277986).** The UniProt entries do not cite the paper. The tie is the locus tag in the abstract: An09g05730 (FwnA), An14g05350 (OlvA), An14g05370 (BrnA) (MEASURED). A2R3S8 and A2R3T0 are unreviewed and have the ORF name only; `gene_name` is the ORF name. A2QUI2 is reviewed as `albA`; UniProt ties it to aurasperone B (a naphthopyrone) through PMID 31067027 (not read). A2R3T0 carries an automatic EC 1.16.3.1 (ferroxidase).
4. **`thn_synthase` P79068.** Unreviewed. The UniProt organism name is "Glomerella lagenarium" with taxid 5462; NCBI Taxonomy names taxid 5462 Colletotrichum lagenaria and has no name "Glomerella lagenarium" or "Colletotrichum lagenarium" (esearch returned nothing; MEASURED). The other Colletotrichum references (P87025, Q00455, F8WTI4) are Colletotrichum orbiculare (taxid 5465, strain 104-T taxid 1213857). The papers call all of them Colletotrichum lagenarium, and PMID 8672814 names the wild type 104-T. The species keys are NOT merged (different taxids). Effect: the tables do not let one Colletotrichum complete route `pks_thn`. The user decides whether P79068 belongs to Colletotrichum orbiculare. The entry cites PMID 7500937 only; the product papers (PMID 10501004, 23261597) are tied by gene name, organism name and the length 2187 (MEASURED in the abstract of PMID 7500937 and in the entry).
5. **`hexaketide_pks` Q9Y7A7 and `ayg1` Q6BCB5.** Both unreviewed. UniProt organism Exophiala dermatitidis; the papers say Wangiella dermatitidis (one NCBI taxid, 5970; alias added). The UniProt gene name of Q6BCB5 is `yg1` (paper: WdYG1).
6. **`t4hn_reductase` F8WTI4 (T4HR1, Colletotrichum orbiculare).** Unreviewed. The entry cites PMID 12967213 for the sequence. The abstract deduces a T4HN-specific reductase and does not name T4HR1. Link tier T3.
7. **Sordaria macrospora A7DWM4 (teh), A7DWM6 (tih), A7DWM5 (sdh).** Unreviewed; the entries cite PMID 17681008; the abstract reports homology only (and RNA interference of sdh without a stated phenotype). All links T3. The Sordaria pks (A7DWM3) is not a reference: the abstract does not say which of the three PKS routes it uses.

### B. Tiers: where they differ from the pathway review, and why

Rule used: a link is T1a only when the abstract that was read states the assay or the product. Deletion with intermediate analysis is T1b. Deletion with colour phenotype only is T1c.

1. **4HNR and 3HNR accessions swapped against the review.** See A2. T4HN to scytalone T1a is on Q9HFV6, T3HN to vermelone T1a is on Q12634 (PMID 10956664).
2. **WdPKS1 (Q9Y7A7), hexaketide synthesis: T1c, review T1b.** The abstract of PMID 11179356 shows an albino disruptant. The hexaketide is shown in the WdYG1 deletion mutant (PMID 18676950), not as the product of the PKS.
3. **Alb1 (Q4WZA8): T1a from PMID 11040426** (product identified in Aspergillus oryzae). The review gave "T1a/T1b" from PMID 11350964.
4. **SCD1 (Q00455), scytalone dehydration: T1c, review T1b.** The abstract of PMID 8953707 reports complementation of a melanin-deficient mutant to the wild phenotype and no intermediate.
5. **SDH1 (P56221), scytalone dehydration: T1a from PMID 9466792** (review: assay papers RECALLED). The abstract reports a kinetic assay of scytalone dehydratase of "Magnaporthea grisea" with scytalone.
6. **Vermelone dehydration: T1a for P56221 from PMID 10320327; T3 for Q00455, O14434, A7DWM5** (review: T3, RECALLED T1a). The abstract names "the two physiological substrates, scytalone and veremelone" and reports relative substrate specificities. It does not name the organism or the product. The tie to P56221 is the UniProt entry, which cites this PMID for catalytic activity and for the mutagenesis of Tyr-30, Asp-31, Tyr-50, His-85, His-110, Ser-129 and Asn-131, the residues of the abstract (MEASURED). If the user requires the organism in the abstract, lower this link to T3; route `pks_thn` and `pks_heptaketide` then end at T3.
7. **THR1 (P87025), T3HN to vermelone: T1b** as in the review. The abstract of PMID 8672814 states the metabolic defect of the complemented mutant; it reports no intermediate measurement of its own. A stricter reading gives T1c.
8. **Arp2 (E9QUT3), T3HN to vermelone: T3** as in the review. Not shown in the abstracts.
9. **PfmaG (W3XAA0): T1a, not in the review's tier table.** Scytalone identified after expression of the PKS and PfmaG in Aspergillus nidulans (PMID 28517364). The reductase was not assayed alone.
10. **Aspergillus niger FwnA, OlvA, BrnA: T1c.** Colour mutants complemented and genes disrupted (PMID 21277986). The abstract names no enzyme product. Whether the Aspergillus niger conidial pigment is DHN-melanin is UNVERIFIED.
11. **Abr1 (Q4WZB4): the quoted sentence of PMID 10515939 is not the sentence that supports this link.** One sentence per paper is quoted. The quoted sentence reports the intermediates of the arp1 and arp2 deletants (needed for two T1b links). The abr1 colour phenotype is in another sentence of the same abstract ("disruption of each gene"); the `assay` field says so. The user may allow a second sentence for this paper.
12. **4HNR (Q9HFV6) on T3HN and 3HNR (Q12634) on T4HN: T1a side-activity links.** The competition experiment of PMID 10956664 shows both enzymes accept both substrates (310-fold and 4.2-fold preference). These links record the substrate overlap. They do not change any route result.

### C. Leads that are not in the tables

1. **BUF1 and BRM2.** The step name "3HNR/THR1/BUF1/BRM2 type" comes from the spec. No reference carries these gene names. Alternaria alternata Brm2 is UniProt Q76K49 (250 residues, unreviewed); its only citation, PMID 15099381, has no abstract in PubMed. PMID 8392512 (ALM, BRM1, BRM2 cluster; abstract read) shows complementation and disruption with colour phenotypes; no UniProt entry cites it. BUF1: a UniProt search `gene:BUF1` in fungi returned only a yeast replication protein (MEASURED).
2. **Botrytis cinerea** (PMID 26514268): bcpks12, bcpks13, bcygh1, bcbrn1, bcbrn2, bcscd1. UniProt has Q6RKI8 (PKS12) and Q6RKI7 (PKS13), both unreviewed, citing PMID 14676319 (not read). No entry cites PMID 26514268. Not used. The lineage restriction is in `lineage_note`.
3. **Pestalotiopsis fici PfmaE (W3X7U2, reviewed), Verticillium dahliae PKS1 (G2WRK7), Zymoseptoria tritici PKS1 (F9XMW3).** Melanin PKS with reviewed entries. The abstracts read do not show which PKS route they use, so they are not references of a PKS step. MEASURED: PfmaE is 71% identical to P79068 (THN synthase) and 44 to 45% to the references of the other two PKS steps.
4. **Pestalotiopsis fici PfmaD (W3X7K0) and Verticillium dahliae LAC1 (G2WRK6).** Reviewed "multicopper oxidase" entries of DHN-melanin clusters. MEASURED: 34 to 37% identical to Abr1 and BrnA, 22 to 23% to Abr2. Not used (abstracts not read for the oxidase function). Note the gene name: LAC1 of Verticillium dahliae is a DHN-melanin cluster oxidase; LAC1 of Cryptococcus (decoy Q55P57) is the diphenol laccase. Same name, different roles.
5. **Second reductase in Aspergillus fumigatus.** Not checked. The cluster has arp2 only.
6. **KEGG.** K17739 is a cross-reference of `t4hn_reductase`, K17740 of `scytalone_dehydratase` (names read from KEGG `list`; no KEGG reaction content copied). KEGG has one KO for the naphthol reductases; which paralogue type its members are was not checked. K19580 (rppA) is bacterial and is not used.
7. **Gene names arp1 and arp2.** In Aspergillus fumigatus they name the scytalone dehydratase and the reductase. In Saccharomyces cerevisiae ARP1 and ARP2 name actin-related proteins (UniProt P38696, P32381; MEASURED). A gene name is valid only with its organism.

### D. Routes that the tables cannot follow, and single organisms

1. **`pks_hexaketide`: no product.** 2-acetyl-1,3,6,8-tetrahydroxynaphthalene (AT4HN) has no ChEBI entry (four name searches and one InChIKey search returned nothing; MEASURED). PubChem has it: CID 86018477, C12H10O5, InChIKey FYPSLERIHMOYRN-UHFFFAOYSA-N (MEASURED). The compound table is keyed by ChEBI id, so the PKS reaction and the Yg1 reaction have no AT4HN participant, and `path_possible` stops after the PKS step. From the Yg1 step onwards the steps reach 1,8-dihydroxynaphthalene (test `test_tail_of_the_hexaketide_route_reaches_dhn`). The user decides whether a compound without a ChEBI id may enter the table.
2. **`polymerisation`: no product.** The polymer has no formula.
3. **No organism completes a PKS route with its own links.** `pks_thn` and `pks_heptaketide` reach 1,8-dihydroxynaphthalene at T1a only with links of three organisms (`single_organism` False). Aspergillus fumigatus completes 4 of 6 steps of `pks_heptaketide` (1,3,8-trihydroxynaphthalene, T1b) and has no reference in `t3hn_reductase`. Colletotrichum orbiculare and Sordaria macrospora complete the steps from T4HN to 1,8-dihydroxynaphthalene at T3 and have no PKS reference under the same species key (see A4, A7).
4. **3-acetylflaviolin** (PMID 18676950) has no ChEBI entry found and no PubChem compound under that name (MEASURED). Not in the table.

### E. Compounds

1. **Flaviolin (CHEBI:42646): the PubChem cross-check fails.** PubChem has no compound with the ChEBI InChIKey RROPNRTUMVVUED-UHFFFAOYSA-N (HTTP 404). The PubChem compound named flaviolin is CID 160478, 4,5,7-trihydroxynaphthalene-1,2-dione, InChIKey XNPCAGMCQDGQKK-UHFFFAOYSA-N: same formula C10H6O5 and mass, another first block (a tautomer; ChEBI draws 2,5,7-trihydroxy-1,4-naphthoquinone) (MEASURED). The table keeps the ChEBI key. A match on the first block of the InChIKey will miss flaviolin in a resource that uses the PubChem structure.
2. **Vermelone (CHEBI:9957)** is a 2-star ChEBI entry without stereochemistry. ChEBI also has (+)-vermelone (CHEBI:198683). UniProt Q12634 writes (-)-vermelone. The first block is the same.
3. **Scytalone.** ChEBI and PubChem (CID 3477029) agree on the full key RTWVXIIKUFSDJB-UHFFFAOYSA-N (no stereochemistry) (MEASURED). The spec's statement that the full keys differ was not reproduced for this ChEBI entry.
4. **KEGG ids** are the KEGG cross-references of the ChEBI records: C00779, C04033, C01173, C01624, C18012, C00024, C00083, C00164 (MEASURED). YWA1 and naphthalene-1,8-diol have none in ChEBI.
5. ChEBI and PubChem licences are UNVERIFIED (spec section 3).
6. `np_class` is empty for all compounds.

### F. Decoys

1. **wA (Q03149) is a decoy of `thn_synthase` and `hexaketide_pks` only.** It makes YWA1 (RHEA:62652), the reference reaction of `heptaketide_pks`, so it cannot be a decoy there. The brief named it for "the PKS clades". Consequence: a hit in the heptaketide clade (Alb1, FwnA, wA; 68 to 70% identical to each other, MEASURED) does not show DHN-melanin.
2. `ayg1`: no decoy. The brief names no decoy family. PMID 38491909 (abstract read) describes a bacterial polyketide-trimming hydrolase (PlAntI) related in function; no fungal non-melanin member of PF06500 was looked for.
3. `scytalone_dehydratase`: no decoy. The brief names no decoy family. No member of PF02982 with another function was looked for.
4. **Versicolorin reductase decoys are close to the 3HNR type.** MEASURED (global alignment, BLOSUM62): P50161, Q00791 and Q8TFD5 are 59 to 63% identical to Q12634, P87025 and A7DWM6, and 44 to 49% to the 4HNR-type references. The 4HNR-type references are 43 to 50% identical to the 3HNR-type references. A hit cannot be assigned by identity alone.
5. **Arp2 and Sordaria teh are not close to 4HNR.** MEASURED: E9QUT3 is 52% identical to Q9HFV6, 44 to 48% to the 3HNR type and 46 to 48% to the versicolorin reductases. A7DWM4 is 45 to 52% to every group. Q9HFV6, W3XAA0 and F8WTI4 are 65 to 75% identical to each other. Whether Arp2 and teh belong to the 4HNR clade needs the tree.
6. **yA (P17489) is the nearest characterised relative of Abr2** (39% identical; all other pairs with Abr2 are 19 to 25%; MEASURED). It is a decoy because its record states another pigment (the wA product). That Aspergillus nidulans makes no DHN-melanin from it is REVIEWED (RECALLED there), not checked here.
7. Decoy functions were read from the UniProt record (name, EC, Rhea). One decoy abstract was read per pmid given (Q03149: PMID 11251292; P17489: PMID 2659435). P50161: the UniProt function text starts with "Cytochrome P450 monooxygenase"; its name, EC and Rhea reaction are those of versicolorin reductase.
8. Identity threshold: none proposed. The decoy data above show that identity does not separate the reductase groups.
