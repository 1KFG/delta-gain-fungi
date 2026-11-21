# Unverified items: pathway `carotenoid`

- Source: `pathways/carotenoid.yaml`
- Pathway-level `verified`: **false**
- `validate()` errors: 0
- Refs: 30 total, 29 verified, 1 unverified.

## Refs with verified: false

| step | accession | species | pmid |
|---|---|---|---|
| carx | - | Ustilago maydis | 19584000 |

## Refs missing an accession or a PMID

- carx: accession=-, pmid=19584000

## Steps with no verified ref

None.

---

## Items for a person to check (hand-written, 2026-11-21, stage 1 task S3)

The section above is the output of `bin/list_unverified_refs.py pathways/carotenoid.yaml`.
The items below were written by the curator on 2026-11-21.
Labels: MEASURED = read from a database record or abstract on 2026-11-21; REVIEWED = stated by the pathway review (`docs/superpowers/specs/2026-11-21-pathway-knowledge-graph-REVIEW_PATHWAYS.md`), not re-checked; PROPOSED = curator's choice; UNVERIFIED = not checked.
The curator read PubMed abstracts only. No full text was read.

### A. References that are not verified, or verified with a gap

1. **`carx`, *Ustilago maydis* Cco1 (PMID 19584000): no accession.** The abstract shows the reaction (MEASURED). No UniProt entry could be tied to the gene name: the UniProt search `gene:cco1 AND taxonomy_id:5270` returned 0 entries and `lit_pubmed:19584000` returned 0 entries (MEASURED). UniProt has two PF03055 proteins for this species: A0A0D1E6L2 (RCO1, reviewed, resveratrol cleavage) and A0A0D1E9J4 (ORF name UMAG_00965 only, unreviewed, "Carotenoid cleavage dioxygenase 1", 787 residues, cites only the genome paper PMID 17080091). A0A0D1E9J4 is a candidate for Cco1. The gene name does not match, so it is not used and `carx.faa` has no *Ustilago* sequence. A person must confirm the locus from the full text.
2. **`cars` G8Z9M3 and `acaa` A0ABR3AZC2 (*Phycomyces blakesleeanus*).** Both UniProt entries are unreviewed and cite no PubMed paper (MEASURED). `verified: true` rests on organism and gene name (carS, acaA) and on the abstract of PMID 21854466. Identity of these sequences with the proteins assayed in the paper is UNVERIFIED.
3. **`astaxanthin_synthase` Q3HR17 and PMID 16416328.** UniProt Q3HR17 cites PMID 25504221 only (MEASURED). Q3HR17 is identical to NCBI CAC16261.1 over 557 of 557 residues (MEASURED, sequence comparison). CAC16261.1 is coded by AX034665.1:33..1706, patent EP 1035206-A1, authors Hoshino, Ojima and Setoguchi (MEASURED, GenPept record). That AX034665 is the cDNA of PMID 16416328 is REVIEWED (full text); the abstract gives no accession. B5LEC8 differs from Q3HR17 at residue 60 (G in Q3HR17, D in B5LEC8; MEASURED).
4. **`tsp3` A5HXP5 (*Blakeslea trispora*).** The abstract of PMID 17822929 does not say which species' tsp3 was expressed in *E. coli* (*B. trispora* or *Rhizopus oryzae*) and names no product. The link to the CarS reaction is T3.
5. **`phytoene_desaturase` Q8X0Z0 (*Fusarium fujikuroi* carB).** The abstract of PMID 12172798 gives 541 residues; UniProt Q8X0Z0 has 570 residues at sequence version 2 (MEASURED). The difference was not investigated.
6. **UniProt EC of P54982.** UniProt gives EC 1.3.99.30 (five-step) for *Phycomyces* carB. The abstract of PMID 9079885 reports four dehydrogenations to lycopene (MEASURED). The step keeps `ec: null`.
7. **Organism names.** `species` is the UniProt organism name. UniProt uses "Phaffia rhodozyma" for *Xanthophyllomyces dendrorhous*, "Gibberella fujikuroi (strain CBS 195.34 / IMI 58289 / NRRL A-6831)" for *Fusarium fujikuroi* IMI 58289 and "Mycosarcoma maydis" for *Ustilago maydis*. The Cco1 reference has no accession and keeps the name used in the paper.

### B. Tiers held low because the abstract does not show the evidence

1. `ggpp_synthase`, *Neurospora* al-3 (P24322): T3. The abstract of PMID 7929398 reports GGPP synthase activity kept in a deletion mutant, but not how it was measured. The abstract of PMID 1826006 reports sequence and expression only. The tier may be T1b after a full-text check.
2. `lycopene_cyclase`, *Neurospora* al-2 (P37295), gamma-carotene to beta-carotene: T3. The abstract of PMID 16928467 states that single cyclisations are preferred.
3. `lycopene_cyclase`, *Mucor* carRP (Q9UUQ6), gamma-carotene to beta-carotene: T3. The abstract of PMID 10951210 states lycopene cyclase activity and names no product.
4. `lycopene_cyclase`, *Fusarium* carRA (Q8X0Z1), all three reactions: T3. The disruptants make no phytoene, so the abstract does not show the cyclase function separately.
5. `lycopene_cyclase`, 3,4-didehydrolycopene to torulene: T1b for *Neurospora* al-2 (cyclase-domain mutants, PMID 11862485; the accumulated carotenoid was not identified), T3 for *Fusarium* carRA and *X. dendrorhous* crtYB.
6. `phytoene_synthase` and `lycopene_cyclase`, *Phycomyces* carRA (Q9P854): T1c. The abstract of PMID 11172012 reports mutant sequences and colour phenotypes.

### C. Astaxanthin route: what the abstract supports and what it does not

1. The abstract of PMID 16416328 states: one P450 does all steps from beta-carotene to astaxanthin, shown by genetics in the fungus, and "The reaction sequence is first 4-ketolation of beta-carotene followed by 3-hydroxylation" (MEASURED). It does not mention *E. coli*, canthaxanthin, 3-hydroxyechinenone, phoenicoxanthin or torulene.
2. T1a for the two ketolation reactions (RHEA:55660, RHEA:55664) comes from PMID 25504221 (crtS expressed in *Mucor circinelloides*; echinenone and canthaxanthin found; MEASURED, abstract). The review's source for T1a (canthaxanthin in *E. coli*, full text of PMID 16416328) is REVIEWED only.
3. The hydroxylation-first branch (echinenone to 3-hydroxyechinenone to phoenicoxanthin) and the torulene side route (torulene to 4-ketotorulene to 3-hydroxy-4-ketotorulene) are REVIEWED only. They are in the tables at T3 with no quote.
4. ChEBI has two isomers: CHEBI:80215 (3-hydroxyechinenone) and CHEBI:80214 (3'-hydroxyechinenone). CHEBI:80215 is used, following the name in the review. Which isomer the paper reports is UNVERIFIED.
5. The 3-hydroxylation reactions of the main line (canthaxanthin to phoenicoxanthin to astaxanthin) are T1b. The intermediate phoenicoxanthin is not named in the abstract.

### D. Compounds and reactions with gaps

1. No ChEBI entry found (ChEBI text search, 2026-11-21, MEASURED as "no hit"): 4-ketotorulene, 3-hydroxy-4-ketotorulene, trisporin, 4-dihydrotrisporin, methyltrisporate, 4-dihydromethyltrisporate. ChEBI has only the generic "trisporic acid" (CHEBI:91017, formula with R groups), which cannot have a mass. The reactions of `tsp2`, `tdh` and the torulene side route therefore have no or incomplete participants.
2. The second fragment of each cleavage by CarS and AcaA is not named in the abstract and is not in the tables.
3. Small molecules (O2, H2O, acceptor, NAD(P), diphosphate) are not listed as reaction participants.
4. `kegg_id` values are the KEGG COMPOUND cross-references of the ChEBI records (MEASURED in ChEBI; not checked in KEGG, except C08607). ChEBI CHEBI:211238 (torularhodin) lists KEGG C08607, but KEGG C08607 is neurosporaxanthin (MEASURED, one KEGG record read). The torularhodin `kegg_id` is empty. CHEBI:53157, CHEBI:53160 and CHEBI:53175 have no KEGG cross-reference in ChEBI.
5. `ko`, `kegg_rn` and `np_class` are empty: not looked up in this task.
6. The `family` column of `kb/carotenoid/step.tsv` holds the curator's descriptive labels, not database names (PROPOSED).
7. The curated cyclisation of 3,4-didehydrolycopene to torulene has no Rhea reaction (Rhea search "torulene", MEASURED). The order (fifth desaturation before cyclisation) is from the abstract of PMID 11017770.

### E. Hypothesis-level (T4) and excluded steps

1. **Torulene to torularhodin: T4.** No enzyme is identified (REVIEWED: the reviewer's PubMed queries found none; the curator ran no literature search because the fetchers have no search function). Rhea has no torularhodin reaction (MEASURED). The step is not in the YAML. The tables have the step `torularhodin_formation` with one T4 link, so torularhodin (CHEBI:211238) is in the compound table and `path_possible` never reports it.
2. **Carotenoid glycosylation or esterification: excluded.** No native fungal enzyme is known to the reviewer (REVIEWED, recalled). Not in the YAML and not in the tables.
3. **Later trisporoid oxidations (from the C18 ketone to trisporic acid, other than `tsp2` and `tdh`): excluded, T4.** No gene is known to the reviewer (REVIEWED, recalled). Not in the YAML and not in the tables.
4. *Fusarium fujikuroi* carS (a regulator, REVIEWED from titles) and carO (an opsin) are not steps.

### F. Decoys

1. All decoy lists are PROPOSED. Each decoy accession was fetched and its organism, gene name and Pfam family were checked (MEASURED). The function comes from the reviewed UniProt record. Only one decoy paper was read (CAO-1, PMID 23893079).
2. `lycopene_cyclase`: no decoy. A UniProt search for reviewed entries with NCBIfam TIGR03462 returned 23 entries (MEASURED). All are carotenoid cyclases (fungal fused proteins, *Halobacterium* CrtY, *Dietzia* C50 cyclase subunits). No characterised member with another function was found.
3. `tsp2` has one decoy (*Myxococcus xanthus* CsgA, named in the abstract of PMID 18931040 as the closest relative). The short-chain dehydrogenase family is very large and no other decoy was chosen. `tdh` has three and `crtr` has four.
4. The six decoys of the carotenoid cleavage oxygenase family are shared by `carx`, `cart`, `cars`, `acaa` and `tsp3`. The references of these five steps must also be used as outgroups of each other.
5. `crtr` is a general cytochrome P450 reductase. Other fungal P450 reductases (for example *S. cerevisiae* NCP1, P16603) do the same reaction and are not decoys. A census hit to this step does not indicate astaxanthin synthesis.
6. G2QIL8 has an ordered locus name only (MYCTH_2309496); it is used as the gene name.

### G. Schema and file layout

1. **Region coordinates (work-around for a limit of `defs.validate`).** `validate` checks a reference `region` against the record in `<step>.faa`. The records of the five fused proteins are cut to the domain, so full-protein coordinates of the C-terminal domain fail that check. The YAML `region` is therefore in record coordinates and covers the whole record. The coordinates on the full protein are in a YAML comment and in the `region` column of `provenance.tsv`. Full sequences are in `full_length/fused.faa`.
2. Cut points are the UniProt feature "Region" ("Lycopene beta-cyclase", "Phytoene synthase") of each entry (MEASURED). For all five entries the Pfam PF00494 match lies inside the phytoene synthase region and both NCBIfam TIGR03462 matches lie inside the cyclase region (MEASURED, InterPro). Pfam PF18916 matches only P37295 (residues 11-94) among the five.
3. Steps with one reference sequence: `carx`, `crtr`, `cars`, `acaa`, `tsp3`, `tsp2`, `tdh`. A profile HMM from one sequence is weak. Whether this is enough for the census is UNVERIFIED.
4. Reference coverage by lineage: Ascomycota (*Neurospora*, *Fusarium*), Mucoromycota (*Phycomyces*, *Mucor*, *Blakeslea*) and Basidiomycota (*Xanthophyllomyces*). No Pucciniomycotina, chytrid or Dothideomycete reference. *Blakeslea trispora* carRA (Q67GH9, reviewed) was fetched but not added.
5. `copy_expect`, `optional` and `role` are PROPOSED. No census data exists.
6. `crtr` is linked to RHEA:55660 (the first CrtS reaction; the crtR deletion accumulates beta-carotene). This is a modelling choice (PROPOSED).
