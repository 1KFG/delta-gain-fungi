# Review: pathway knowledge graph amendment (carotenoid and melanin)

Reviewed 2026-11-21. Reviewer: Claude (Fable), independent read.
Design reviewed: `2026-11-21-pathway-knowledge-graph-design.md` (not edited).
Also read: the approved census spec and my review of it (2026-11-20), the pilot lessons, the carotenoid real-data plan, and three documents from the EBC agent (`quant_v5_handoff.md`, `REVALIDATION.md`, `FOLLOWUP_integration_plan.md`). I treated the EBC documents as reported data.

Labels in this review:
- **VERIFIED**: I ran a query or read a source today. The source is given.
- **REPORTED**: a number from the EBC documents. I did not re-check it.
- **RECALLED**: from my memory of the literature. Not checked today. Treat as a lead.
- **OPINION**: my judgement.

Literature checks used PubMed. Each cited abstract has its DOI link.

## Verdict

The compound layer has real value. The three comparison models, as written, do not. Summary:

1. The compound, reaction and evidence tables are worth building. A graph store and a query layer are not needed for about three pathways.
2. No detectability data exist yet. The amendment plans three models before it knows if one target compound has one library hit in the corpus. A detectability census must come first and must be a go/no-go gate.
3. M1 (mass shifts) cannot test the pathway. Drop it.
4. M2 has about 13 fungal taxa as its true sample size. Its permutation null is degenerate as written. Its negative controls differ from the capable taxa in growth form and extraction. A negative result from M2 cannot be interpreted.
5. M3 repeats the label-based design that the EBC revalidation retracted, and the spectra it needs are not allowed to leave the EBC host.
6. The scope grew (melanin, a sterol branch, a graph, three models) while the pilot still has no real result (pilot lessons: "No real results exist").

OPINION: approve stages 1 and 2 as tables. Replace stage 4 with a staged gate (finding 1). Do not approve M1.

---

## Findings, ranked

### 1. Detectability is not measured, and it decides everything (highest importance)

**What the documents say.**
- REPORTED: 0.64% of features have a top-1 library hit. 16.2% have MS2. Negative mode has a lower MS2 rate (11.8%) than positive mode (17.3%).
- REPORTED: fungal groups are few. 247 groups touch fungi. The name-level tests used 82 fungal accessions (57 discovery, 25 held out). *S. cerevisiae* and *Komagataella* dominate.
- REPORTED (follow-up plan, section 4): "Ergosterol ionizes badly in ESI". The sterol branch proposed as a control has the same detection problem as the carotenoids.
- The amendment states the detectability request only under melanin (6.7). Section 8 lists "Carotenoids may not be detectable in the corpus at all" as a risk, with no gate.

**VERIFIED today (PubChem PUG REST, XLogP3 values).** Carotenoid targets are very lipophilic: beta-carotene 13.5, lycopene 15.6, torulene 14.6, neurosporaxanthin 11.4, astaxanthin 10.3. For comparison ergosterol is 7.4. Melanin intermediates are small and polar: scytalone 0.9, vermelone 1.2, flaviolin 0.7, L-DOPA -2.7.

**VERIFIED (PubMed abstract).** According to PubMed, carotenoids and xanthophylls form the radical molecular ion M+. and the protonated molecule [M+H]+ in ESI, in a ratio that depends on the structure ([DOI](https://doi.org/10.1002/jms.874)). The amendment's adduct list (6.5) has [M+H]+, [M+Na]+, [M+NH4]+ and [M+2H]2+. It has no M+. entry. A feature-finding tool that assumes [M+H]+ will assign a wrong neutral mass (off by 1.0073) and a wrong formula to a radical cation. The BUDDY formula layer is then wrong for the main carotenoid ion.

**RECALLED.** Standard carotenoid methods use APCI, C30 columns and non-polar solvents. A typical methanol/water extract with a C18 gradient that ends at high acetonitrile will recover little carotene. I did not check what methods the fungal groups used.

**VERIFIED (PubMed abstracts), melanin.** DHN intermediates accumulate when the pathway is blocked. According to PubMed, tricyclazole made pentaketide metabolites accumulate in *Wangiella dermatitidis* ([DOI](https://doi.org/10.1007/BF00693396)), and a *Sporothrix schenckii* mutant "accumulated and extruded" scytalone and flaviolin ([DOI](https://doi.org/10.1128/IAI.68.6.3696-3703.2000)). RECALLED: wild-type cultures polymerise the intermediates in the cell wall, and 1,3,6,8-tetrahydroxynaphthalene oxidises in air to flaviolin. OPINION: the melanin targets with the best chance in wild-type extracts are the shunt products (flaviolin, 2-hydroxyjuglone), not the on-pathway intermediates. The step list in 6.3 has no shunt compound. Section 1.3 ("melanin is the better test") is a hypothesis with no data. The amendment labels it UNVERIFIED, which is correct.

**Not checked.** Whether the GNPS libraries used by the corpus hold spectra for any target compound. My MoNA query failed (API syntax changed). This is the first thing to measure. A compound with no library spectrum is a structural zero in M2, not a negative observation.

**Recommended change: a Stage D0 gate before any model.** One request to the EBC agent, counts only:
- For each graph compound (InChIKey first block, and name): the number of library hits, by cohort (fungal reference, plant, animal, bacteria, algae, unlabeled), by polarity, at accession level, mirrors collapsed.
- Whether the library holds a spectrum for the compound at all.
- Positive-control cohorts: plant and algal groups must contain beta-carotene, lutein and zeaxanthin if the methods can see carotenoids. If plant groups show no carotenoid hits, the corpus cannot test carotenoid capability. Stop there for carotenoids.
- Write the pass rule before the counts arrive. Example (PROPOSED by me, no data behind it): at least 3 target compounds, each with hits in at least 5 independent fungal accessions.

### 2. M1 (reaction mass shifts) cannot test the pathway. Drop it.

**VERIFIED (formulas from PubChem and Rhea today; mass arithmetic is mine).** The reaction shifts of these pathways are the most common shifts in any LC-MS data set:

| Pathway step | Formula change | Mass shift (Da) |
|---|---|---|
| Desaturation (phytoene to lycopene, gamma-carotene to torulene) | -H2 | -2.0157 |
| Cyclisation (lycopene to beta-carotene) | none | 0 |
| Hydroxylation | +O | +15.9949 |
| 4-ketolation | +O -H2 | +13.9793 |
| THN reduction (T4HN to scytalone) | +H2 | +2.0157 |
| Dehydration (scytalone to 1,3,8-THN) | -H2O | -18.0106 |

OPINION, with reasons:
- Pair differences of H2, O and H2O occur between lipids, between in-source fragments, and between oxidation artefacts in every group, fungal or not. The proposed null is "random shifts of matched size". Real chemical shifts will beat random shifts in every group. The test will reject its null everywhere and say nothing about carotenoid or melanin biology.
- The cyclisation step has a zero shift. M1 cannot see it.
- M1 uses mass alone. This conflicts with the amendment's own anti-goal ("Mass alone read as identity") and with the collaboration rule in 6.8.
- REPORTED: the EBC bin-level screen (m/z and RT) found nothing above chance (FDR 0.44 to 1.0) and nothing replicated. REPORTED: 72% of features are isolated in the case-study network.
- A rescue is possible only if one member of each pair is a level-2 identified pathway compound. That is then M2 with a weaker second observation. It is not a separate model.

The amendment itself calls M1 "expected to be weak". The dependency note ("M1 needs only the graph and the feature tables") makes it the easiest model to run first. That is a risk: the easiest model is the one most likely to give a confident wrong result.

### 3. M2 (capability against observation): the sample size is taxa, and the null is degenerate

**Unit of independence.** Capability is a property of a taxon. REPORTED: the fungal taxa in the corpus are *S. cerevisiae*, *Komagataella*, three *Aspergillus*, three *Penicillium*, *S. pombe*, *Yarrowia*, *Fusarium*, *Ophiocordyceps sinensis*, *Basidiobolus*, *Agaricus*, *Pleurotus*. That is about 13 to 15 taxa. Accessions that share a taxon share its capability value. An AUC interval from accession-level resampling is too narrow. The six Eurotiales taxa are also not independent of each other.

**VERIFIED (my arithmetic, `math.comb`).** With 13 taxa and k capable taxa, the number of distinct labelings limits the smallest attainable permutation p-value: k = 2 gives 78 labelings (min p 0.013); k = 3 gives 286 (min p 0.0035); k = 4 gives 715 (min p 0.0014). No multiple-testing correction over compounds can pass with k = 2.

**The stated null does not work.** "Permute capability across taxa within study, instrument and polarity." REPORTED: the fungal groups are mostly pure cultures. A culture study usually has one taxon. A permutation inside one study then changes nothing. The null must permute taxon labels across studies inside strata, and the strata will often hold one or two taxa. State the number of taxa per stratum before the test. If most strata have one taxon, there is no valid stratified null.

**Power (illustration only; my assumptions, not corpus data).** One-sided Fisher exact test at accession level, 82 fungal accessions, zero hits in non-capable accessions:

| Capable accessions | Hits needed for p < 0.01 |
|---|---|
| 5 | 2 of 5 (p = 0.003) |
| 10 | 3 of 10 (p = 0.0014) |
| 20 | 4 or more of 20 (3 of 20 gives p = 0.013) |

These p-values are optimistic because accessions are not independent (see above). REPORTED: the genus-restricted natural products that the EBC work did detect appeared in only 3 to 8 accessions each, and "0 / 17" discovery hits replicated on held-out accessions in contrast A. Carotenoids are harder targets than those compounds. The number of capable accessions is not known; the amendment must count it before M2 is designed.

**Negative controls.**
- *Saccharomyces* and *Komagataella* are yeasts in liquid defined or rich media. The capable taxa (*Fusarium*, possibly others) are filamentous. REPORTED (follow-up plan section 4): extraction solvent drives the lipid contrast, and "the current signal is a yeast signal". Capability is then confounded with growth form, medium and extraction class. A positive M2 result can be a filamentous-versus-yeast result.
- VERIFIED (PubMed abstract): engineered *S. cerevisiae* strains produce carotenoids. According to PubMed, one strain made 8.10 mg/g DCW of astaxanthin ([DOI](https://doi.org/10.1016/j.enzmictec.2017.02.006)). RECALLED: *Komagataella* and *Yarrowia* are also common carotenoid production hosts. A public corpus can contain such strains. The reference genome says "cannot make", the sample can. Screen study titles and metadata for engineered strains before a group is used as a negative control. I did not check if such studies are in the corpus.
- The genome is a species reference. The corpus sample is a strain that the label does not name. State this limit.
- Better controls (OPINION): (a) non-capable taxa of the same growth form and extraction class as the capable taxa; (b) decoy compounds that no fungus makes but that have similar chemistry and library coverage, for example lutein (VERIFIED: same formula C40H56O2 and mass 568.4280 as zeaxanthin, different InChIKey first block). A lutein hit in a "fungal" group measures plant-medium carry-over. RECALLED: potato, carrot, tomato and maize media contain carotenoids. Beta-carotene in a culture on such a medium is not evidence of fungal synthesis.

**The response variable mixes two evidence levels.** "Observed at identification level 2 (library name) or by formula plus MS2 class." REPORTED: BUDDY and library formulas agree in 37.3% of cases (top 1) in the case study. Formula plus class is level 3 at best. Keep the two responses as separate analyses. Headline claims need level 2 (the follow-up plan says the same).

**Expression.** VERIFIED (PubMed abstract, review): according to PubMed, in the investigated *Fusarium* species "the synthesis of carotenoids is induced by light" ([DOI](https://doi.org/10.3390/jof3030039)). RECALLED: DHN-melanin is made in conidia, appressoria and other differentiated cells. Liquid cultures often lack these stages. So capable taxa will often not make the compound in the sampled condition. Consequence: M2 is one-sided. A positive association is informative. A null result has at least four explanations (not extracted, not ionised, not in the library, not expressed) plus a wrong graph. State before the run that a null result will not be read as evidence against the graph or the census.

**A positive control that the corpus can support (OPINION).** REPORTED (revalidation section 7): roquefortine C/D, andrastin A, averufin, yanuthone D and cyclopenin are detected at library level and are restricted to *Penicillium* and *Aspergillus* accessions. Averufin is a pathway intermediate (RECALLED: aflatoxin and sterigmatocystin pathway). Run the M2 method first on these: genome capability for the known clusters against observation. If M2 cannot recover these known pairs, it cannot be trusted for carotenoids or melanin. This control is stronger than the sterol branch (finding 8).

### 4. M3 (class anchors): same design that was retracted, and the spectra cannot move

- The amendment's constraint (section 3): "No spectra and no raw feature tables" cross the repository. SIRIUS is installed on HPCC, where the spectra are not. M3 can run only if the EBC host runs SIRIUS, or the constraint changes. The amendment does not say which.
- REPORTED: 267.7 M gap-filled spectra, about 80% of MGF blocks are empty placeholders. CANOPUS on the whole corpus is not realistic. The run must be limited to the fungal groups and matched controls. The cost is not measured.
- "Test enrichment of the matching classes in groups containing capable taxa against matched groups without" is a group-label contrast. REPORTED: the label-based marker panel failed on held-out studies (median sensitivity 0% on holdout; AUC 0.5 to 0.6 against plant and animal). M3 has fewer positive groups than that analysis had.
- The null ("permute class labels within instrument") permutes the wrong thing. Features inside one group and one molecular family are correlated. Permute the taxon or accession capability labels, as the revalidation did.
- Class resolution. VERIFIED (NPClassifier public API, `https://npclassifier.gnps2.org/classify?smiles=...`, scytalone SMILES): pathway "Polyketides", superclass "Naphthalenes", class "Naphthalenones". These classes cover many fungal polyketides that are not melanin intermediates. L-DOPA will fall into an amino-acid class. M3 has some meaning for carotenoids only. RECALLED: plants, algae and many bacteria make carotenoids, so a "carotenoid class" signal is not fungal.
- The same check shows the NPClassifier web API answers today. A local install with old TensorFlow (risk in section 8) is not needed for a few dozen graph compounds. RECALLED: CANOPUS in SIRIUS 5 and 6 reports NPClassifier classes directly; check this in the SIRIUS output before any mapping between class systems is built.
- The false-positive rate of CANOPUS for the carotenoid class on these spectra is unknown. Measure it on spectra with a level-2 identity first.

OPINION: keep M3 as an optional, descriptive follow-up for carotenoids, only after D0 and M2. Do not build it for melanin.

### 5. Value of a graph: keep the tables, drop the graph layer

OPINION.
- Three additions have real value over the census step table: (a) compounds with InChIKey, formula and mass, which make the corpus join possible; (b) a source and an assay type on each enzyme-to-reaction link; (c) explicit branch points and shared steps.
- All three fit in three flat tables: `compounds`, `reactions` (with substrate and product rows), and `enzyme_reaction_evidence`. The expected size is tens of compounds and tens of reactions. `can_produce` is a reachability function of about 50 lines over those tables. A graph database, a query layer and a storage decision (section 10) add work and no test.
- "A language model can then answer questions over the graph and cite each answer" (1.2) has no acceptance test. Remove it, or define the question set and the pass rule.
- G6 ("diversity of potential biology across the kingdom") is the part that the census can support. It does not need the corpus comparison. Keep G6 independent of stage 4 so that a failed comparison does not block the chapters.
- Over-scope. The pilot lessons state that no real results exist, that `copy_filtered` is 0 everywhere without a contamination table, and that every real copy is `unknown` in triage. My earlier review (finding 5) showed that the DHN steps sit in families of 2.6 M (SDR), 374 k (PKS) and about 220 k (multicopper oxidase) proteins, where orthology alone is unlikely to assign a step. The amendment adds melanin and a sterol branch to the graph before one carotenoid step call exists on real data. The graph's `present_in` edges for melanin will be no better than those step calls.

### 6. Evidence tiers and `can_produce`

**Tiers (6.2).**
- The tier is attached to the EnzymeFamily-to-Reaction edge. The uncertainty that matters is at the genome copy: how far is this genome's protein from the characterised protein? A T1 edge from *X. dendrorhous* CrtS does not make a P450 in another class T1. Add a second field on the capability row: placement support and identity to the nearest characterised reference (the census already computes both).
- "T2: shown in another organism for the same enzyme family." "Family" is not defined. The amendment's own anti-goal says carX and carT are close relatives with different substrates. Define the unit as the census step clade, not the Pfam family.
- Only T4 is excluded from "possible". T3 (orthology) is allowed. Most fungal steps outside a few model species will be T3. The capability score will then be close to "the census found the marker genes". State this. It means the graph adds little to M2 beyond the census matrix.
- T1 mixes assay strength (in vitro, heterologous, complementation, mutant phenotype) with organism. The schema records assay type; the tier ignores it. A mutant phenotype with no product identification is weaker than a heterologous assay with product identification.
- A second check by the same language model is not independent (6.4 step 6). Define "independent": a different source, or a person.

**`can_produce` logic.** The amendment leaves the rule to the builder (section 10). These points decide validity, so they belong in the design:
- *Start set.* Reachability needs declared starting compounds (GGPP, malonyl-CoA, L-tyrosine) and a rule for co-substrates. Reactions with two substrates need both. This is a hypergraph, not a simple path.
- *Multi-step enzymes.* RECALLED: one phytoene desaturase does four or five desaturations depending on the species; one cyclase makes one or two rings; one P450 does ketolation and hydroxylation (the amendment cites this for CrtS). Rhea stores these as separate reactions. If presence of the family enables every reaction, the product set is an upper bound. Lycopene, gamma-carotene, beta-carotene and torulene will all be "possible" in every genome with the two core genes. Capability then cannot separate compounds, only pathways.
- *Shared and ubiquitous steps.* GGPP synthase, multicopper oxidases and short-chain reductases are in almost every genome (my earlier review, finding 5, verified family sizes). The path call reduces to one or two marker steps: the fused synthase/cyclase and the desaturase for carotenoids; the PKS and scytalone dehydratase for DHN. Report which steps carry the call.
- *Scd1 acts twice; reductases overlap.* The path has repeated enzymes. A count of "steps present" double-counts one gene.
- *PKS product.* VERIFIED (Rhea): two routes to T4HN exist as separate reactions, RHEA:41524 (5 malonyl-CoA to T4HN, EC 2.3.1.233) and RHEA:62700 (YWA1 hydrolysis, the Ayg1 route). Presence of a non-reducing PKS in the clade does not say which product it makes (my earlier review, finding 2). The graph needs an OR between routes and an "unknown product" state.
- *Non-enzymatic steps.* Dopaquinone to dopachrome, DHN polymerisation and T4HN autoxidation have no enzyme. The schema needs a `spontaneous` reaction type. Without it, shunt products can never be "possible".
- *Regulation and compartment.* Not modelled, and correctly excluded as a non-goal. The consequence must be stated in every capability table: the table is an upper bound on production.
- *Roll-up to taxon.* "Per taxon ... the best path" is an any-genome rule. With 2,174 *Fusarium* genomes (verified in my earlier review), one contaminated or mis-annotated genome makes the genus capable. Report the fraction of genomes (and of species) with the path, and set the threshold before results.
- *Contamination.* The pilot found no per-gene taxonomy table in BFD. Bacterial or algal carotenoid genes in a lichen assembly give a false `can_produce`. The graph inherits this unsolved problem.

### 7. Identifiers, isomers and the open backbone

**VERIFIED today.**
- The KEGG trap is real. KEGG C08580 (astaxanthin) lists "PubChem: 10773" and "ChEBI: 40968" (KEGG REST). PubChem CID 10773 has formula C9H8NNaO3 (PUG REST; consistent with sodium hippurate). PubChem SID 10773 maps to CID 5281224 (PUG REST). The amendment's rule (join through ChEBI or InChIKey) is correct.
- Masses in section 4 agree with PubChem monoisotopic masses: beta-carotene 536.4382, torulene 534.4226, echinenone 550.4175, canthaxanthin 564.3967, astaxanthin 596.3866, zeaxanthin 568.4280, retinal 284.2140. Lycopene and gamma-carotene are also C40H56, 536.4382.
- InChIKey first blocks separate constitutional isomers: beta-carotene `OENHQHLEOONYIE`, lycopene `OAIJSZIZWZSQBC`, gamma-carotene `HRQKOYFGHJYEFS`; zeaxanthin `JKQXZKUSFCKOGQ`, lutein `KBPHJBAIARWVSC`. They merge stereoisomers: L-DOPA and D-DOPA are both `WTDRDQBEARUVNC`.
- Full-key matching fails between databases for one compound. Scytalone is `RTWVXIIKUFSDJB-UHFFFAOYSA-N` in ChEBI 16945 (no stereo) and `RTWVXIIKUFSDJB-SSDOTTSWSA-N` in PubChem CID 162567. First-block matching is therefore required, not optional. The amendment's proposal in 6.5 is sound.
- **Charge-state trap (new).** Rhea uses the major microspecies at pH 7.3. RHEA:31515 uses CHEBI:63069, "neurosporaxanthin(1-)": formula C35H45O2, charge -1, monoisotopic mass 497.3425, key `UGJYMKZYSUMAKJ-ZGMBEONKSA-M`. The neutral acid in PubChem is C35H46O2, 498.3498, same first block, last character `N`. If the compound table takes formula and mass from the Rhea ChEBI participant, every adduct m/z for neurosporaxanthin is wrong by 1.0073. The same applies to torularhodin and to any other acid or amine. Rule: resolve each Rhea participant to its neutral parent (ChEBI `is conjugate base of` / `is conjugate acid of`) before formula, mass and adducts are computed. The zwitterion entries (L-DOPA CHEBI:57504, L-tyrosine CHEBI:58315) have charge 0 and the correct mass; they are safe.
- **Rhea gaps (open backbone, option B).** Rhea has scytalone dehydration (RHEA:24396) and T4HN reduction (RHEA:21908). A query for "vermelone" returns no reaction. So the second reduction and the second dehydration of the DHN pathway are not in Rhea by that name. Queries for "torularhodin" and "trisporic" also return nothing. About half of the DHN pathway and all named downstream fungal carotenoid branches need hand-curated reactions. The effort for option B is larger than the table in section 5 suggests. (Limit: I searched by name only, not by ChEBI ID.)
- **The example reaction is not linked to a fungal enzyme.** RHEA:56772 (the reaction quoted in section 3) is annotated in UniProtKB to one entry, Q4VKB4, a *Chlamydomonas* ketolase (UniProt REST, `rhea:56772`). Build steps 1 and 2 (UniProt fungal entries, then Rhea by EC) will not connect the fungal astaxanthin synthase to its reactions. The link must be curated from the paper.
- **Stereochemistry of that reaction.** RHEA:56772 gives (3S,3'S)-astaxanthin. According to PubMed, astaxanthin has three optical isomers and their distribution differs by source organism, yeast included ([DOI](https://doi.org/10.1080/02652030600838399)). RECALLED: *Xanthophyllomyces* makes the (3R,3'R) isomer, the opposite of the algal product. If so, the Rhea reaction is the wrong stereoisomer for the fungal enzyme. First-block matching hides the error for the corpus join, but the graph edge would cite a reaction the fungal enzyme does not do. LC-MS on an achiral column cannot separate the isomers, so the corpus cannot say if an astaxanthin hit is fungal or algal.
- KEGG has KOs the amendment did not find: K19580 (`rppA`, 1,3,6,8-tetrahydroxynaphthalene synthase, EC 2.3.1.233) and K17739 (`THNR`, tetrahydroxynaphthalene reductase, EC 1.1.1.252) (KEGG REST `find/ko`). RECALLED: `rppA` is the bacterial type III PKS, not the fungal type I PKS; do not use K19580 as the fungal step.

**Isomer handling, further points (OPINION).**
- First-block matching merges E/Z isomers and enantiomers. That is correct for MS. Record in the compound table which ChEBI entries collapse to one first block, so that one observation is not counted for two graph compounds.
- Constitutional isomers with the same formula (beta-carotene, gamma-carotene, lycopene; zeaxanthin, lutein) are separable only at level 2 with good spectra. RECALLED: carotene MS2 spectra are similar across isomers. A library "beta-carotene" label may be any C40H56 carotene. Treat a hit as "C40H56 carotene" unless the cosine and matched peaks support more. This weakens any compound-level claim and argues for pathway-level (not compound-level) capability in M2.
- The adduct list needs M+. for carotenoids, water-loss ions for hydroxy-xanthophylls, and [M-H]- for the naphthol intermediates. [M+2H]2+ is unlikely for these compounds (RECALLED).

**Licensing.**
- VERIFIED (kegg.jp/kegg/legal.html, read today): "Academic users may freely use the KEGG website"; providers of services "are requested to obtain an academic service provider license"; non-academic use needs a commercial licence. This matches the amendment.
- VERIFIED: Rhea states CC BY 4.0 (rhea-db.org licence page). UniProt states CC BY 4.0 (UniProt help, licence page). ChEBI: the page I fetched mentions a Creative Commons licence; I did not capture the version (RECALLED: CC BY 4.0). PubChem: not checked; RECALLED: it has no single licence, and depositor records keep their own terms. The amendment labels all of these UNVERIFIED; two can now be marked checked.
- OPINION: take KEGG cross-references from ChEBI and Rhea records, which carry them under CC BY, not from the KEGG REST API. Then no KEGG-derived file enters the repository. Section 4 of the amendment already holds KEGG counts and a module's KO list; that is a factual summary in a design note, but do not push KEGG-derived tables to the shared `EBC-KT` repository. Whether that repository is public is not stated; check before pushing.
- The corpus library labels come from mixed sources (REPORTED: 21.5% are NIST14 labels). NIST content is licensed. Only counts cross the repository, which is fine; do not move spectra or names lists in bulk.

### 8. The sterol control branch

OPINION.
- Every fungus makes ergosterol (RECALLED). Capability does not vary across the fungal taxa. A constant predictor cannot be a control for M2 among fungi.
- It can serve as a detectability control only, and REPORTED evidence says it is a poor one ("ionizes badly"; the reported hit is a NIST14 "ergosterol-region" C28H40 label, a tentative label).
- The census cost is not small: the sterol branch has many steps, and squalene synthase shares the SQS_PSY family with phytoene synthase. It is useful as a decoy clade for the census (already planned in R1). It is not useful as a comparison control.
- Use the natural-product positive control in finding 3 instead.

### 9. Smaller points

- **Section 4, astaxanthin synthase.** "B5LEC8 differs from CAC16261 at one of 557 residues" is labelled MEASURED. The pilot lessons say B5LEC8 cites a genome-organisation paper, not a catalysis paper. The T1 edge for this step must cite PMID 16416328 with the quoted sentence, not the UniProt entry. I did not re-check the alignment.
- **Section 6.3, crtR.** A cytochrome P450 reductase is a general partner of all microsomal P450s. As a required step it is present in every genome and adds nothing; as a marker it is wrong. Model it as a co-factor requirement, not a pathway step (same point as my earlier review, finding 6).
- **L-DOPA.** The amendment says L-DOPA is not specific. Add: RECALLED, L-DOPA and dopachrome oxidise quickly in extracts, and tyrosine in media plus any tyrosinase gives L-DOPA. A hit is weak evidence for a melanin pathway. Keep L-DOPA melanin as a family census (as approved) and do not run a comparison model on it.
- **Mirrored duplicates.** REPORTED: 174 mirrored duplicate groups exist. Any count from the EBC side must state that mirrors were collapsed.
- **Multiple testing.** Compounds x taxa x polarity x two response definitions x three models. No rule is given. With about 13 taxa, fix one primary test per pathway before the counts arrive.
- **Section 7 order.** Stage 4 depends on stage 3 (real census), which depends on the pilot checklist items still open (decoys, contamination table, placement tool). D0 (finding 1) has no such dependency. It needs only the compound table from stage 2. Run D0 right after the compound table exists.
- **"Report AUC with intervals."** With a binary predictor (capable or not) and a rare binary response, AUC carries no more information than the two detection rates. Report the two rates with exact intervals at taxon level.

---

## Recommended order of work

1. Stage 1 and the compound, reaction and evidence tables (flat files; neutral-parent rule; first-block keys; shunt products and spontaneous steps included).
2. **D0 detectability census** through `collab/` (finding 1), with a written pass rule. This costs one request and decides the rest.
3. Census on real data (the existing real-data plan). Independent of the corpus.
4. Chapters A and B from the census and the tables (G6). Independent of the corpus.
5. Only if D0 passes: the positive-control M2 on known *Penicillium*/*Aspergillus* products, then M2 on the pathway compounds, at taxon level, pathway-level capability, one pre-registered primary test per pathway.
6. M3 only for carotenoids, only if the EBC host can run CANOPUS, descriptive.
7. M1: drop.

## What I checked and what I did not

Checked today (VERIFIED):
- KEGG C08580 fields; PubChem CID 10773 formula; SID 10773 to CID 5281224.
- PubChem formula, monoisotopic mass, InChIKey and XLogP3 for 22 compounds (carotenoids, DHN intermediates, flaviolin, L-DOPA, D-DOPA, dopachrome, ergosterol).
- ChEBI formula, charge, mass and InChIKey for 9 entries, including the neurosporaxanthin anion.
- Rhea reactions by name for astaxanthin, torulene, neurosporaxanthin, scytalone, T4HN, vermelone, torularhodin, trisporic acid, and EC 1.14.18.1; UniProt annotation of RHEA:56772.
- KEGG KO search for T4HN synthase and reductase; the KEGG legal page; Rhea and UniProt licence statements.
- The NPClassifier web API (one compound).
- Five PubMed abstracts (ESI ions of xanthophylls; two on accumulation of scytalone and flaviolin; engineered *S. cerevisiae* astaxanthin; *Fusarium* carotenoid light induction) and one on astaxanthin optical isomers by source.
- Fisher and permutation arithmetic (my own assumptions for the capable-accession counts).

Not checked (RECALLED or open):
- Whether any target compound has a spectrum in the libraries the corpus used. This is the most important open fact.
- LC and extraction methods of the fungal groups; presence of engineered carotenoid strains in the corpus.
- The (3R,3'R) configuration of *Xanthophyllomyces* astaxanthin; which *Aspergillus*, *Penicillium* and other corpus taxa carry carotenoid or DHN genes (the census must answer this).
- CANOPUS output classes in the installed SIRIUS versions; CANOPUS accuracy for carotenoids.
- ChEBI licence version; PubChem terms; whether `EBC-KT` is a public repository.
- Rhea coverage by ChEBI ID (I searched by name only).
- Any number in the EBC documents.
