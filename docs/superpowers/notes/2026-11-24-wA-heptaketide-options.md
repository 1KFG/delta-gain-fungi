---
title: "wA (Q03149) and the step heptaketide_pks: options to keep its information"
subtitle: "DHN melanin, pathway census step 5 (joint PKS tree). Options note for the owner."
date: 2026-11-24
author: "Claude agent (read-only analysis; no reference, decoy or result changed)"
---

::: {.callout-important}
## What this is, and where it stops being reliable

- An options note. It changes nothing. No reference, decoy, HMM, `min_score` or result file was edited. The choice is the owner's.
- All census results are PROVISIONAL (`pathways/melanin_dhn.yaml` has `verified: false`).
- Labels: **MEASURED** = read from a file named here, or from a record fetched on 2026-11-24. **REPORTED** = stated in a PubMed abstract or a UniProt record that was read on 2026-11-24 (id given). **PROPOSED** = a choice. **UNVERIFIED** = not checked.
- Only abstracts were read. No full text was read.
- The tree has 2 raxml-ng starts and no bootstrap. All tree statements depend on one topology (`results/step5/melanin_dhn/pks_joint/ref.raxml.bestTree`).
- The sub-region counts in section 4 use the single best epa-ng edge of each protein. For 49 of 2,200 proteins that edge has LWR under 0.5 (MEASURED). The clade-level calls of `variant_*` do not depend on this.
- The scripts for section 4 ran from a scratch folder and are not in the repository. The method is described in section 4 so it can be repeated.
- Nothing here shows the function or the product of any census protein.
:::

Paths are relative to `analysis/pathway_census/` unless they start with `docs/`.

## 0. Short answer

- **wA is not a decoy of `heptaketide_pks` in the curated files.** It is a decoy of `thn_synthase` and `hexaketide_pks` only. The decoy provenance row says so in words (MEASURED, `refs/melanin_dhn/decoys/provenance.tsv`; git history below). It became a boundary for the heptaketide target because the joint tree uses one union decoy set for all three steps (MEASURED, `bin/step5_joint.py`, `define_joint_targets`).
- wA, Alb1 and AlbA carry the same reaction (RHEA:62652, YWA1) and the same domain order on their UniProt records (MEASURED). They are 67.7 % to 69.9 % identical to each other (MEASURED). On the census data the three regions follow the species: A. nidulans 19 of 19 genomes on the wA side, A. fumigatus 392 of 393 on the Alb1 side, A. niger 122 of 122 on the AlbA leaf (MEASURED).
- Recommendation (PROPOSED): option B, a step-specific decoy rule in the joint tree, so that wA is a labelled member inside the heptaketide target. The result already exists as `variant_wA`: 2,123 of 2,200 proteins assigned, in 2,091 genomes. No HMM, `min_score` or search re-run is needed.

## 1. What wA, Alb1 and AlbA make

### 1a. UniProt records (MEASURED, fetched from rest.uniprot.org on 2026-11-24)

| Accession | Organism (taxid) | Gene (ORF) | Name | Length | Entry / seq version | Rhea | Function text (short) |
|---|---|---|---|---|---|---|---|
| Q03149 | Aspergillus nidulans FGSC A4 (227321) | wA (AN8209) | Conidial yellow pigment biosynthesis polyketide synthase | 2,157 | 163 / 2 | RHEA:62652, cited PMID 11251292 | NR-PKS that forms the heptaketide naphthopyrone YWA1, a yellow pigment of conidia |
| Q4WZA8 | Aspergillus fumigatus Af293 (330879) | alb1, syn. pksP (AFUA_2G17600) | Conidial pigment polyketide synthase alb1 | 2,146 | 132 / 1 | RHEA:62652, cited PMID 11040426 | first step of DHN-melanin: YWA1; Ayg1 then makes 1,3,6,8-THN |
| A2QUI2 | Aspergillus niger CBS 513.88 (425011) | albA (An09g05730) | Non-reducing polyketide synthase albA | 2,153 | 109 / 1 | RHEA:62652, cited PMID 31067027 | first step of aurasperone B: YWA1 |

All three records have the same reaction text: "6 malonyl-CoA + acetyl-CoA + 6 H(+) = naphtopyrone YWA1 + 6 CO2 + 7 CoA + H2O" (MEASURED).

The AlbA record ties AlbA to aurasperone B, a naphthopyrone, not to DHN melanin (MEASURED). So one of the two current references already has a UniProt function that is not DHN melanin. Its DHN-melanin link comes from PMID 21277986 (FwnA, colour phenotype, T1c in `kb/melanin_dhn/step_reaction.tsv`).

### 1b. Abstracts (REPORTED)

- PMID 11251292 (wA): "Only the 32 amino acid deletion from the C-terminus of WA PKS caused product change to heptaketide isocoumarins from heptaketide naphthopyrone, YWA1 1, a product of intact WA PKS." The abstract also says the authors previously expressed wA in a heterologous fungal host and identified a heptaketide naphthopyrone synthase. This sentence names the protein and the product, so it meets the strict reading of rule 1 (PROPOSED reading; the owner decides).
- PMID 11040426 (Alb1): "The Alb1p PKS product in this heterologous expression system was identified as heptaketide naphthopyrone instead of pentaketide T4HN." The abstract also says Alb1p "shows higher sequence similarity to a naphthopyrone synthase than to a T4HN synthase". The abstract does not name that naphthopyrone synthase.
- PMID 11350964: Ayg1p "enzymatically shortened the heptaketide product of Alb1p to 1,3,6,8-THN".
- PMID 21277986 (A. niger): FwnA (An09g05730) complements the fawn mutant; FwnA, OlvA and BrnA are "the respective homologs of alb1/pksP, ayg1 and abr1"; naphtho-γ-pyrones depend on fwnA.
- PMID 31067027 (A. niger and other Aspergilli): biaryl naphthopyrone dimers; the tailoring genes are in a cluster separate from the PKS gene. The abstract does not name albA (the UniProt record cites this paper for the albA function).
- PMID 38392776 (A. nidulans, 2024): the pigment backbone "is a polyketide synthesized by a polyketide synthase encoded by the wA gene"; the yA laccase modifies it; chaA (a P450, AN10028) and fwA (an EthD domain protein, AN1088) are further spore colour genes. The abstract does not name an Ayg1, Arp1 or Arp2 homologue.
- PMID 1465094 (A. nidulans, 1992): wA is needed for the green conidial pigment; wA mutants have white conidia.

### 1c. What this means for the step

- The step `heptaketide_pks` is defined by its product. Its YAML name is "heptaketide naphthopyrone YWA1 synthase" and its only reaction is RHEA:62652 (MEASURED, `pathways/melanin_dhn.yaml`). The YAML `lineage_note` already says: "Aspergillus nidulans wA (Q03149) does the same reaction (YWA1) for a pigment that is not DHN-melanin, so a hit in this clade does not show DHN-melanin without Ayg1 and the downstream steps" (MEASURED).
- wA does the same biochemical step as Alb1 and AlbA (REPORTED: PMID 11251292, 11040426; UniProt RHEA:62652 on all three). The difference is downstream: Alb1 feeds Ayg1 and the DHN route (REPORTED, PMID 11350964); wA feeds the yA laccase route (REPORTED, PMID 38392776).
- Therefore a "decoy" label for wA on the heptaketide step would be a misclassification. The curated files do not make that error (section 2). The joint tree does, by its union rule.
- Whether wA is an ortholog of Alb1 is UNVERIFIED (no species-tree reconciliation was done). Section 4 gives MEASURED data that fit one copy per genome that follows the species.
- That the A. nidulans conidial pigment is not DHN melanin is stated by the YAML note and by the decoy row of yA P17489 (curator text, 2026-11-21). No abstract read here says it in those words (UNVERIFIED as a literature statement).

## 2. Why wA is a decoy, and whether the reason holds

### 2a. The provenance (MEASURED)

- `refs/melanin_dhn/decoys/provenance.tsv` lists Q03149 twice: for `thn_synthase` and for `hexaketide_pks`. The reason text is: "UniProt reviewed entry: heptaketide naphthopyrone YWA1 synthase (RHEA:62652) of the conidial yellow pigment; not a T4HN or 2-acetyl-T4HN synthase. **It is NOT a decoy of heptaketide_pks: it does that reference reaction**."
- `refs/melanin_dhn/decoys/heptaketide_pks.faa` holds 4 decoys: aflC Q12053, stcA Q12397, CTB1 Q6DQW3, bik1 S0DZM7. It never held wA: the two commits that touch it (5832a03, 9e02589) both have the same 4 ids (MEASURED, `git log`, `git show`).
- `pathways/melanin_dhn.yaml` lists wA as a decoy under `thn_synthase` and `hexaketide_pks`, not under `heptaketide_pks` (MEASURED).
- `kb/melanin_dhn/` has no row for wA (MEASURED, grep).

So the reason "wA makes YWA1, not T4HN and not 2-acetyl-T4HN" holds for the two other PKS steps. It does not apply to the heptaketide step, and the curator wrote that down.

### 2b. How wA still came to bound the heptaketide step (MEASURED)

1. Step 2 scored every decoy of the pathway against every HMM. `results/step2/melanin_dhn/target_scores_at_min.PROVISIONAL.tsv` has wA on the heptaketide HMM with `target_step` = `hexaketide_pks` and `thn_synthase` (3,904.0 against `min_score` 3,656.7). The step-2 note (`docs/superpowers/notes/2026-11-23-melanin-mycosporine-step2.md`, line 48) says Q03149 is "in three PKS decoy files". That is not correct: it is in two.
2. Step 5 joint: `inputs/decoys_all.faa` is "the union of the three decoy files" (`results/step5/melanin_dhn/pks_joint/RUN.md`). `define_joint_targets` in `bin/step5_joint.py` builds one set `labelled = decoys | outgroups | members` and uses it as the boundary for every step. A decoy of any step is a boundary for all steps.
3. `run_clades` stops if one leaf has two kinds ("a leaf has two kinds (reference, decoy, named member, outgroup)"). So wA cannot be a heptaketide reference and a thn/hexa decoy in the same run without a code change.

### 2c. Domains and distances (MEASURED)

- Domain order on all three records: SAT (8-244), KS3, MAT, PT (about 1,290-1,603), Carrier 1, Carrier 2, Claisen cyclase (C-terminal; "TE" in older papers). Pfam and NCBIfam sets on the records: Alb1 and AlbA have the same set (TIGR04532, PF00698, PF00109, PF02801, PF00550, PF16073, PF00975). wA has this set plus PF22621. The Alb1 DOMAIN text says "an acyl-carrier protein" but its features list Carrier 1 and Carrier 2.
- Global identity (BLOSUM62, gap open -10, extend -0.5; full-length sequences from `results/step5/melanin_dhn/pks_joint/inputs/`): wA–AlbA 69.9 % (1,517 of 2,170 columns); wA–Alb1 67.7 % (1,468 of 2,167); Alb1–AlbA 69.9 % (1,514 of 2,167). For comparison: P9WET9 to the three 65.0 % to 66.2 %; I1RF58 46.4 % to 46.8 %; F9XMW3 44.8 % to 45.5 %; stcA (decoy) 39.1 % to 39.5 %.
- Patristic distance on the joint tree (`joint_report/labelled_distances.tsv`): wA–AlbA 0.43, wA–Alb1 0.54, Alb1–AlbA 0.45; P9WET9 0.53 and 0.60; I1RF58 1.33 and 1.40; F9XMW3 1.94 and 2.00 substitutions per site.
- Topology (tree rooted on decoy bik1 S0DZM7 for display): ((AlbA, (wA, an Aspergillus filifer representative)), (Alb1, two A. fumigatus representatives)). **wA is the sister of an Aspergillus representative and these two are the sister of AlbA.** wA is not an outlier of the group. It sits inside it.

The HMM score and the tree both put wA inside the heptaketide group. Neither separates it.

## 3. Options

Each option states cost, what changes in the census, risk, and the owner decisions involved. Counts are MEASURED from `results/step5/melanin_dhn/pks_joint/` unless marked.

### A. Make wA a third heptaketide reference

- What changes: `refs/melanin_dhn/heptaketide_pks.faa` and `provenance.tsv` get Q03149. A new `kb` row (RHEA:62652, PMID 11251292). The heptaketide HMM is rebuilt. `min_score` becomes the lowest leave-one-out score of 3 references (owner rule 3). Step 4 full search for heptaketide is re-run. The joint tree is re-run (wA becomes a reference leaf).
- Cost: step 2 (minutes), step 4 heptaketide search over 23,376 genomes (the full step 4 ran in 6 chunks for three pathways; a heptaketide-only re-run time is UNVERIFIED), step 5 prepare and place (49 min 37 s and 42 min 52 s last time). Code: `run_clades` rejects a leaf that is a reference of one step and a decoy of another, so either wA leaves the thn/hexa decoy files or the tool must accept two roles.
- Effect on other PKS steps: thn and hexa keep their references. If wA is dropped from their decoy files, their `min_score` floors change: both floors are set by wA (thn 1,908.6, hexa 1,839.0; step-2 note). Their candidate sets and targets would then change too. This is a real cost of option A.
- Risk: the new `min_score` is unknown (UNVERIFIED). wA scores 3,904.0 on the current HMM, above the current LOO floor 3,656.7, so a 3-reference LOO floor may rise or fall. The candidate set may differ from 2,200.
- Owner decisions: rule 1 strict (the PMID 11251292 sentence names WA PKS and YWA1); gene-name rule is not involved (`wA` is the UniProt gene name); whether a reference may have a pigment pathway other than DHN melanin (AlbA already does, section 1a).

### B. Step-specific decoys in the joint tree: wA labelled inside the heptaketide target (recommended)

- What changes: the joint clade rule uses, for step S, the decoys of S plus the named leaves plus the references of the other steps. Decoys of other steps are not boundaries for S. For this pathway that differs from today only by wA: the thn and hexa decoy files are the 4 heptaketide decoys plus wA (MEASURED, decoy files). So the result equals `variant_wA` (`--include heptaketide_pks=Q03149`).
- wA keeps its role as decoy of thn and hexa. Their targets do not change (99 and 27 leaves). Target overlap stays 0 (`variant_wA/target_overlap.tsv`).
- Result (MEASURED, `variant_wA/`): target 7 leaves (Alb1, AlbA, wA, 4 representatives: Aspergillus filifer, 2 A. fumigatus, Evansstolkia leycettana). 2,123 of 2,200 proteins assigned, 77 unassigned, 0 on the other side; 2,091 genomes. Minimum LWR sum in target 0.986.
- Boundary of the target: P9WET9 (Paecilomyces variotii, a named leaf with YWA1 on its record). It is the next labelled leaf when the side grows from 7 to 8 leaves (MEASURED, sides of the tree). So the boundary is another named YWA1 synthase, not a non-YWA1 NR-PKS. The 77 unassigned proteins are the Paecilomyces (72) and Thermoascus (5) proteins near P9WET9 (MEASURED, section 4).
- wA could also be moved to a `tree_members` file with a role such as "same reaction, other pigment", as for the carotenoid named members (`refs/carotenoid/tree_members/README.md`). For the clade rule this is the same as B. The README states that a named-member option in the tools is an owner decision.
- Cost: a small change to `define_joint_targets` (per-step decoy sets) plus a test; re-run `clades`, `assign`, `filter-placements`, `assign_steps` from the saved `epa_result.jplace.gz`. No HMM, search or tree re-run. The `variant_wA/` folder already holds these outputs.
- Risk: the target rule is changed for a reason found after the run. State that in the report. The 7-leaf target depends on one tree with no bootstrap.
- Owner decisions: the joint clade rule (per-step decoys instead of union decoys); whether a target may hold a labelled leaf of another pigment pathway.

### C. Widen to all leaves with YWA1 on the record (38-leaf variant), bounded by an outer NR-PKS group

- Result (MEASURED, `variant_YWA1/`): 38 leaves (2 references, wA, P9WET9, I1RF58, 33 representatives). All 2,200 proteins assigned, 2,168 genomes.
- Boundary: when the side grows from 38 to 45 leaves it adds 7 labelled leaves and no representative: pgmA A0A1W7M1U5, elcA Q0UI00, aflC Q12053, stcA Q12397, CTB1 Q6DQW3, fsr1 S0DTP6, bik1 S0DZM7 (MEASURED). These are the nearest non-YWA1 NR-PKS on this tree. Whether they are a good outgroup for a "YWA1-type" group in a larger tree is UNVERIFIED.
- Problem 1: "has YWA1 on the record" is not a clade on this tree. F9XMW3 (Zymoseptoria PKS1) has RHEA:62652 on its record, but by similarity, with no cited PMID (MEASURED). It lies 1.94 and 2.00 from the references and is outside the 38 leaves.
- Problem 2: when the side grows from 13 leaves (the side that first holds I1RF58, plus one representative) to 38 leaves, it adds 25 representatives and no labelled leaf. They are mostly Lecanoromycetes, with Leotiomycetes, Dothideomycetes and other classes (MEASURED, sides of the tree and subtree listing). Only the 2,200 proteins with a heptaketide HMM row are called, and all extra calls are 77 Eurotiomycetes proteins. A later change of the HMM cutoff would let the lichen and other proteins in.
- Cost: as B, plus the owner must accept I1RF58 (aurofusarin, Fusarium) and P9WET9 inside the target.
- Owner decisions: as B; and whether the step is "any YWA1 synthase" or "YWA1 synthase of the Aspergillus pigment type".

### D. Split the step in two (YWA1 for DHN melanin; YWA1 for other pigments)

- Not supported. The three records and abstracts give the same product, YWA1 (section 1). A split by product has no basis.
- A split by downstream route belongs to the route logic (Ayg1 and the reductases), not to the PKS step. The YAML already does this: route `pks_heptaketide` needs `ayg1`, `t4hn_reductase`, `scytalone_dehydratase`, `t3hn_reductase` after `heptaketide_pks`.
- Cost if done anyway: two steps, two HMMs from 2 and 1 references, a single-reference floor for the wA step. The tree would still not separate them by a gap (section 2c). PROPOSED: do not do this.

### E. Leave wA out and report heptaketide as unresolved

- Cost: none in compute.
- What is lost: the step has no target. 0 heptaketide calls instead of 2,123 (B) or 2,200 (C) proteins; 0 genomes instead of 2,091 or 2,168. Route `pks_heptaketide` cannot be completed in any genome.
- Risk: a known biochemical fact (same reaction) is reported as a failure of the method.

### Recommendation (PROPOSED)

Option B. Reasons:

1. It applies the curator's own record (2026-11-21): wA "is NOT a decoy of heptaketide_pks". The union rule of the joint tree broke that record without a decision.
2. It changes no reference, HMM, `min_score` or search. thn and hexa keep wA as decoy and keep their targets.
3. It keeps a named boundary (P9WET9) and leaves 77 proteins near it unassigned and visible, so the target does not grow past a leaf with a different organism-pathway name.
4. The cost is one rule change and a re-run from the saved placement; the result already exists in `variant_wA/`.

Option A can follow later if the owner wants wA in the HMM. It costs steps 2, 4 and 5 and it moves the thn and hexa floors if wA leaves their decoy files.

### What would change the recommendation

- A bootstrap or second tree in which the 7-leaf group is not stable, or in which P9WET9 falls inside the Alb1/AlbA/wA group. Then C or a new named-leaf set is better.
- A paper that shows wA and Alb1 give different products in the same host. None was found in the abstracts read here.
- An owner decision that the step means "the PKS of DHN melanin" (a pathway role), not "YWA1 synthase" (a reaction). Then the YAML step name and reaction must change, and wA-side genomes need evidence from Ayg1 and the reductases.
- More named YWA1 synthases from Penicillium or Talaromyces (UNVERIFIED which exist in UniProt). 671 genomes (mostly Penicillium) place on the stem of the 6-leaf group (section 4). A named leaf there could split the target.

## 4. Check in the data (MEASURED)

Method: each of the 2,200 proteins with a heptaketide row in `results/search/pks_joint/pks_scores.tsv` gets its best epa-ng edge (highest LWR) from `results/step5/melanin_dhn/pks_joint/epa_result.jplace.gz`. The edge side that does not hold S0DZM7 is compared with these leaf sets of `ref.raxml.bestTree`: Alb1 side = {Alb1, 2 A. fumigatus representatives}; wA side = {wA, A. filifer representative}; AlbA leaf; AlbA+wA stem; 6-leaf stem; Evansstolkia / 7-leaf stem; rest of the 38-leaf YWA1 variant. Genome class, genus and species are from `results/manifests/genomes.tsv.zst`. All 2,200 proteins mapped to a genome.

| Sub-region | Proteins | Genomes | Main species (proteins) |
|---|---|---|---|
| Alb1 side | 448 | 448 | A. fumigatus 392, A. lentulus 12, A. udagawae 6 |
| wA side | 242 | 217 | Aspergillus sp. 63, A. latus 52, A. nidulans 19, A. spinulosporus 11, A. sydowii 10 |
| AlbA leaf | 218 | 218 | A. niger 122, A. tubingensis 23, A. luchuensis 14 |
| AlbA+wA stem | 538 | 537 | A. flavus 324, A. oryzae 112, A. parasiticus 48 |
| 6-leaf stem (Alb1, AlbA, wA) | 674 | 671 | Penicillium 605 (P. chrysogenum 88, P. commune 65), Aspergillus 61 |
| Evansstolkia / 7-leaf stem | 3 | 3 | Hamigera, Evansstolkia, Penicilliopsis (1 each) |
| rest of the 38-leaf variant | 77 | 77 | Paecilomyces 72, Thermoascus 5 |
| total | 2,200 | 2,168 | |

Genome counts per row add to more than 2,168 because 3 genomes have proteins in two sub-regions.

- By class (proteins): Eurotiomycetes 2,192 of 2,200; Dothideomycetes 5; Sordariomycetes 3. Variant wA and variant YWA1 agree with this table: `variant_wA` calls all proteins in the first six rows `target` (2,123) and the 77 in the last row `unassigned`; `variant_YWA1` calls all 2,200 `target`.
- One protein per genome is the rule: 2,200 proteins in 2,168 genomes.
- Species follow the sub-regions. A. nidulans: 19 of 19 genomes on the wA side only. A. fumigatus: 392 of 393 on the Alb1 side, 1 on the 6-leaf stem. A. niger: 122 of 122 on the AlbA leaf.
- **Genomes with a candidate on the wA side and none on the Alb1 side or AlbA leaf: 216 of 217 wA-side genomes** (214 Aspergillus, 1 Zalaria, 1 Penicillium). One genome (Aspergillus) has candidates on both the wA side and the Alb1 side. 447 genomes have an Alb1-side candidate and no wA-side candidate.
- Under the default rule these 2,200 proteins have no call. Under B, the 216 genomes above get their target call from a wA-side protein. In 215 of them it is the only heptaketide candidate; 1 also has a 6-leaf-stem candidate.
- Aspergillus genomes with proteins: 1,594; with a heptaketide candidate: 1,466 (92.0 %). Top species without one: A. terreus 24, A. fumigatus 23, A. oryzae 17. The YAML says A. terreus lacks the DHN-melanin genes (PMID 27133313, not re-read here). Why 23 A. fumigatus genomes have no candidate is UNVERIFIED (gene-model or assembly gaps are possible causes; not checked).

Downstream HMM rows in the same genomes (MEASURED, `results/search/melanin_myco_full/`; these are HMM rows at `min_score`, not tree calls, so they do not show orthology):

| Sub-region | Genomes | ayg1 | t4hn_reductase | scytalone_dehydratase | mco_abr1 |
|---|---|---|---|---|---|
| Alb1 side | 448 | 448 | 448 | 447 | 448 |
| wA side | 217 | 195 | 96 | 137 | 98 |
| AlbA leaf | 218 | 218 | 23 | 23 | 217 |
| AlbA+wA stem | 537 | 537 | 0 | 2 | 0 |
| 6-leaf stem | 671 | 667 | 664 | 661 | 660 |

A. nidulans (19 genomes): ayg1 row 19 (scores 501.9 to 644.8; `min_score` 374.2; A. fumigatus 681.3 to 807.6), t4hn_reductase 0, mco_abr1 0, scytalone_dehydratase 8. This fits the reported A. nidulans route (wA, yA, chaA, fwA; PMID 38392776) and not the full DHN route. Whether the A. nidulans ayg1 row is an Ayg1 orthologue is UNVERIFIED (no ayg1 tree placement was read).

This table also shows why the PKS step alone does not call DHN melanin: AlbA-leaf and AlbA+wA-stem genomes (A. niger, A. flavus group) mostly lack t4hn_reductase rows (23 of 218 and 0 of 537). The route logic, not the PKS target, has to carry that.

## 5. Negative and unexpected results

- The step-2 note says wA is in three PKS decoy files. It is in two (section 2b).
- The joint union rule made a curated non-decoy into a boundary. The default rule then gave the step no target. This was not a property of the data.
- "YWA1 on the UniProt record" does not define a clade on this tree (F9XMW3 is outside, by similarity only).
- No literature read here separates the product of wA from that of Alb1 or AlbA.

## 6. Owner decisions involved

1. Joint clade rule: union decoys (today) or per-step decoys (option B). Code change in `bin/step5_joint.py`.
2. Definition of the step: by reaction (YWA1, RHEA:62652, as the YAML says today) or by pathway role (DHN-melanin PKS).
3. If option A: rule 1 strict reading for the wA link (PMID 11251292 sentence); `min_score` by 3-reference LOO (rule 3); whether wA leaves the thn and hexa decoy files (this moves their floors) or the tool accepts two roles for one leaf.
4. If wA goes to a named-member file: the named-member option in the tools (same open decision as `refs/carotenoid/tree_members/README.md`).
5. Pathway-level `verified: false` stays in all options. The gene-name rule is not involved (`wA` is the UniProt gene name).

## 7. Not done here

- No bootstrap or second tree. The test that would settle stability: 100 bootstrap replicates (or a second tree with 0.8 identity, 532 representatives, drawn to 300) and the support of the 7-leaf and 38-leaf groups.
- No tree with more named YWA1 synthases (for example from Penicillium or Talaromyces). Test: add the reviewed UniProt entries with RHEA:62652 that pass the alignability screen, and see whether the 671 6-leaf-stem genomes split.
- No re-run of step 2 with wA as a reference (option A). Test: build the 3-reference HMM, compute LOO scores, compare with 3,656.7 and with the thn and hexa floors.
- No ayg1 or reductase tree placement for wA-side genomes.
- No full-text reading. The A. nidulans pigment structure was not checked beyond the abstracts.
- Files that this note found to be in error (the step-2 note) were not edited.
