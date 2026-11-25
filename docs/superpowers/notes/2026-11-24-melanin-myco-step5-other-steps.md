---
title: "Census step 5: tree placement for the non-PKS melanin and mycosporine steps"
subtitle: "12 of 15 steps placed, 3 stopped, every target is a clade bounded by labelled leaves"
date: 2026-11-24
---

::: {.callout-important}
## Status and limits

- All results are PROVISIONAL. The three pathway files have `verified: false`. Nothing here was checked at the bench or against known orthologs. No biology claim is made.
- **A `target` call is not an ortholog call.** `target` is the largest tree side that holds all references of the step and no labelled leaf (decoy, outgroup or named member). It depends on which labelled leaves were put in the tree. For 7 of 12 placed steps the target holds every candidate or all but 4 or fewer (`laccase_diphenol`, `hmga`, `mco_abr1`, `ayg1`, `ddgs`, `t3hn_reductase`, `tyrosinase`; the last 21 of 50,447; section 3). Those placements show that the candidates lie near the references; they do not separate step enzymes from other family members.
- Floors are the owner floors of step 2 (`omt` 100.6). No floor was lowered. No placement floor was set; no candidate row was dropped on purpose. 12 candidate rows were lost in domain extraction (section 6).
- One tree per step, one seed, 2 starting trees (`pars{2}`), no bootstrap. "Support" in this note is the placement weight (LWR) of the epa-ng call, not tree support.
- The trees hold at most 300 representatives. For `tyrosinase` and `hppd` the representatives cover only part of the candidates (sections 3 and 5).
- The leaves were chosen by a rule from reviewed UniProt entries (section 2). The function of a leaf is the UniProt annotation; it was not checked in the literature.
- The genome table is not corrected for the 5 re-annotated genomes (step 4 note, section 2). Genome counts use the locus-tag prefix.
- Labels: MEASURED is read from a file named in the text. PROPOSED is a choice or an estimate. UNVERIFIED is not checked.
:::

Record with commands, jobs, parameters and failures: `analysis/pathway_census/results/step5/other_steps_RUN.md`. Tables: `analysis/pathway_census/results/step5/other_steps_summary.tsv` and the files beside it. The PKS steps (`thn_synthase`, `hexaketide_pks`, `heptaketide_pks`) are not in this note.

## 1. Result in one table

Steps: 15. Placed: 12 (`t4hn_reductase` and `t3hn_reductase` on one tree). Stopped: 3 (`scytalone_dehydratase`, `omt`, `atp_grasp_ligase`; section 4). All numbers MEASURED (`other_steps_summary.tsv`). "Refs clade free" is `True` when the references form a clade free of labelled leaves. Leaves = references + decoys + named members + outgroups + representatives. "Other side" = calls on the non-target side; "unassigned" = LWR below 0.5.

| Step | Floor | Identity | Candidates placed (of rows) | Leaves (refs + decoys + named + outgroups + reps) | Refs clade free | Target size (smallest rule) | Target / other side / unassigned | Assigned proteins | Genomes |
|---|---:|---:|---:|---|---|---:|---|---:|---:|
| `aspulvinone_e_synthase` | 1205.5 | 0.95 | 168 of 168 | 55 (1+3+8+5+38) | yes | 18 (1) | 83 / 85 / 0 | 83 | 49 |
| `laccase_diphenol` | 1127.7 | 0.95 | 402 of 402 | 25 (3+1+10+6+5) | yes | 8 (8) | 402 / 0 / 0 | 402 | 215 |
| `hmga` | 983.2 | 0.95 | 2,313 of 2,313 | 76 (2+0+0+1+73) | yes | 75 (65) | 2,313 / 0 / 0 | 2,313 | 2,280 |
| `asp_melanin_tyrosinase` | 320.5 | 0.8 | 2,471 of 2,471 | 112 (1+4+2+0+105) | yes | 5 (1) | 42 / 2,429 / 0 | 42 | 42 |
| `hppd` | 164.3 | 0.7 | 20,098 of 20,103 | 305 (1+1+0+3+300) | yes | 38 (1) | 5,850 / 14,216 / 32 | 5,850 | 5,685 |
| `tyrosinase` | 63.6 | 0.5 | 50,447 of 50,452 | 310 (7+0+0+3+300) | yes | 307 (306) | 50,426 / 0 / 21 | 50,426 | 14,212 |
| `laccase_abr2` | 464.3 | 0.7 | 9,631 of 9,631 | 267 (1+1+5+0+260) | yes | 82 (1) | 3,809 / 5,822 / 0 | 3,809 | 3,411 |
| `mco_abr1` | 638.2 | 0.8 | 1,918 of 1,918 | 124 (2+2+10+6+104) | yes | 106 (26) | 1,918 / 0 / 0 | 1,918 | 1,897 |
| `ayg1` | 374.2 | 0.7 | 11,835 of 11,836 | 283 (3+4+6+6+264) | yes | 267 (266) | 11,835 / 0 / 0 | 11,835 | 8,525 |
| `t4hn_reductase` (joint SDR tree) | 356.0 | 0.8 | 10,010 of 10,010 | 281 (5+4+12+10+250) | yes | 192 (192) | 8,457 / 666 / 887 | 8,457 | 7,852 |
| `t3hn_reductase` (joint SDR tree) | 550.4 | 0.8 | 2,240 of 2,240 | 281 (3+4+14+10+250) | yes | 5 (4) | 2,236 / 0 / 4 | 2,236 | 2,199 |
| `ddgs` | 463.5 | 0.8 | 8,865 of 8,866 | 156 (3+6+10+6+131) | yes | 133 (49) | 8,864 / 1 / 0 | 8,864 | 8,760 |

- Genomes: of 23,376 searched. Target size under the largest-side rule, with the default smallest-side rule in brackets. A smallest-side target of 1 to 5 leaves (for example `aspulvinone_e_synthase` 1, `laccase_abr2` 1, `hppd` 1) is the reference stem: no candidate could be called `target` under that rule. That is why the largest-side rule was used (section 2).
- Support in the `target` call (LWR): minimum 0.97 or more for 8 steps; lower minimum for `tyrosinase` 0.794, `hppd` 0.531, `t3hn_reductase` 0.508, `ddgs` 0.95 (`other_steps_summary.tsv`, `support_min`).
- Long pendant branches (flag only): 1,526 of 50,447 for `tyrosinase`, 15 of 20,098 for `hppd`, 0 for the other 10 steps.

## 2. Method

1. **Floors.** `min_score` of `results/step3/melanin_myco/min_scores.combined.omt100.6.PROVISIONAL.tsv`. Candidates = all rows of the step in `results/search/melanin_myco_full/`. The distribution of each step is in `results/step5/other_inputs/candidate_score_distribution.tsv`. No valley was found near any floor, and no step has more than 50,452 rows, so no placement floor was needed (rows dropped: 0). PROPOSED, UNVERIFIED that any floor keeps all true orthologs.
2. **Target rule.** `bin/step5_multi_ref_clades.py` (new, tests in `tests/test_multi_ref_clades.py`): the rule of `step5_single_ref_clades.py` for one or more references. Owner-accepted PROVISIONAL for one reference (2026-11-23). I use it for steps with 2 to 7 references too. Reason: the default smallest side gave 1 to 5 leaves for 5 of 12 steps (the table above), which puts only near-identical candidates on the stem.
3. **Labelled leaves.** For each step: the references; the decoys of `refs/<pathway>/decoys/` that pass the alignability rule (best domain 10 bits or more and 50% of the HMM or more; `bin/filter_decoys_for_tree.py`, `best` mode); named members and outgroups from reviewed UniProt entries of the Pfam families of the references. Pfam ids: PF06500 (ayg1), PF02982 (scytalone dehydratase), PF13561 (the two SDR steps), PF00394 (three laccase and MCO steps), PF00264 (two tyrosinase steps), PF04209 (hmga), PF00903 (hppd), PF00501 (aspulvinone), PF01761 (ddgs), PF01596 (omt), PF02655 and five more ATP-grasp families (atp_grasp_ligase). Every entry was fetched from `rest.uniprot.org` (reviewed entries, fungal and non-fungal sets; table and FASTA with the same query; the lengths in the two agree for all entries; 16 Pfam families, 32 sets, `fetch_log.tsv`). Organism, gene, name, length and EC in the tables below are from the fetched table.
   - Screen (`bin/screen_outgroups.py`): best domain 10 bits or more and 50% of the HMM or more. The tree uses `best` mode. `union` mode counts are in `leaf_summary.tsv` and are not used in a tree (not run).
   - Excluded from the leaf pool (`bin/select_family_leaves.py`, tests in `tests/test_select_family_leaves.py`): any reference; any entry with a best score at or above the step floor (the search would call it a candidate); any entry with a gene name of a reference of the step (for example Arp1 of Metarhizium for the scytalone dehydratase references); and, for `hppd`, `t4hn_reductase`, `t3hn_reductase`, `scytalone_dehydratase`, `ddgs` and `hmga`, any entry with the same complete EC number as a reference. The EC rule is not used for the laccase, tyrosinase and `omt` steps, because their EC numbers cover the enzyme class of the family (laccase 1.10.3.2, tyrosinase 1.14.18.1) or come from automatic annotation (`omt` reference B2J6X8: 2.1.1.104, UNVERIFIED).
   - Named members = up to 10 fungal entries, best scores first, one per organism. Outgroups = up to 6 other entries at evenly spaced score ranks, non-fungal first. References of the sibling steps that pass the screen are added as named members (`mco_abr1`: Q55P57 and J9VQZ4 of `laccase_diphenol`; `laccase_diphenol`: A2R3T0 of `mco_abr1`). All leaves are cut to the domain (`bin/cut_leaves_to_domain.py`, pad 20, tests in `tests/test_cut_leaves_to_domain.py`).
4. **Representatives.** mmseqs2 identity = the lowest of 0.5, 0.7, 0.8, 0.9, 0.95 that gives at least 100 clusters; 0.95 when none does (`aspulvinone_e_synthase` 38 clusters, `laccase_diphenol` 5, `hmga` 73). More than 300 clusters are capped at 300 (stratified by phylum). PROPOSED. Cluster counts at each identity are in the run record. At 0.5, 8 full runs gave 1 to 32 clusters (kept as `<step>_id50`).
5. **Placement.** `step5_prepare_ident.sbatch` (mafft, clipkit, raxml-ng LG+G4, 4 threads, `pars{2}`), `step5_multi_ref_clades.sbatch`, `step5_place.sbatch` (epa-ng, `upo_place assign`, `assign_steps` with minimum support 0.5). `PC_CAND_DIR` was exported for both parts. After each run the rows of `candidates_step.tsv`, `placements.tsv` and the extraction summary were compared (section 6).
6. **SDR joint tree.** One tree built from the union of the candidates of `t4hn_reductase` and `t3hn_reductase` (10,010 + 2,240 rows; 583 proteins hit both; 11,667 proteins). Each protein is cut at its best domain across both HMMs. The references of the other step are named members of each step. Each step has its own target.

## 3. What the placement shows, step by step

For every step: the target is the clade bounded by the labelled leaves listed in `other_steps_leaves.tsv`. Nearest-leaf summary: `other_steps_nearest_summary.tsv` (median distance in substitutions per site from each representative to its nearest labelled leaf). Assigned genomes by class: `other_steps_assigned_by_class.tsv`.

**Placements where almost every candidate is called target (no separation).**

- `ayg1`: 11,835 of 11,835 `target`. All 264 representatives are nearest to a reference (median distance 1.39). The 6 outgroups are bacterial esterases (FrsA) and the 4 decoys are alpha/beta hydrolases; the 6 named members are fungal cluster hydrolases. None attracts a candidate. Genomes with an assigned protein by class: Dothideomycetes 2,691 of 2,724, Eurotiomycetes 2,478 of 3,201, Sordariomycetes 2,744 of 5,981, Leotiomycetes 310 of 429, Lecanoromycetes 258 of 478, no Saccharomycetes, Pichiomycetes, Agaricomycetes or Tremellomycetes (0 of 7,889). The result says the candidates sit near the three references, with a family not found in the yeasts and basidiomycetes. It does not show that they have the Ayg1 role.
- `hmga`: 2,313 of 2,313 `target`. The only alignable non-HGD leaf is one archaeal protein of unknown function (Q97Y16, 36.5 bits, coverage 0.596); the human and fungal homogentisate dioxygenases were excluded as the same enzyme. So the tree has one outgroup and the call carries no information (section 4).
- `mco_abr1`: 1,918 of 1,918 `target` (104 representatives; the 2 decoys and 6 outgroups and 10 named leaves sit outside). Eurotiomycetes 1,638 of 3,201 genomes, Lecanoromycetes 207 of 478, Dothideomycetes 35 of 2,724.
- `laccase_diphenol`: 402 of 402 `target`; all in Cryptococcus (Tremellomycetes, 215 of 556 genomes; same as the search). The 402 domains give 5 clusters at 0.95 identity, so the tree holds 5 representatives, all nearest to a reference (median distance 0.015).
- `ddgs`: 8,864 of 8,865 `target` (1 on the decoy side). 130 of 131 representatives nearest to a reference (median 2.29). Classes with assigned genomes: Sordariomycetes 4,454 of 5,981, Dothideomycetes 2,539 of 2,724, Lecanoromycetes 416 of 478, Eurotiomycetes 372 of 3,201 and basidiomycete classes (Tremellomycetes 131 of 556, Ustilaginomycetes 118 of 133, Exobasidiomycetes 82 of 91, Agaricomycetes 70 of 1,844). The 6 decoys include the pentafunctional AROM (P07547); the 10 named members are fungal pentafunctional AROM entries (ARO1; 108.1 to 122.4 bits); the 6 outgroups are one 2-epi-5-epi-valiolone synthase (317.1 bits), four bacterial 3-dehydroquinate synthases and one fungal AROM entry.
- `t3hn_reductase`: 2,236 of 2,240 `target`, but the target has 5 leaves (3 references and 2 representatives) and 2,236 of 2,240 candidates are placed inside it. Of the 250 representatives, 52 are nearest to a reference and 196 to a named member (`nearest_summary`). Sordariomycetes 2,187 of 5,981 genomes.
- `tyrosinase`: 50,426 of 50,447 `target`. The tree has 3 outgroups (cat, gorilla and human tyrosinase, 33.7 to 36.4 bits, coverage 0.501, at the edge of the alignability rule) and no decoy: the 5 tyrosinase decoys score 19 to 36 bits and cover less than half the HMM. `long_pendant` is set for 1,526 queries. The 300 representatives are a draw from 2,460 clusters; they cover the clusters of 5,585 candidates (11%). 14,212 of 23,376 genomes carry an assigned protein (Sordariomycetes 5,950 of 5,981; Dothideomycetes 2,713 of 2,724; Orbiliomycetes 67 of 67; Agaricomycetes 1,450 of 1,844; Pichiomycetes 2 of 2,929). The placement does not split the 50,447 proteins into tyrosinase and other members of the family.

**Placements with a split.**

- `t4hn_reductase`: 8,457 `target`, 666 on the decoy side, 887 unassigned (of 10,010). The target has 192 leaves with all 5 references; the 3 `t3hn_reductase` references and 9 other named SDRs are outside. Score separates the calls (median 436.1 for `target`, 363.4 for the decoy side, 358.3 for unassigned; floor 356.0). Sordariomycetes 2,889 of 5,981, Dothideomycetes 2,638 of 2,724, Leotiomycetes 346 of 429, Lecanoromycetes 316 of 478, Eurotiomycetes 1,624 of 3,201.
- `laccase_abr2`: 3,809 `target` and 5,822 on the decoy side (of 9,631). The tree is rooted on the one aligned decoy (P17489); the target holds the reference and 81 representatives. Median score 635.6 (target) against 510.5 (decoy side). Eurotiomycetes 2,153 of 3,201 genomes, Sordariomycetes 1,110 of 5,981, Lecanoromycetes 76 of 478.
- `hppd`: 5,850 `target`, 14,216 decoy side, 32 unassigned (of 20,098). Median score 721.2 (target) against 534.3 (decoy side). Eurotiomycetes 3,167 of 3,201, Dothideomycetes 2,131 of 2,724, Pezizomycetes 160 of 171, Sordariomycetes 36 of 5,981. **Reading trap:** the decoy side holds 14,216 candidates and the tree has only 4 labelled leaves (3 bacterial outgroups and the 4-hydroxymandelate synthase decoy O52791). So the target is one lineage of HPPD-like proteins around the *A. fumigatus* reference, not "the HPPDs". 14,216 candidates outside the target do not show that they are not HPPD. The 300 representatives cover the clusters of 13,811 of 20,098 candidates.
- `asp_melanin_tyrosinase`: 42 `target`, 2,429 decoy side (of 2,471). The target has 5 leaves. The floor is the score of the decoy ivoB (320.5), so the candidates include many ivoB-like proteins; the 42 target proteins score 466.3 to 878.6 (median 833.6) and the decoy-side calls score 320.5 to 422.8 (median 346.9). All 42 are in Eurotiomycetes (42 of 3,201 genomes). 56 of 105 representatives are nearest to the reference (median distance 1.85) and the median distance to the nearest labelled leaf is 1.2 to 1.8 for all roles, so the tree is long-branched.
- `aspulvinone_e_synthase`: 83 `target`, 85 decoy side (of 168). Median score 1,455.8 against 1,285.1. 49 genomes, all Eurotiomycetes (49 of 3,201). 3 decoys (the *A. terreus* and *Tapinella* atromentin synthetases and the *A. nidulans* terrequinone-pathway synthetase tdiA) and 8 named NRPS-like synthetases bound the target (18 leaves with 17 representatives). 21 of 38 representatives are nearest to a decoy (median 0.87).

**Overlaps (MEASURED; `tyrosinase_laccase_target_overlap.tsv`, `sdr_target_overlap.tsv`).**

- `tyrosinase` against `laccase_abr2`, `mco_abr1`, `laccase_diphenol` and `asp_melanin_tyrosinase`: 0 proteins shared, in the target sets and in the raw hits. A protein is not assigned to both.
- `laccase_abr2` and `mco_abr1` share 4 proteins in 3 genomes, in the target sets and in the raw hits (the 4 proteins hit both steps). `laccase_diphenol` shares 0 proteins with either.
- Correction of the step 4 note: it says the 4 proteins that hit two of `laccase_abr2`, `laccase_diphenol` and `mco_abr1` are also the only ones that hit `tyrosinase` plus a laccase or MCO step. In the candidate files, 0 proteins hit `tyrosinase` and a laccase or MCO step; the 4 proteins hit `laccase_abr2` and `mco_abr1` (`target_overlap.py` with `--hit-floor 0`).
- `t4hn_reductase` and `t3hn_reductase`: 0 proteins in both target sets (8,457 and 2,236). Of the union of 11,667 candidate proteins, 10,693 are in exactly one target set and 974 are in none. The 583 proteins that hit both HMMs are one tree placement each; a protein can be `target` for one step only, because the two targets are disjoint.

## 4. Steps stopped, and one that carries no information

These are results. No step was forced.

- **`scytalone_dehydratase`: the references do not form a clade free of named members.** Tree of 219 leaves (4 references, 8 named members, no decoy, no outgroup, 207 representatives at identity 0.8). The smallest side that holds all 4 references has 97 leaves and holds one named member: A0A2I2F2K6 (*Aspergillus candidus* cfoI, a "Dehydratase" of the chlorflavonin cluster; 73.9 bits, coverage 0.929; branch length 5.06 substitutions per site). `step5_multi_ref_clades.py` exit 3 (job 32435511). A first run at identity 0.5 (`scytalone_dehydratase_id50`) with 6 Metarhizium Arp1 entries as leaves (3 named members and 3 outgroups) also stopped at exit 3: the references split, with the *A. fumigatus* Arp1 reference O14434 inside the Arp1 group. I then excluded entries with the gene name of a reference (Arp1) as possible orthologs. 6 Arp1 entries and 4 entries at or above the floor were excluded. The decoys (ketosteroid isomerase, SnoaL, AcmA, NUTF2, limonene epoxide hydrolase) score 0 on the HMM, so none is in a tree. Decision for the owner: drop the one long-branch leaf (cfoI) and rerun, or keep the stop (section 8).
- **`omt`: the references do not form a clade free of decoys and outgroups.** Tree of 311 leaves (2 references, 3 decoys, 0 named members, 6 outgroups, 300 representatives at identity 0.7; 420 clusters). The outgroups are one zebrafish COMT, three plant caffeoyl-CoA O-methyltransferases, one *Dictyostelium* entry and one daffodil entry. The smallest side that holds both references has 210 leaves and holds 6 outgroups and 3 decoys (the 3 are P21964 human COMT, O42898 *S. pombe* COMT1 and Q8NKC1). The two references are the Phaffia OMT and the *Nostoc* MysB (17% identical, step 2 note) and sit far apart. Job 32435495, `rc=3`. The 14,441 candidates were not placed. Whether the fungal reference alone (single-reference rule, the *Nostoc* protein as a named leaf) is a valid target is a curation decision (section 8).
- **`atp_grasp_ligase`: no alignable non-target leaf.** The 6 decoys score 0 on the HMM (0 of 6 pass). The 11 reviewed entries of PF02655 pass 0 (highest score 10.7 bits). 1,100 reviewed entries of five more ATP-grasp Pfam families (PF02222, PF02786, PF07478, PF08443, PF13535) pass 0 (highest 21.5 bits, coverage 0.216). Without a leaf, `target` would be the whole tree (as for `carx`). Not placed; the 9,167 rows stay a search-hit step.
- **`hmga`: placed, but the call carries no information.** See section 3: one archaeal outgroup, all 2,313 candidates `target`. The 85 family entries that pass the screen are excluded: 75 by the gene name `hmgA` (these are homogentisate dioxygenases), 8 by the EC number 1.13.11.5 (same enzyme) and 2 are the references.

## 5. Leaf screen (MEASURED; `leaf_summary.tsv`, `decoy_alignability.tsv`, `leaves_selected.tsv`)

| Step | Family entries screened | Pass `best` / `union` | Passing but excluded | Decoys alignable (`best` / `union`, of all) | Named selected | Outgroups selected |
|---|---:|---|---:|---|---:|---:|
| `ayg1` | 81 | 79 / 79 | 2 | 4 / 4 of 4 | 7 | 6 |
| `scytalone_dehydratase` | 21 | 21 / 21 | 13 | 0 / 0 of 5 | 8 | 0 |
| `t4hn_reductase` | 465 | 438 / 444 | 15 | 4 / 4 of 4 | 10 | 6 |
| `t3hn_reductase` | 465 | 420 / 423 | 14 | 4 / 4 of 4 | 10 | 6 |
| `laccase_abr2` | 186 | 22 / 87 | 16 | 1 / 5 of 5 | 5 | 0 |
| `mco_abr1` | 186 | 40 / 101 | 3 | 2 / 3 of 4 | 10 | 6 |
| `laccase_diphenol` | 186 | 25 / 136 | 7 | 1 / 5 of 5 | 10 | 6 |
| `tyrosinase` | 85 | 13 / 14 | 10 | 0 / 0 of 5 | 0 | 3 |
| `asp_melanin_tyrosinase` | 85 | 7 / 14 | 3 | 4 / 4 of 4 | 3 | 0 |
| `hmga` | 86 | 86 / 86 | 85 | 0 / 0 of 6 | 0 | 1 |
| `hppd` | 195 | 35 / 35 | 32 | 1 / 1 of 7 | 0 | 3 |
| `aspulvinone_e_synthase` | 860 | 28 / 47 | 4 | 3 / 3 of 6 | 10 | 5 |
| `ddgs` | 561 | 550 / 553 | 2 | 6 / 6 of 6 | 10 | 6 |
| `omt` | 90 | 52 / 52 | 6 | 3 / 3 of 5 | 1 | 6 |
| `atp_grasp_ligase` | 11 | 0 / 0 | 0 | 0 / 0 of 6 | 0 | 0 |

- Families: Pfam sets above. "Passing but excluded" = passes `best` and is a reference, scores at or above the floor, has a reference gene name or the same EC (reasons in `leaf_summary.tsv`). After removing the ids that are also decoys, the leaf counts in the trees are those in the table of section 1 (members include the sibling-step references that pass).
- A leaf in the tree is a UniProt reviewed entry (organism, gene, name, length, EC, best score, coverage in `leaves_selected.tsv`). Known weak leaves: `hppd` outgroups are one Vibrio protein named "Hemolysin VllY" (O06695, 137.1 bits), one *Shewanella* "Protein MelA" (P23996, 121.1) and one cyanobacterial "Uncharacterized protein slr0090" (Q55810, 105.9); their relation to HPPD is not stated in the records (UNVERIFIED). `tyrosinase` has only vertebrate outgroups, at the edge of the rule. `omt` has one zebrafish COMT and plant and other eukaryote O-methyltransferases as outgroups.
- `union` mode screens pass many more leaves for the laccase steps (87 to 136 against 22 to 40), because the three cupredoxin domains of the HMM are scattered. No tree was built from `union` passes. A `union` run is not done.

## 6. Failures and problems (MEASURED unless marked)

- **12 candidate rows were not placed.** `extract_failures.tsv`: `protein_file_missing` for prefixes FB9275E9 (ayg1 1, ddgs 1, hppd 1, tyrosinase 5), F7CB232C (hppd 1) and FAD1592B (hppd 1); `bad_coordinates` for F75533D7 and FADCCC61 (hppd 2). The three prefixes are 3 of the 5 new locus tags of the re-annotated genomes in the step 4 note; the 2 `bad_coordinates` rows are in genomes whose protein length no longer matches the row (UNVERIFIED cause: changed protein files). The check "`candidates_step.tsv` rows equal `placements.tsv` rows" fails for exactly these 4 steps; in every case the difference equals the extraction failures.
- **Stalled dependency.** A job with `afterany` on a job that stays `DependencyNeverSatisfied` does not start. Twice I cancelled the stuck job and removed the dependency with `scontrol update JobId=<id> Dependency=` (job 32435498, 32435513).
- **The `short` partition was at its user limit.** The first submission (38 jobs) stayed in `AssocGrpCpuLimit` and was cancelled; the run used `intel` (the user CPU limit on `short` is 384).
- **Round 1 at identity 0.5.** 8 steps gave 1 to 32 clusters; their runs are kept as `<step>_id50` and the placements are listed in `other_steps_summary.tsv` (they are diagnostic results, not the placement). Example: `laccase_diphenol_id50` has 1 representative.
- **Tyrosinase place job took 1 h 17 min (limit 2 h).** `hmmalign` of 50,447 domains is single-threaded. A second run with `hmmalign` split over 15 CPUs (`slurm/step5_place_parallel.sbatch`, job 32438154, 10 min 7 s) gave the same result for 50,446 of 50,447 queries; the one query that differs is `F5E150A2_004610-T1` (target with LWR 1.000 in one run, unassigned in the other; the cause is not known, UNVERIFIED). The main run is the first; the second is `tyrosinase_par`.
- **Test suite:** `analysis/pathway_census/tests`: 924 passed, 1 failed, 1 skipped. The failure is `test_step5_place.py::test_extract_end_to_end_with_failures`; `bin/step5_place.py` is not changed, and the carotenoid single-reference record (888 passed, 1 failed) lists the same test. 23 of the passed tests are new.

## 7. What the placements do and do not show

- Shown (MEASURED): for 8 steps (`laccase_diphenol`, `mco_abr1`, `ayg1`, `ddgs`, `t4hn_reductase`, `t3hn_reductase`, `hmga`, `tyrosinase`) at least 84% of the candidates are on the target side, and for 4 of them (`laccase_diphenol`, `mco_abr1`, `ayg1`, `hmga`) all of them are. The candidates lie close to the references in a tree that has few labelled neighbours.
- Shown: the candidate sets of the two SDR steps are disjoint in the joint tree (0 shared proteins), and the tyrosinase set shares 0 proteins with every laccase or MCO set.
- Not shown: that a `target` protein has the function of the step. `target` is not an ortholog set. A fall in the `decoy` side is not shown to be a non-ortholog (see `hppd`).
- Not shown: any absence. A genome without an assigned protein is not shown to lack the step (floors are PROVISIONAL, representatives are a draw).

## 8. Decisions for the owner

All are PROPOSED. None has data beyond the numbers above.

1. **`scytalone_dehydratase`.** One long-branch named leaf (cfoI, A0A2I2F2K6) breaks the clade of the 4 references. Drop it and rerun, or keep the stop. A run without it is a selection of leaves after seeing the tree.
2. **`omt`.** Use the fungal Phaffia reference alone as target (the *Nostoc* MysB as a named leaf), or keep the stop. The two references are 17% identical.
3. **`atp_grasp_ligase`.** Stays a search-hit step. A placement needs an alignable non-target leaf; none exists in 1,111 reviewed ATP-grasp entries. A fungal ATP-grasp family member (not reviewed) would have to be chosen.
4. **`hmga`.** Keep as "family present". The family is homogentisate dioxygenase, and every reviewed entry that aligns is the same enzyme.
5. **Target rule for steps with more than one reference.** The largest-side rule gave a target of 5 to 18 leaves for 4 steps (`aspulvinone_e_synthase`, `laccase_diphenol`, `asp_melanin_tyrosinase`, `t3hn_reductase`). Accept the rule for steps with 2 or more references, or choose another.
6. **`tyrosinase` and `hppd`.** The 300 representatives cover 11% and 69% of the candidates. Choices: a higher cap (cost: raxml time; the `tyrosinase` tree took 49 min with 310 leaves at 4 threads), or a second draw with another seed. The labelled leaves are the weak part (3 vertebrate outgroups; 3 bacterial outgroups and 1 decoy).
7. **Identity rule.** The lowest identity with at least 100 clusters is a choice. A step whose 402 candidates fall into 5 clusters (`laccase_diphenol`) has a 5-representative tree.
8. **Genome table.** The 12 rows not placed and 13 candidate rows without a genome row need the genome table fix (other agent).

## 9. Not done here

- No `union`-mode tree, no second tree or seed for any step, no bootstrap, no test with known orthologs.
- The PKS steps (other agent). No change to floors, references or decoys. No lowered cutoff. No second round on `omt`.
- No literature check of any leaf. No fusion check. No tree-based check of the 4 proteins that hit two laccase steps.
- The target rule was not tested for stability against the number of labelled leaves.
