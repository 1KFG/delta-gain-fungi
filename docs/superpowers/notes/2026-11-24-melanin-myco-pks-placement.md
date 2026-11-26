---
title: "Census step 5: one joint tree for the three PKS steps of DHN melanin"
subtitle: "thn_synthase, hexaketide_pks, heptaketide_pks on one reference tree"
date: 2026-11-24
---

::: {.callout-important}
## Status and limits

- All results are PROVISIONAL (`pathways/melanin_dhn.yaml` has `verified: false`). Labels: MEASURED = read from a file named here; PROPOSED = a choice; UNVERIFIED = not checked.
- The tree has 264 leaves, of which 237 are cluster representatives at 0.7 identity. The 9,374 candidates were placed on it with epa-ng. They were not put in the tree.
- A target clade is a tree region bounded by the labelled leaves that were put in the tree (4 references, 5 decoys, 18 named leaves). It is not an ortholog set. A different set of named leaves would give a different target.
- **The two references of heptaketide_pks (Alb1 and AlbA) do not form a clade without the decoy wA (Q03149).** The default rule gives that step no target (section 3). The two widened targets in section 4 are PROPOSED variants, not the rule.
- Two raxml-ng starts only. Tree stability (bootstrap) was not measured (UNVERIFIED).
- No outgroup is in the tree. All four outside candidates failed the alignability rule (section 2).
- Nothing here shows the function or the product of any protein.
:::

Record, commands, jobs and checks: `analysis/pathway_census/results/step5/melanin_dhn/pks_joint/RUN.md`. Code: `bin/step5_joint.py`, `bin/pks_joint_candidates.py`, `slurm/step5_prepare_joint.sbatch`, `slurm/step5_place_joint.sbatch`, tests `tests/test_step5_joint.py`.

## 1. Candidate set

- Input: `results/search/melanin_myco_full/` (rows at the `min_scores` of step 4). Rows: thn_synthase 7,286; hexaketide_pks 8,604; heptaketide_pks 2,200 (MEASURED).
- Union of proteins: 9,374. Hit on 1 HMM: 994; on 2: 8,044; on all 3: 336 (MEASURED, `results/search/pks_joint_union/union_summary.json`).
- One row per protein for the tree: the row with the highest raw score over the three HMMs. Best row is thn_synthase for 5,407 proteins, hexaketide_pks for 1,767, heptaketide_pks for 2,200. Raw scores of different HMMs are on different scales (PROPOSED rule, as in `step5_place best_rows`).
- All 9,374 domains were cut (0 failures). Cut length 942 to 2,250 residues, median 2,160.
- The three step tools (`step5_place.py`) make one tree per step. Joint use needed a new script (`bin/step5_joint.py`) for the clade rule, per-step placement and the report. The old tools were not changed. The new sbatch files are copies with `PC_ROOT` and `PC_CAND_DIR` required.

## 2. Leaves of the tree and alignability

| Kind | Leaves | Alignability (MEASURED, `inputs/screen_alignability.tsv`) |
|---|---|---|
| References | P79068 (thn_synthase), Q9Y7A7 (hexaketide_pks), Q4WZA8 and A2QUI2 (heptaketide_pks) | 4,850.3 to 5,390.1 bits on its own HMM; 1,872.0 to 2,099.0 on the other two; best-domain coverage 0.987 to 1.0 |
| Decoys (5) | wA Q03149, aflC Q12053, stcA Q12397, CTB1 Q6DQW3, bik1 S0DZM7 | best score 1,187.0 to 3,904.0 bits; coverage 0.598 to 1.0 on all three HMMs |
| Named leaves (18) | UniProt reviewed fungal non-reducing PKS (list: `inputs/named_leaves.tsv`) | best score 407.8 to 2,230.4 bits; coverage 0.526 to 0.998; `pass_best` on all 3 HMMs |
| Representatives (237) | clusters at 0.7 identity, 0.8 coverage | not screened (they are hits) |

- Every UniProt accession was fetched from rest.uniprot.org on 2026-11-24. For each of the 24 records the FASTA length equals the record length and OX equals the record organism id (MEASURED, `inputs/fetched_records.tsv`). The product text in `named_leaves.tsv` is the text of the UniProt record, not a claim of this note.
- Three named leaves have the YWA1 reaction (Rhea) on the UniProt entry: I1RF58 (PKS12, aurofusarin, Gibberella zeae), F9XMW3 (PKS1, Zymoseptoria tritici; the entry text says DHN-melanin) and P9WET9 (Paecilomyces variotii pigment cluster 24).
- Excluded by the rule: F1DBA9 (mpaC) and Q5ATJ7 (ausA) fail `pass_best` on the hexaketide HMM (coverage 0.407).
- **No outside outgroup passes.** LovB Q0C8M3 (best-domain coverage 0.402 to 0.411), EryA1 Q03131 (0.402 to 0.418), yeast FAS1 P07149 (under 12 bits) and FAS2 P19097 (no hit) all fail best and union mode. The tree is therefore not rooted on an outgroup. The display root is the 5 decoys; the clade rule does not use the root.
- Alb1 (Q4WZA8) is a reference of heptaketide_pks and is not used as an outgroup or named leaf.

## 3. Tree and target clades

- 264 leaves; alignment 3,658 columns after clipkit (10.562% removed), 44.63% gaps; LG+G4; best log-likelihood -630,903.14 (second start; first start -630,909.43); 2 starts.
- Clustering identity: 0.5 gives 21 clusters, 0.6 gives 93, 0.7 gives 237, 0.8 gives 532, 0.9 gives 1,129 (MEASURED). 0.7 was used (PROPOSED; under the cap of 300, so no draw). All 237 representatives are from Ascomycota (the only phylum among the candidates, `reps_by_phylum.tsv`).

| Step | References | Isolated? | Target leaves | Representatives in target | Named leaves / decoys in target |
|---|---|---|---|---|---|
| thn_synthase | 1 | yes | 99 | 98 | 0 / 0 |
| hexaketide_pks | 1 | yes | 27 | 26 | 0 / 0 |
| heptaketide_pks | 2 | **no** | none | none | the smallest side with both references has 6 leaves: Alb1, AlbA, wA (decoy) and 3 representatives |

(MEASURED, `clade_report.json`.) The two targets share no leaf (`target_overlap.tsv`: 0 of 99 and 27). By construction they cannot share a leaf while a decoy exists.

The three references are far apart: patristic distance P79068 to Q9Y7A7 1.90, P79068 to A2QUI2 2.13, Q9Y7A7 to A2QUI2 2.06 substitutions per site (MEASURED, `joint_report/labelled_distances.tsv`). Alb1 and AlbA are 0.45 apart. wA is 0.43 from AlbA and 0.54 from Alb1; the YWA1-entry leaf P9WET9 is 0.53 and 0.60 away.

## 4. Calls per step (MEASURED, `joint_report/calls_per_step.tsv`; `--min-lwr 0.5`, `--min-support 0.5`)

| Step | Proteins with a row | target | other side | unassigned | Assigned (assign_steps) | Genomes with an assigned protein (of 23,376) |
|---|---|---|---|---|---|---|
| thn_synthase | 7,286 | 3,163 | 3,677 | 446 | 3,163 | 3,124 |
| hexaketide_pks | 8,604 | 609 | 7,505 | 490 | 609 | 600 |
| heptaketide_pks, default rule | 2,200 | no target | | | | |
| heptaketide_pks, variant wA | 2,200 | 2,123 | 0 | 77 | 2,123 | 2,091 |
| heptaketide_pks, variant YWA1 | 2,200 | 2,200 | 0 | 0 | 2,200 | 2,168 |

- Variant wA (PROPOSED): the target may hold wA. It has 7 leaves (2 references, wA, 4 representatives). Variant YWA1 (PROPOSED): the target may hold wA, P9WET9 and I1RF58. It has 38 leaves (2 references, wA, 2 named leaves, 33 representatives). Both variants widen the target past leaves that the default rule uses as its boundary, so a call in a variant does not separate the candidate from wA.
- Support (LWR sum in target): thn_synthase minimum 0.981; hexaketide_pks 0.990; variant wA 0.986; variant YWA1 0.996. No placement has a long pendant (0 in all four).
- The 490 unassigned proteins (LWR under 0.5 for every clade) are the same proteins for the two step clade files.

### Overlap of the HMM hits with the tree calls (MEASURED, `joint_report/hit_pattern_vs_target.tsv`)

| HMM hits | In thn target | In hexa target | In neither |
|---|---|---|---|
| thn + hexa (6,316) | 2,741 | 523 | 3,052 |
| thn only (554) | 422 | 7 | 125 |
| hexa only (304) | 8 | 86 | 210 |
| hept + hexa (1,648) | 0 | 0 | 1,648 |
| hept only (136), hept + thn (80) | 0 | 0 | 216 |
| all three (336) | 0 | 0 | 336 |

The row sums add to 9,374. No protein is `target` for two steps (0 of 9,374; also 0 for both variants). With the YWA1 variant, 2,496 proteins fall in the heptaketide target: 2,200 with a heptaketide HMM row and 296 without one (153 thn + hexa, 137 hexa, 6 thn).

So an HMM hit of one PKS step is not a placement in that step's clade: 3,052 of the 6,316 proteins that hit both thn_synthase and hexaketide_pks HMMs lie in neither target.

## 5. Nearest labelled leaves (MEASURED, `joint_report/reference_neighbours.tsv`, `nearest_labelled.tsv`)

- P79068: nearest representatives are Colletotrichum (0.088, 0.093, 0.099), then Verticillium (0.326 to 0.343). Nearest labelled leaves: A7UMW1 (Elsinoe PKS1) 1.73, F9XMW3 (Zymoseptoria PKS1) 1.78, Q9Y7A7 1.90.
- Q9Y7A7: nearest representatives are Exophiala (0.014, 0.137) and Hypoxylon (0.147). Nearest labelled leaves: A7UMW1 1.57, F9XMW3 1.62.
- Q4WZA8 (Alb1): nearest representatives are Aspergillus (0.009, 0.013). A2QUI2 0.45, wA 0.54.
- A2QUI2 (AlbA): nearest representatives are Evansstolkia (0.405), Aspergillus (0.427, 0.463, 0.468).
- Nearest labelled leaf of each of the 237 representatives: P79068 for 92, P9WET9 for 44, A7UMW1 for 43, Q9Y7A7 and Q5BH30 for 26 each, F9XMW3 for 2, Q4WZA8 for 2, A2QUI2 for 1, wA for 1. By candidates in their clusters: P79068 3,161, A7UMW1 2,260, A2QUI2 1,643, F9XMW3 641, Q9Y7A7 616, Q4WZA8 448, P9WET9 345, Q5BH30 137.
- The 24 Dothideomycetes representatives (2,756 candidates) have as nearest labelled leaf A7UMW1 (13), Q5BH30 (4), P9WET9 (4), F9XMW3 (2), P79068 (1). Named leaves lie outside every target by the rule, so these proteins are on the other side of both targets by construction.

## 6. Genomes by class (MEASURED, `joint_report/genomes_by_class.tsv`; denominator: genomes in the class; all classes with at least 20 genomes are in the file)

| Step | Class | Genomes | With a hit | With an assigned call |
|---|---|---|---|---|
| thn_synthase | Sordariomycetes | 5,981 | 2,477 | 2,477 |
| | Leotiomycetes | 429 | 349 | 322 |
| | Lecanoromycetes | 478 | 374 | 308 |
| | Dothideomycetes | 2,724 | 2,550 | 5 |
| | Eurotiomycetes | 3,201 | 740 | 0 |
| | all | 23,376 | 6,535 | 3,124 |
| hexaketide_pks | Eurotiomycetes | 3,201 | 2,350 | 271 |
| | Leotiomycetes | 429 | 344 | 241 |
| | Lecanoromycetes | 478 | 361 | 79 |
| | Dothideomycetes | 2,724 | 2,556 | 5 |
| | Sordariomycetes | 5,981 | 2,160 | 1 |
| | all | 23,376 | 7,819 | 600 |
| heptaketide_pks, variant YWA1 | Eurotiomycetes | 3,201 | 2,160 | 2,160 |
| | Dothideomycetes | 2,724 | 5 | 5 |
| | Sordariomycetes | 5,981 | 3 | 3 |
| | all | 23,376 | 2,168 | 2,168 |
| heptaketide_pks, variant wA | Eurotiomycetes | 3,201 | 2,160 | 2,083 |
| | all | 23,376 | 2,168 | 2,091 |

The 5 genomes whose protein files were rewritten after the genome table was built (step 4 record) are not mapped by prefix. 13 candidate rows are affected; the effect on these counts was not measured.

## 7. Negative and unexpected results

- The default clade rule gives heptaketide_pks no target. wA lies inside the smallest clade of Alb1 and AlbA. The step 2 record already showed that wA scores above the heptaketide cutoff, so no cutoff separates it (3,904.0 against 3,656.7). The tree does not separate it either.
- The hexaketide target (609 assigned proteins in 600 genomes) holds 7% of the 8,604 proteins with a hexaketide HMM row. 7,505 (87%) lie on the other side.
- The thn_synthase and hexaketide targets are small against the hits. In Dothideomycetes 2,550 genomes have a thn_synthase hit and 5 have an assigned call.
- The three PKS HMMs score each other's references 1,872 to 2,099 bits (heptaketide reference on the thn HMM 1,888.5; Q9Y7A7 on thn 2,008.0; P79068 on hexa 2,030.7). The HMM hit sets overlap strongly (8,380 of 9,374 proteins hit two or three HMMs).
- First prepare job (32434626) failed on a path (section 1 of RUN.md). The identity 0.5 run gave 21 clusters, too few.
- The two raxml-ng starts differ by 6.3 log-likelihood units. The topology of the two starts was not compared.

## 8. What the placement does and does not show

It shows, for each protein: the clade of the 264-leaf tree where epa-ng puts the cut domain, with LWR sum at or above 0.5 (all target calls have 0.98 or more). It shows that the thn_synthase and hexaketide_pks references each have a region bounded by the 18 named leaves, 5 decoys and the other references. It shows that Alb1 and AlbA cannot be separated from wA on this tree.

It does not show:

- the product or function of any candidate. Reference identity of a named leaf is a label;
- that a protein outside a target lacks that function. A protein on the "other side" is next to a decoy, named leaf or other-step reference, or in a part of the tree with none of them;
- that the target boundaries are stable. They depend on the 18 named leaves chosen, the 0.7 identity, 2 raxml-ng starts and one alignment;
- a placement of the 9,374 candidates inside the tree. They were placed on the reference tree of 237 representatives, so a candidate next to a representative inherits that region;
- the effect of the rewritten genomes on class counts (UNVERIFIED).

## 9. Not done

- Bootstrap or a second tree with other named leaves.
- A larger identity (0.8, 532 representatives) or the cap of 300.
- A rooted tree on a real outgroup (none passed).
- A check of the 490 unassigned proteins by hand, and a check of why 2,550 Dothideomycetes genomes with a hit are outside both targets.
- A decision on which heptaketide target (none, wA, YWA1) to carry forward. This is the owner's.
- Per-genome tables of the five mapping-affected genomes.


# Update 2026-11-25: identity 0.8, bootstrap, and per-step decoys

::: {.callout-important}
## Status and limits of this update

- All results are PROVISIONAL (`verified: false`). MEASURED = read from a file named here; PROPOSED = a choice; UNVERIFIED = not checked.
- **Bootstrap is not standard.** One standard raxml-ng replicate takes 4.1 to 4.7 h (2 replicates measured). 100 would take about 70 to 80 h with 12 CPUs. The 100 replicates are NNI-only searches (453 s on average). On the same two resampled alignments the NNI tree differs from the standard tree by relative RF 0.113 and 0.129, and has a lower log-likelihood (by 727 and 522). The effect on each support value is UNVERIFIED (section 3).
- **Sampling label: "sampling not yet adequate"**, by the criterion written before the comparison (`SAMPLING_CRITERION.txt`): C1 to C3 pass, C4 fails by a small margin (section 1).
- No outgroup (section 6). The targets are tree regions bounded by labelled leaves, not ortholog sets.
- Nothing here shows the function or the product of any protein.
:::

Record, commands, jobs, times, failures: `analysis/pathway_census/results/step5/melanin_dhn/pks_joint_id0.8/RUN.md`. Code: `bin/step5_joint.py` (decoy scope), `bin/step5_bootstrap.py`, `slurm/step5_{tree,bootstrap,support,rarefaction}_joint.sbatch`, `slurm/step5_reapply_decoy_scope.sh`; tests `tests/test_step5_joint.py`, `tests/test_step5_bootstrap.py`.

**Correction (2026-11-25).** The statement in the step-2 and earlier PKS notes, and in section 3 of this note above, that wA (Q03149) is a decoy of `heptaketide_pks` was wrong. The curated decoy files never held wA for `heptaketide_pks` (`refs/melanin_dhn/decoys/heptaketide_pks.faa` has 4 decoys: Q12053, Q12397, Q6DQW3, S0DZM7; `provenance.tsv` says wA is "NOT a decoy of heptaketide_pks"). wA is a decoy of `thn_synthase` and `hexaketide_pks`. The first run used the union of the three decoy files, so wA bounded the heptaketide target. Older notes are not edited.

## 1. Sampling check before the bootstrap

Criterion (PROPOSED, written first, `SAMPLING_CRITERION.txt`): C1 under 2% of protein-step pairs change target/non-target; C2 at most 10% of old target representatives flip, and the new target is not smaller; C3 references still isolated; C4 every domain has a representative within identity 0.8 and coverage 0.8.

| Check | Result (MEASURED) |
|---|---|
| (a) clusters | 532 clusters for 9,374 domains (0.7: 237). All 532 are leaves; none was drawn (`capped` false). 531 are Ascomycota, 1 has no phylum label (`reps_by_phylum.tsv`). |
| (b) cluster sizes | 263 of 532 are singletons (49.4%); median size 2; largest 774; 81 clusters have 10 or more domains and hold 8,448 of 9,374 domains (90.1%); the 10 largest hold 4,372 (46.6%). |
| (b) representative within the rule | Re-measured by `mmseqs easy-search` of each non-representative domain against its own representative (alignment mode 3). Of 8,842 non-representative domains, 8,642 (97.7%) have identity >= 0.8 and coverage >= 0.8 both ways. 194 have identity under 0.8 (minimum 0.758; 1st percentile 0.798, 5th 0.802; median 0.884). The other 6 fail coverage. The 532 representatives are in their own cluster. The clustering identity of mmseqs and this re-measure are not the same estimate (the reason for the 194 is UNVERIFIED). |
| (c) old calls under the new tree, thn_synthase | 7,286 proteins with a row: target 3,163 in both runs; 0 target to not-target; 0 not-target to target (0.00%). Other classes: 402 proteins change between "other side" and "unassigned" (37 other to unassigned; 365 unassigned to other). |
| (c) hexaketide_pks | 8,604 proteins: target 609 in both; 0 and 0. 407 change between other side and unassigned (37 and 370). |
| (c) heptaketide_pks, per-step rule | 2,200 proteins: target 2,123 in both; 0 and 0; 7 change from unassigned to other side. |
| (c) C1 over the two default steps | 0 of 15,890 protein-step pairs (7,286 + 8,604) change target status. |
| (c) the 237 old representatives as queries | 183 (thn) and 175 (hexa) of them have a row of the step: 0 flips (old target 94 and 25). Of the 56 old representatives that are not representatives now (46 and 48 with a row): 0 flips (old target 19 and 10). 181 of 237 old representatives are also representatives now, so their own call is trivially stable. |
| (c) target size | thn_synthase target 99 leaves (0.7) to 240 (0.8); hexaketide 27 to 56; all old-target representatives that were re-placed stay target. The target grows because the tree holds more representatives, not because the set of proteins in it changes: the assigned proteins are the same numbers (3,163; 609). |
| (d) rarefaction (PROPOSED, approximate) | 50% (266 representatives) and 75% (399) of the 532, 3 seeds each, all labelled leaves kept; the 0.8 tree and alignment pruned, branch lengths of removed leaves summed, model and other lengths not re-optimised; epa-ng, clades by the same rule (145 s per placement). Result: 0 target/non-target changes against the full tree for every step and every one of 6 subsets (denominators 7,286, 8,604, 2,200). The only class changes are "other side" to "unassigned": 80 to 476 proteins for thn_synthase, 80 to 487 for hexaketide (smaller at 75%). Limit: this tests the density of representatives, not the topology, and target calls are bounded by labelled leaves by construction. |

Decision (stated criterion): C1 passes (0 of 15,890). C2 passes (0 flips in 237 old representatives; new targets hold 239 and 55 representatives against 98 and 26). C3 passes. **C4 fails**: 200 of 8,842 non-representative domains (2.3%) lie outside the identity and coverage rule when re-measured (194 on identity, 6 on coverage). By the written rule the label is "sampling not yet adequate". The call stability evidence (C1, C2, rarefaction) is clean, and the C4 failure is at the threshold edge (1st percentile 0.798). The criterion is not changed after the fact. What would improve it: run at identity 0.9 (1,129 clusters at the 0.9 count of the first run) or add the representatives for the 200 domains; the owner can accept the label with C4 noted.

## 2. Tree and targets at identity 0.8 (MEASURED, `clade_report.json`, per-step decoys)

- 559 leaves; 4,118 columns after clipkit (13.25% removed); gaps 49.86%; LG+G4; best log-likelihood -1,098,828.08 (start 1; start 2: -1,098,842.13); two starts only. raxml-ng took 25,358 s.

| Step | Target leaves (0.8) | Representatives in it | Target leaves (0.7, same rule) |
|---|---|---|---|
| thn_synthase | 240 | 239 | 99 |
| hexaketide_pks | 56 | 55 | 27 |
| heptaketide_pks, per-step rule | 22 (Alb1, AlbA, wA, 19 representatives) | 19 | 7 (Alb1, AlbA, wA, 4 representatives) |
| heptaketide_pks, union rule | none | | none |

Per-step targets share no leaf (22 and 56, 22 and 240, 56 and 240: 0 shared). Proteins that are `target` for more than one step: 0.

## 3. Calls per step (MEASURED, `joint_report/calls_per_step.tsv`, `--min-lwr 0.5`, `--min-support 0.5`)

| Step, rule | Proteins with a row | target | other side | unassigned | Assigned | Genomes with an assigned protein (of 23,376) |
|---|---|---|---|---|---|---|
| thn_synthase, 0.8 (0.7) | 7,286 | 3,163 (3,163) | 4,005 (3,677) | 118 (446) | 3,163 (3,163) | 3,124 (3,124) |
| hexaketide_pks, 0.8 (0.7) | 8,604 | 609 (609) | 7,838 (7,505) | 157 (490) | 609 (609) | 600 (600) |
| heptaketide_pks, per-step, 0.8 (0.7) | 2,200 | 2,123 (2,123) | 7 (0) | 70 (77) | 2,123 (2,123) | 2,091 (2,091) |
| heptaketide_pks, union rule, 0.8 (0.7) | 2,200 | no target | | | | |

- thn_synthase and hexaketide_pks give identical `assign` files under the per-step and the union rule, in both runs (cmp, MEASURED): their decoy sets equal the union (5 decoys), so the per-step rule does not change them. The 0.7 per-step result for heptaketide (7 leaves, 2,123 of 2,200) equals the earlier "variant wA" (same `assign` file).
- The three PKS HMM hit patterns against the target are as in the 0.7 run (thn+hexa 6,316 proteins: 2,741 thn target, 523 hexa target, 3,052 neither). With the heptaketide target, 2,126 proteins lie in it: 2,123 with a heptaketide HMM row and 3 without (hexaketide HMM only). Heptaketide genomes by class: Eurotiomycetes 2,083 of 3,201 genomes with an assigned protein; Dothideomycetes 5; Sordariomycetes 3 (`joint_report/genomes_by_class.tsv`).
- Nearest labelled leaves outside the heptaketide target (`bootstrap_report/boundary.heptaketide_pks.tsv`, patristic distance from the nearest target leaf): P9WET9 (named, YWA1 entry) 0.466, I1RF58 (named) 1.250, S0DZM7 (decoy) 1.560, Q12053 (decoy) 1.811, A0A1W7M1U5 (named) 1.848. At 0.7: P9WET9 0.408, I1RF58 1.206, S0DZM7 1.538. Hexaketide target boundary: F9XMW3 1.484, A7UMW1 1.488, P9WET9 1.657. thn_synthase boundary: F9XMW3 1.546, A7UMW1 1.550, Q9Y7A7 (reference) 1.556.

## 4. Bootstrap support (MEASURED, 100 NNI replicates; `bootstrap_report/`)

FBP is a percent (0 to 100) and TBE a fraction (0 to 1), on the best ML tree. The values are for NNI replicates (see limits).

| Node | Leaves | FBP | TBE |
|---|---|---|---|
| thn_synthase target | 240 | 88 | 0.997 |
| hexaketide_pks target | 56 | 99 | 0.999 |
| heptaketide_pks target (per-step) = the node that bounds it | 22 | 100 | 1.000 |
| Alb1 + AlbA + wA group (smallest side with the three) | 13 | 100 | 1.000 |
| inside it: Alb1 + wA with 8 other leaves | 10 | 50 | 0.836 |
| the four references of the three steps (smallest side holding all 4) | 486 | 66 | 0.965 |

Clade-rule stability over the 100 bootstrap trees (`stability_per-step.tsv`, `stability_union.tsv`):

| Measure | Trees with it, of 100 |
|---|---|
| thn_synthase: a target exists | 100 |
| hexaketide_pks: a target exists | 100 |
| heptaketide_pks: a target exists, per-step rule | 100 |
| heptaketide_pks: a target exists, union rule | 2 |
| thn / hexa / hepta target leaf set identical to the best tree | 88 / 99 / 81 |
| Alb1 and AlbA form a clade that excludes wA (smallest side holding both does not hold wA) | 2 |

- The wA question, measured: Alb1 and AlbA are separated from wA in 2 of 100 NNI trees and in 0 of 2 standard trees. In the best tree the node of Alb1 + wA (with 8 other leaves) has FBP 50 and TBE 0.836 and sits inside the 13-leaf node of all three (FBP 100). So the tree does not resolve wA against Alb1 and AlbA. wA stays inside the smallest clade of the heptaketide references. Under the per-step rule this is not a boundary, so the target exists in 100 of 100 trees.
- This supports "wA is UNRESOLVED as an Alb1/AlbA outgroup"; it does not show that wA makes YWA1 in the tree sense or in any other sense.
- thn_synthase target support (FBP 88) is the lowest of the three target nodes.
- Limit: the NNI search. If the standard search resolved wA away from Alb1 and AlbA more often, the numbers above would change (UNVERIFIED; 2 standard trees agree with NNI on this point).

## 5. Negative and unexpected results

- The standard bootstrap did not finish in a useful time (4.1 and 4.7 h per replicate; 4 of 6 started tasks cancelled at 1 h 45 min). A `--fast` preset pilot was cancelled at about 58 min without a finished replicate.
- The first prepare job was cancelled to set CPUs apart; the tree job ended FAILED because the clades call in its spooled script had no `--step-decoys`; support job 32453906 read the wrong folder; rarefaction job 32455492 failed on a missing folder (all in `RUN.md`; all repaired, no result affected).
- The C4 check fails for 200 of 8,842 domains.
- Rarefaction shows no target change at all. This is expected from the rule (labelled leaves bound the targets) and is not independent evidence that the topology is stable.
- Support for the node joining the four references is low (FBP 66). The position of the three PKS references relative to each other is not well supported (UNVERIFIED reading: the deep nodes of a 559-leaf tree of PKS domains).

## 6. Outgroup: return to this

No outgroup is in the tree. No outside candidate passed the alignability rule (LovB Q0C8M3, EryA1 Q03131, yeast FAS1 P07149 and FAS2 P19097 failed; `inputs/screen_alignability.tsv` of the 0.7 run). The tree is unrooted for the clade rule; the display root is the decoys. A separate investigation will define how to choose an outgroup. **Return to this:** until an outgroup is chosen, no statement about the order of divergence of the three PKS clades is supported.

## 7. Not done here

- Standard (non-NNI) bootstrap; a test of the NNI bias beyond 2 replicates.
- A tree on a different named-leaf set; a rooted tree.
- Identity 0.9 or a set of representatives for the 200 domains outside the rule.
- A decision on whether to define the step by pathway role instead of by reaction (YAML). The step definition by reaction stays; the owner has not answered.
- Bootstrap of the 0.7 tree (not run).
- A hand check of the 157 + 118 unassigned proteins.

## 8. Decisions for the owner

1. Accept the label "sampling not yet adequate" (C4 fails for 2.3% of domains), or ask for identity 0.9.
2. Accept NNI-only bootstrap support, or pay about 70 to 80 h of CPU time for standard replicates (or give more CPUs).
3. Heptaketide headline is the per-step rule (target 22 leaves, 2,123 of 2,200 proteins). Confirm; the union result (no target) stays as the labelled comparison.
4. Define the PKS steps by pathway role instead of reaction (open).
