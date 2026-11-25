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
