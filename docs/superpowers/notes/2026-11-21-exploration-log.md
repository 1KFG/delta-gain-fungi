# Exploration log (pathway knowledge base and the EBC metabolomics side)

Kept by the BFD-HPCC-master agent. Last updated 2026-11-21 11:30 local. Newest status first within each section. Labels: **MEASURED** = I ran or fetched it; **REPORTED** = from an EBC agent message (named); **PROPOSED**; **UNVERIFIED**. Nothing here is a biological finding.

## 1. Where each exploration stands

| Exploration | Status | One-line result | Source |
|---|---|---|---|
| Carotenoids (genome side) | Stage 1 curation in progress | Step set v2 curated: 29 verified references, 61 decoy rows, 36 evidence links; under independent review | plan `2026-11-21-stage1-verified-step-sets.md`, Task 3 report |
| Carotenoids (corpus comparison) | Closed as "not testable" | 0 of 84 fungal-reference accessions show beta-carotene, zeaxanthin or canthaxanthin; astaxanthin 1, retinal 2, torulene 0 | REPORTED, `ebc-006` |
| DHN-melanin, L-DOPA melanin (genome side) | Not started (Tasks 4 and 5) | Step sets corrected after review (see spec 6.3) | spec v2 |
| Melanin (corpus comparison) | Closed as "not testable" | Five of seven library names have zero hits; L-DOPA 0 fungal-reference groups; formulas are background (10-35% of groups) | REPORTED, `ebc-006` |
| Fucose | Dropped (weak target) | Free fucose and rhamnose: 17 library features in 13 accessions; none in fungal-reference or milk accessions; positive controls mostly fail | REPORTED, `ebc-008`, `ebc-012` |
| Mycosporines and MAAs (genome side) | Added; curation is Task 6, not started | Leads only: PMID 36983003, 40036870; cyanobacterial DDGS entries Q3M6C3, B2J6X9 | spec addendum 2.1 |
| MAAs (corpus) | Done by EBC | A real marine positive control; no fungal signal (below) | REPORTED, `ebc-010`, `ebc-013` |
| Amphibians | Done by EBC (census only) | No skin, mucus, gland or swab samples in the corpus; compound panel not run | REPORTED, `ebc-011`, `ebc-014` |
| SIRIUS on the 734 EBC features | Running (4 serial container jobs) | Container pilot clean; chunk 1 running | MEASURED |
| Methods check on a detectable pathway (plan B) | Open | Fucose fails; MAA has no fungal signal; fungal natural products (roquefortine, andrastin) still untested | spec 6.6, request 006 |

## 2. Mycosporine-like amino acids in the corpus (REPORTED, `ebc-010`, `ebc-013`)
- 36 of 87 marine accessions carry an MAA library name, against 216 of 5,379 others (odds ratio 16.9, permutation p 0.0005). Unchanged after removing the top three contributing accessions, the largest marine study, all fungal-labelled accessions, and after restricting to precursors within 5 mDa (odds ratio 15.7).
- Eight MAA species have library names: palythine (127 accessions), palythene (98), usujirene (87), porphyra-334 (76), mycosporine-methylamine-threonine (18), mycosporine-glycine-alanine (12), mycosporine glutaminol (9), mycosporine serinol (7). Shinorine, asterina-330, palythinol, mycosporine-glycine and mycosporine-2-glycine have **zero** library hits.
- Fungi: 0 of 82 fungal-reference accessions. Two MIXED_FUNGAL groups carry names (dominant taxon *Pleurotus ostreatus*; and *Penicillium*), but those groups hold non-fungal runs, so the hits cannot be assigned to the fungus. One group with no cohort label (dominant taxon *Agaricus*) carries three names; its runs were not checked.
- Isobars: palythene and usujirene share m/z 285.14; mycosporine-methylamine-threonine and glutaminol share 303.154. Formula-only hits are background in every cohort (MSI 4).
- **Gate reading:** fails for fungi. Marine enrichment shows presence, not producer (cyanobacteria, algae, fungi and animals that take MAAs up all occur; UNVERIFIED).

## 3. Amphibians (REPORTED, `ebc-011`, `ebc-014`)
- 1,778 amphibian-labelled runs in 24 accessions: 1,586 runs in 23 accessions are host samples (mostly laboratory *Xenopus* eye, brain and testis; also bile, serum, feces, liver, axolotl forelimb, *Ambystoma maculatum* eggs with the alga *Oophila*); 99 runs in one accession (MSV000090464) are *Batrachochytrium dendrobatidis* and *B. salamandrivorans* cultures.
- Skin, mucus, gland or swab labels: 0 runs. A skin-community screen is not possible in this corpus. Usable feature tables: 12 of 24 accessions (17 groups, 900 runs). One study (ST002414, *Xenopus* eye) holds 25% of runs.
- Narrower options not started (need the operator's approval on the EBC side): the Bd/Bsal culture accession; the *Ambystoma*–*Oophila* egg accessions. Tetrodotoxin hits (14 non-marine accessions) not checked for an amphibian host.

## 4. BFD genome coverage for the taxa named by the EBC notes (MEASURED 2026-11-21, `species.parquet`)
| Genus | Genomes | Species |
|---|---|---|
| Batrachochytrium | 15 | 2 (*B. dendrobatidis* and *B. salamandrivorans*; per-species counts not tallied) |
| Pleurotus | 53 | 17 (*P. ostreatus*: 10 genomes) |
| Agaricus | 20 | 6 |
| Penicillium | 655 | 192 |
| Aspergillus | 1,632 | 273 |
| Naganishia | 28 | 17 |
| Rhodotorula | 220 | 26 |
| Phaffia | 12 | 3 |
| Fusarium / Neurospora | 2,174 / 194 | 291 / 35 |
| Saccharomyces / Komagataella | 1,566 / 188 | 22 / 7 |
These counts say a genome exists. They say nothing about pathway genes. Whether any of these genomes carry mycosporine genes is the question for Task 6 and the later census.

## 5. Answers owed to, or from, the EBC side
- Owed to EBC: a fungal MAA gene list and the taxa that carry it (waits on Task 6 and the census); whether BFD covers *Batrachochytrium* biosynthetic genes (genomes exist, 15; no pathway curated). BFD has fungal genomes only, so bacterial skin-microbe pathways are not covered.
- Owed from EBC: nothing blocking. Their SIRIUS template was changed to fail when no summary tables appear and not to pre-create the output path (commit 38f9a15, REPORTED).

## 6. Decisions on record
- 2026-11-21 user: reduced plan accepted (flat tables, M1 dropped, thresholds in spec 6.6); melanin included; mycosporines added; SIRIUS full run approved; message `006` pushed on the user's go.
- Gate (spec 6.6) applied as written: carotenoids and melanin are reported as "not testable in this corpus"; thresholds not lowered.

## 7. Open questions
1. Which detectable pathway, if any, serves as the methods check (fungal natural products, MAAs outside fungi, or none)?
2. Whether to release the tables outside the group (licence rule).
3. Which fungal taxa serve as negative controls (candidates only; confirm from the census).

## 8. Jobs and files (for restart)
- SIRIUS chain: jobs 32355517 to 32355520 (serial, `afterany`); run folder `/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/analysis/ms_comparison/sirius_full/` (project output in `outputs_v2/`); container image `/bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif`. Two empty folders from a failed first attempt (`outputs/orbitrap_pos_a`, `outputs/orbitrap_pos_b`) await manual removal.
- Collab repo: `stajichlab/EBC-KT`, folder `collab/`; local clone in the session scratchpad; guarded sync script `sync_collab.sh` there. Sign as BFD-HPCC-master agent.
- Branch `pathway-census-realdata` (worktree `DeltaGain_Fungi_wt_realdata`), not pushed.
