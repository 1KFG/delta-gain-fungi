# Exploration log (pathway knowledge base and the EBC metabolomics side)

Kept by the BFD-HPCC-master agent. Last updated 2026-11-21 (evening; Tasks 1-5 closed, Task 6 running, branch pushed). Newest status first within each section. Labels: **MEASURED** = I ran or fetched it; **REPORTED** = from an EBC agent message (named); **PROPOSED**; **UNVERIFIED**. Nothing here is a biological finding.

## 1. Where each exploration stands

| Exploration | Status | One-line result | Source |
|---|---|---|---|
| Carotenoids (genome side) | Stage 1 curation in progress | Step set v2 curated: 29 verified references, 61 decoy rows, 36 evidence links; under independent review | plan `2026-11-21-stage1-verified-step-sets.md`, Task 3 report |
| Carotenoids (corpus comparison) | Closed as "not testable" | 0 of 84 fungal-reference accessions show beta-carotene, zeaxanthin or canthaxanthin; astaxanthin 1, retinal 2, torulene 0 | REPORTED, `ebc-006` |
| DHN-melanin (genome side) | Task 4 closed | Review found the two DHN reductases swapped in the pathway review (Q9HFV6 = 4HNR, Q12634 = 3HNR); 6 user decisions in `refs/melanin_dhn/UNVERIFIED.md` | Task 4 report, ledger |
| L-DOPA and other melanins (genome side) | Task 5 closed | 15 refs, 29 links: 2 T1a, 3 T1b, 2 T1c, 20 T3, 2 T4 after review applied the quote-names-the-protein rule; melB reference is A0A1S9DK56 (PDB 3W6W), not Q2UP46; tyrosinase_dopa T1c in Neurospora crassa only | Task 5 report, review, re-review |
| Melanin (corpus comparison) | Closed as "not testable" | Five of seven library names have zero hits; L-DOPA 0 fungal-reference groups; formulas are background (10-35% of groups) | REPORTED, `ebc-006` |
| Fucose | Dropped (weak target) | Free fucose and rhamnose: 17 library features in 13 accessions; none in fungal-reference or milk accessions; positive controls mostly fail | REPORTED, `ebc-008`, `ebc-012` |
| Mycosporines and MAAs (genome side) | Task 6 running (Fable implementer) | Leads only: PMID 36983003, 40036870; cyanobacterial DDGS entries Q3M6C3, B2J6X9 | spec addendum 2.1 |
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

## 9. Decisions and lessons recorded 2026-11-21 (evening)
- Branch `pathway-census-realdata` pushed to origin on the user's instruction (Stage 1 through Task 5, commit e1e82da). Not merged; merging is the owner's call.
- Rule applied from review of Task 5: a quote that does not name the link's own protein and a substrate or product gives T3 (T1c for a colour phenotype naming the protein). A "substrate class" reading was tried and withdrawn.
- Review lesson: a UniProt accession can be the right locus but not the assayed sequence (Q2UP46 vs A0A1S9DK56). Check sequence against the deposited structure.
- Review lesson: reviewer-recalled accessions and PMIDs are leads only; each is re-fetched.
- Seven Task 5 decisions for the user are in the ledger and `refs/melanin_ldopa/UNVERIFIED.md` (Q2TL94, melB ORF-name reference, P00440 tie, LAC1 T1c without substrate, TyrP T3, no asp_melanin product, ORF-name-only rule).
- Next: Task 6 review, Task 7 verifier, Task 8 packet, then the user review gate before any real-data run.
- The shareable synthesis lives at `analysis/pathway_census/REPORT.md` (page `reports/pathway_corpus.qmd`, linked from the reports index). Update it when Task 6 closes or the gate counts change.

## 10. Earlier MAA/HGT work by the owner, found 2025 folder `/bigdata/stajichlab/shared/projects/Mycosporine` (read 2026-11-21, MEASURED by listing and reading the Rmd; no results re-run)
- Queries (`query.pep`): ATPG B2VZE5 (*Pyrenophora tritici-repentis*, ATP-grasp, PTRG_02785), DDGS B2VZE7 (same species, UniProt name "2-epi-5-epi-valiolone synthase", PTRG_02787), OMT A1C9L2 (*Aspergillus clavatus*, "catechol O-methyltransferase"). phmmer (E 1e-10) against the 1KFG 5k fungal proteome set (56.7 M sequence lines in `db/fungi_5k.fasta`), then per-query alignments, HMMs and FastTree trees.
- Hit counts in `results/parallel_phmmer_MSA/ATPG_DDGS_OMT.Rmd` comments: species with a hit for ATPG 2,347, DDGS 5,633, OMT 3,184; all three 2,118; ATPG with DDGS 2,341. Hit files: 2,604 ATPG, 8,123 DDGS, 4,999 OMT proteins.
- Not found: any written HGT conclusion, any source/target call, any notes. The folder's git repo has no commits. The Rmd makes presence/absence and Venn/phylum summaries only.
- Reading (PROPOSED): the three queries are family-level, not step-specific. The DDGS query is an EEVS-type protein and the OMT query a catechol OMT; both families are broad (Task 6 lists the same families as decoys and found the Phaffia OMT closer to catechol OMTs than to MysB). So "2,118 species with all three" is a family census, not a mycosporine capability call. ATPG (2,347 species) is the rarest and the most informative axis; ATPG plus DDGS adjacent in the genome would be a stronger criterion (the folder has a GFF-pull step for the all-three species, but its output holds almost nothing: 3 lines).
- Useful for the census: the Phaffia reference proteins (WEV88679-81) can be searched against the same db, and the existing trees could be used to place Phaffia and the Pyrenophora proteins relative to the decoys. Not done.
