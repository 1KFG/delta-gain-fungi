# Tree members for `phytoene_synthase` and `lycopene_cyclase` (carotenoid pathway)

Written 2026-11-22. Curator: Claude (agent). No person reviewed this list.

Labels: MEASURED = read from a fetched record or computed on 2026-11-22. PROPOSED = the choice of the curator. UNVERIFIED = not checked.

## What this folder is

This folder holds extra NAMED sequences for the step 5 reference tree. They are TREE MEMBERS ONLY.

- They are not references. They do not build the step HMM.
- They are not decoys. They do not set a score cutoff.
- They are not evidence for any step call.
- The absence check and the decoy-calibrated `min_score` must never read them.

Why they exist: the step 5 pilot for `phytoene_synthase` put the 5 references in one target clade of 186 leaves (MEASURED, `results/step5/carotenoid/phytoene_synthase/SUMMARY.tsv`). The tree had no other named sequence in that clade and no non-fungal outgroup. More named leaves let a later run test whether the clade splits into groups with a name.

## Limits

- **12 of the 13 named fungal members have no experiment in their UniProt record** (MEASURED). Protein existence is "3: Inferred from homology". The EC numbers, the function text and both Region features carry ECO:0000250 (by similarity to P37295, *Neurospora crassa* AL-2). The records cite genome papers only. The name of such a leaf is a UniProt curator transfer from P37295. It is not an independent observation of function.
- Only Q67GH9 (*Blakeslea trispora* carRA) has experimental evidence codes on its EC numbers (ECO:0000269, PMID 15345447). The paper was not read. Only the UniProt record was read.
- "Protein existence 1: Evidence at protein level" is the UniProt label. It does not say that an enzyme assay exists.
- The cut regions of the 12 by-similarity members are UniProt Region features with ECO:0000250. They are positions transferred by similarity, not measured domain borders. For Q67GH9 the Region features carry ECO:0000305 (curator inference).
- No tree was built with these files. The statement that they can split the target clade is a hypothesis.
- `bin/step5_place.py` and `slurm/step5_prepare.sbatch` have no input for named members today (read on 2026-11-22). `PC_OUTGROUP_FA` takes one FASTA file and gives every sequence in it the role `outgroup` and the rooting role. Do not pass the named-member files through `PC_OUTGROUP_FA`. `outgroups_phytoene_synthase.faa` fits `PC_OUTGROUP_FA`. The named members need a new option in the tool. This is an owner decision.
- The choice of the 7 outgroups is PROPOSED.
- `verified` in `provenance.tsv` has a narrow meaning here: the record was fetched again, and the organism, the protein name (and the gene name for outgroups) and the EC number matched the fetched record. It does not mean that the function is experimentally shown. Read `protein_existence` and `evidence_codes` for that.

## How the list was made (MEASURED, 2026-11-22)

1. UniProt REST search `(ec:2.5.1.32) AND (taxonomy_id:4751) AND (reviewed:true)`: 19 entries. 5 are the references of the step (P37295, Q7Z859, Q8X0Z1, Q9P854, Q9UUQ6). 14 remain.
2. UniProt REST search `(ec:1.3.99.31) AND (taxonomy_id:4751) AND (reviewed:true)`: 0 entries. This EC (phytoene desaturase) adds nothing to this list.
3. Each of the 14 accessions was fetched again as a full UniProt JSON record. `build_tree_members.py` checks the accession, the organism, the protein name and both EC numbers (2.5.1.32 and 5.5.1.19).
4. Cut rule: the same rule as for the references (`curation_scripts/carotenoid/build.py`). The synthase cut is the UniProt feature Region "Phytoene synthase". The cyclase cut is the UniProt feature Region "Lycopene beta-cyclase". The script also runs the same containment check as `build.py`: every InterPro match of Pfam PF00494 lies inside the synthase region, and every InterPro match of NCBIfam TIGR03462 lies inside the cyclase region. All 13 members pass (MEASURED).
5. 13 of 14 entries have both Region features. 1 entry does not (see "Excluded").

## Named fungal members (13)

All values are MEASURED from the fetched UniProt records. "Function PMID" is a PubMed id that the record cites with the position FUNCTION or CATALYTIC ACTIVITY. The two HMM columns are the best domain score of the cut domain against `results/hmms/carotenoid/<step>.PROVISIONAL.hmm` (HMMER 3.4, `hmmsearch -E 10`).

| Accession | Organism | Class | Phylum | Gene | PE | Score | EC evidence | Function PMID | Synthase region | Cyclase region | Synthase HMM | Cyclase HMM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Q67GH9 | Blakeslea trispora | Mucoromycetes | Mucoromycota | carRA | 1 | 5.0 | ECO:0000269, experimental | 15345447 | 247-608 | 1-240 | 431.5 | 311.2 |
| B2ATB0 | Podospora anserina | Sordariomycetes | Ascomycota | Pa_1_15240 | 3 | 4.0 | ECO:0000250, by similarity | none | 249-593 | 1-242 | 410.6 | 304.8 |
| D1Z4K7 | Sordaria macrospora | Sordariomycetes | Ascomycota | SMAC_01277 | 3 | 4.0 | ECO:0000250, by similarity | none | 248-602 | 1-241 | 476.1 | 348.5 |
| E3Q717 | Colletotrichum graminicola | Sordariomycetes | Ascomycota | GLRG_02475 | 3 | 4.0 | ECO:0000250, by similarity | none | 249-587 | 1-242 | 350.3 | 284.6 |
| B2WAQ3 | Pyrenophora tritici-repentis | Dothideomycetes | Ascomycota | PTRG_07366 | 3 | 4.0 | ECO:0000250, by similarity | none | 250-583 | 1-243 | 359.3 | 300.5 |
| E4ZUB5 | Leptosphaeria maculans | Dothideomycetes | Ascomycota | Lema_P114090.1 | 3 | 4.0 | ECO:0000250, by similarity | none | 250-581 | 1-243 | 352.3 | 303.2 |
| Q0V6M5 | Parastagonospora nodorum | Dothideomycetes | Ascomycota | SNOG_00339 | 3 | 4.0 | ECO:0000250, by similarity | none | 250-585 | 1-243 | 373.1 | 287.8 |
| D4D802 | Trichophyton verrucosum | Eurotiomycetes | Ascomycota | TRV_03236 | 3 | 4.0 | ECO:0000250, by similarity | none | 256-595 | 1-249 | 301.7 | 251.6 |
| E4UPP6 | Nannizzia gypsea | Eurotiomycetes | Ascomycota | MGYG_02096 | 3 | 4.0 | ECO:0000250, by similarity | none | 256-594 | 1-249 | 315.2 | 249.0 |
| Q2U4X9 | Aspergillus oryzae | Eurotiomycetes | Ascomycota | AO090020000159 | 3 | 4.0 | ECO:0000250, by similarity | none | 243-587 | 1-236 | 235.3 | 237.3 |
| A2QM49 | Aspergillus niger | Eurotiomycetes | Ascomycota | An07g00800 | 3 | 3.0 | ECO:0000250, by similarity | none | 268-582 | 1-261 | 131.2 | 237.6 |
| Q4P0S6 | Mycosarcoma maydis | Ustilaginomycetes | Basidiomycota | UMAG_06287 | 3 | 4.0 | ECO:0000250, by similarity | none | 262-686 | 1-255 | 137.7 | 208.2 |
| E6ZZ11 | Sporisorium reilianum | Ustilaginomycetes | Basidiomycota | sr13176 | 3 | 3.0 | ECO:0000250, by similarity | none | 262-685 | 15-255 | 133.6 | 209.4 |

Full strain names are in `provenance.tsv`. Gene: the gene name when the record has one (Q67GH9 only), else the ordered locus name, else the ORF name.

Things to know before use (all MEASURED):

- The by-similarity source is P37295 in all 12 cases. P37295 is a reference of both steps.
- The synthase HMM hits Q4P0S6 and E6ZZ11 in two pieces (HMM positions 35-208 and 209-344; second piece 85.1 and 91.8). The table gives the best piece. Their synthase regions are 425 and 424 residues, longer than the HMM (350 match states).
- The synthase HMM hits A2QM49 in one piece that covers HMM positions 25-223 only. Its synthase region is 315 residues, the shortest of the 13.
- 12 of the 13 synthase scores are below the PROVISIONAL `proposed_min_score` of the step (457.9, `results/step2/carotenoid/step_scores.PROVISIONAL.tsv`). Only D1Z4K7 (476.1) is above it. All 13 cyclase scores are at or below the cyclase value in `LYCOPENE_CYCLASE_DECOYS.md` (347.5) except D1Z4K7 (348.5). These members did not set those values and must not. The numbers are given so the owner can see them. What they mean for the cutoff is not decided here.
- Class coverage: Sordariomycetes 3, Dothideomycetes 3, Eurotiomycetes 4, Ustilaginomycetes 2, Mucoromycetes 1.

## Excluded

| Accession | Organism | Gene | Why |
|---|---|---|---|
| S0EQ07 | Gibberella fujikuroi (strain CBS 195.34 / IMI 58289 / NRRL A-6831) | carRA | The fetched record has no Region feature (MEASURED). The cut rule gives no cut, and no cut was guessed. |

More on S0EQ07 (MEASURED): PE 1, annotation score 4.0, EC evidence ECO:0000269 with PMID 12172798, length 580. PMID 12172798 is also the paper of the reference Q8X0Z1 (*Fusarium fujikuroi* carRA, length 612). A plain string comparison (Python `difflib`, not an alignment) finds 578 of the 580 residues of S0EQ07 in matching blocks of Q8X0Z1; the longest identical block is 532 residues. The full-length protein scores 522.3 on the synthase HMM (alignment at residues 243-580). It is very probably the same gene as the reference from another strain, so it would add almost no new information to the tree. That last sentence is an inference from the numbers, not a statement in a record.

## Tree outgroups for `phytoene_synthase` (7)

File: `outgroups_phytoene_synthase.faa`. The sequences are FULL LENGTH. The records have no Region feature for the enzyme domain, so the cut rule does not apply. The two plant proteins include a transit peptide (UniProt feature, ECO:0000255, predicted: residues 1-70 in P37271, 1-74 in P37272).

All record values are MEASURED. HMM score: best domain score of `hmmsearch -E 10` with `results/hmms/carotenoid/phytoene_synthase.PROVISIONAL.hmm` (350 match states, 5 references). All 7 give a hit.

| Accession | Organism | Class, phylum | Gene | Type | EC | PE | Score | EC evidence | Function PMID | Length | Best domain score | Aligned residues |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P21683 | Pantoea ananas | Gammaproteobacteria, Pseudomonadota | crtB | bacterial CrtB | 2.5.1.32 | 1 | 4.0 | ECO:0000269 | 9593819 | 309 | 154.7 | 7-284 |
| P54975 | Paracoccus sp. (strain N81106 / MBIC 01143) | Alphaproteobacteria, Pseudomonadota | crtB | bacterial CrtB | 2.5.1.32 | 1 | 3.0 | ECO:0000269 | 7592436 | 304 | 102.7 | 88-271 (second piece 29.3 at 8-52) |
| P37269 | Synechococcus elongatus (strain PCC 7942) | Cyanophyceae, Cyanobacteriota | crtB | cyanobacterial CrtB | 2.5.1.32 | 1 | 3.0 | ECO:0000269 | 1537409 | 308 | 107.1 | 25-296 |
| P37294 | Synechocystis sp. (strain PCC 6803) | Cyanophyceae, Cyanobacteriota | crtB | cyanobacterial CrtB | 2.5.1.32 | 1 | 3.0 | ECO:0000269 | 8018713 | 337 | 89.0 | 43-324 |
| P37271 | Arabidopsis thaliana | Magnoliopsida, Streptophyta | PSY1 | plant PSY | 2.5.1.32 | 1 | 5.0 | ECO:0000269 | 25675505 | 422 | 102.7 | 141-403 |
| P37272 | Capsicum annuum | Magnoliopsida, Streptophyta | PSY1 | plant PSY | 2.5.1.32 | 1 | 5.0 | ECO:0000269 | 16578835 (also 8250898, 20582146) | 419 | 104.1 | 140-402 |
| Q2FV59 | Staphylococcus aureus (strain NCTC 8325 / PS 47) | Bacilli, Bacillota | crtM | bacterial CrtM (4,4'-diapophytoene synthase) | 2.5.1.96 | 1 | 4.0 | ECO:0000269 | 12426357, 16269684 (also 22535955) | 287 | 59.2 | 10-276 |

For comparison (MEASURED, same search): the 4 squalene synthase decoys of the step score 49.3, 43.4, 37.1 and 30.8. The 5 references score 457.9 to 537.9. All 7 outgroups score above every decoy. This is one more reason why they must not be read as decoys: `max_decoy_score` would rise from 49.3 to 154.7.

Outgroup candidates that were fetched and not added (MEASURED):

| Accession | Organism | Gene | Best domain score | Why not added |
|---|---|---|---|---|
| D5KXJ0 | Enterobacter agglomerans (Pantoea agglomerans) | crtB | 149.4 | Second Pantoea CrtB. P21683 covers the genus. |
| P08196 | Solanum lycopersicum | PSY1 | 99.8 | The EC number is by similarity to P37272 (ECO:0000250). P37272 itself was taken. |
| A9JQL9 | Staphylococcus aureus | crtM | 59.9 | Second CrtM of the same species. This record has PDB cross-references. |
| O07854 | Staphylococcus aureus (strain Newman) | crtM | 59.2 | Third CrtM of the same species. |

No outgroup was added for `lycopene_cyclase`. The existing outgroup is Q9HNE5 in `refs/carotenoid/outgroups/`. No search for a better one was made.

## Who reads this folder (checked by grep on 2026-11-22)

`grep -rn "tree_members"` over `*.py`, `*.sbatch`, `*.sh`, `*.nf`, `*.yaml`, `*.toml` under `analysis/pathway_census/`, outside this folder: 0 lines before `tests/test_tree_members.py` was added.

- `bin/score_refs_decoys.py` lists the steps with `refs_dir.glob("*.faa")`. The glob is not recursive. It then reads `refs_dir/<step>.faa` and `refs_dir/decoys/<step>.faa`. With `--refs-dir refs/carotenoid` (as in `slurm/score_refs_decoys.sbatch`) it does not read this folder.
- `bin/absence_check.py` reads `refs_dir/<step>.faa` only.
- `bin/build_step_hmms.py` reads `refs_dir/<step>.faa` only.
- `bin/make_min_scores.py` reads the output table of `score_refs_decoys.py`, no FASTA file.
- `pathways/carotenoid.yaml` has no accession of this folder.
- `curation_scripts/carotenoid/build.py` creates its output directories with `exist_ok` and writes named files. It has no delete call. It does not read or write this folder.

Risk: `phytoene_synthase.faa` and `lycopene_cyclase.faa` in this folder have step names. A run with `--refs-dir refs/carotenoid/tree_members` would read them as references. Never pass this folder as `--refs-dir`. `tests/test_tree_members.py` checks that the scoring and absence scripts do not name the folder, that no member is a reference, decoy or outgroup of the same step, and that no member is in the YAML.

## Files

| File | Content |
|---|---|
| `phytoene_synthase.faa` | 13 cut synthase domains. Header: accession, `region=START-END`, then the UniProt header. |
| `lycopene_cyclase.faa` | 13 cut cyclase domains of the same 13 proteins. |
| `full_length/named_members.faa` | The 13 full-length sequences. |
| `outgroups_phytoene_synthase.faa` | 7 full-length outgroup sequences. |
| `provenance.tsv` | 33 rows: 26 `named_member` (13 per step) and 7 `tree_outgroup`. `pmid` is empty when the record cites no paper for the function. |
| `build_tree_members.py` | Writes all files above from the fetched records. A second run gave the same bytes (md5, MEASURED). Fetches use the cache in `analysis/pathway_census/.cache/`. InterPro returned HTTP 408 once in the first run; the second run passed. |
| `README.md` | This file. |

Not done here: no tree, no alignment, no change to any reference, decoy, YAML, HMM or cutoff, and no full-text paper was read.
