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

---

# Addition: tree members for `phytoene_desaturase` and `cart`

Written 2026-11-22. Curator: Claude (agent). No person reviewed this list. The labels MEASURED, PROPOSED and UNVERIFIED have the same meaning as above. Everything in "What this folder is" applies to these files too: TREE MEMBERS ONLY, not references, not decoys, not evidence for a step call.

Why: step 5 for these two steps has no leaf that can root the tree. The decoys of both steps cover less than 20% of the step HMM with their best domain (MEASURED, `bin/filter_decoys_for_tree.py`), so they cannot be tree leaves.

## Limits of this addition

- **`phytoene_desaturase` has 1 named fungal member, and it has no experiment in its record. `cart` has 0 named fungal members** (MEASURED). UniProt has no other reviewed fungal entry for these two enzymes. These files cannot split the target clade into lineage groups. One labelled leaf per class is below the "at least 2 labelled leaves" rule of `bin/step5_place.py`.
- **All 5 `cart` outgroups are vertebrate proteins (BCO1 and RPE65).** No plant, cyanobacterial, bacterial or other fungal member of the carotenoid oxygenase family passes the alignability rule (MEASURED, table below). The outgroup of `cart` is one lineage. It can root the tree. It cannot show where plant or bacterial cleavage oxygenases attach.
- **No plant or cyanobacterial PDS or ZDS passes for `phytoene_desaturase`.** Their best domain covers 9% of the HMM (MEASURED). The outgroups are bacterial CrtI-type enzymes and two distant proteins (mouse RetSat, tomato CRTISO).
- Q8S4R4 (tomato CRTISO) passes with coverage 0.501. The limit is 0.5. The *Arabidopsis* CRTISO Q9M9Y8 fails with 0.496. The two best domains cover 290 and 287 of the 579 HMM positions (MEASURED); 290 is the smallest count that passes.
- The alignability rule reads the best single HMM domain only. In the carotenoid oxygenase family the HMM hit breaks into 2 or 3 domains. Example (MEASURED): human BCO2 Q9BYV7 fails with coverage 0.187 for its best domain, and its three domains of 10 bits or more together cover 0.691 of the HMM. The rule was applied as written. A change of the rule is an owner decision.
- A mafft alignment with the references does not separate passing from failing candidates (MEASURED, see "mafft check"). It is a global alignment and places residues of every sequence. It is not evidence of homology over the full length.
- No tree was built with these files. The choice of outgroups is PROPOSED.
- Only UniProt records were read. No paper was read. `verified` has the same narrow meaning as above.
- The numbers depend on the PROVISIONAL HMMs of 2026-11-22: `phytoene_desaturase.PROVISIONAL.hmm` (579 match states, 4 references, md5 a6d90200a813da2121bac04c082ad653) and `cart.PROVISIONAL.hmm` (621 match states, 2 references, md5 c9306e38ff00ffc46192fb17974b0420). A new HMM needs a new measurement.

## How the lists were made (MEASURED, 2026-11-22)

1. `pathways/carotenoid.yaml` gives `ec: null` for `phytoene_desaturase` and `ec: 1.13.11.59` for `cart`.
2. UniProt REST searches, all with `taxonomy_id:4751` and `reviewed:true`:
   - `ec:1.3.99.31`: 0 entries. `ec:1.3.99.30`: 2 entries (P48537, P54982).
   - `protein_name:"phytoene desaturase" OR gene:carB OR gene:al-1 OR gene:crtI`: 5 entries (P21334, Q8X0Z0, S0EPU6, P54982, P48537).
   - `family:"carotenoid/retinoid oxidoreductase family"`: the same 5 entries.
   - `ec:1.13.11.59` (no taxon filter): 3 entries (A1KQY3, A7UXI1, S0E9A6).
   - `family:"carotenoid oxygenase family"`: 11 entries.
3. Outgroup candidates: the reviewed non-fungal entries of the two UniProt families (71 and 58 entries), plus plant and cyanobacterial ZDS and PDS entries from a gene-name search.
4. Alignability of every candidate: `hmmsearch -E 1000 --domE 1000` (HMMER 3.4) with the step HMM, run through `bin/filter_decoys_for_tree.py --min-score 10 --min-coverage 0.5`. Coverage = (hmm_to - hmm_from + 1) / HMM length for the best domain.
5. Each accession that was written was fetched again as a full UniProt JSON record. `build_tree_members.py` checks the organism, the gene name, the protein name, the EC number and the family comment. For outgroups it also checks that the EC number has experimental evidence (ECO:0000269).

## Named fungal members

### `phytoene_desaturase` (1)

| Accession | Organism | Class | Phylum | Gene | PE | Score | Evidence | Function PMID | Length | Best domain score | HMM coverage |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P48537 | Cercospora nicotianae | Dothideomycetes | Ascomycota | PDH1 | 3 | 3.0 | ECO:0000250, by similarity | none | 621 | 776.2 | 0.801 |

All values are MEASURED. Things to know before use:

- The record is 'by similarity' only. The FUNCTION text carries ECO:0000250 and names no source protein. The EC number (1.3.99.30) and the catalytic activity line carry no evidence code. The record cites no paper for the function.
- The sequence is FULL LENGTH (1-621). Nothing was cut. The header has `region=1-621` to keep the format of the other member files. The record has a signal peptide feature at 1-23 (ECO:0000255, sequence model); it was not removed.
- The HMM hit has a second domain (46.5 bits, HMM positions 445-515). The record has a "Disordered" region at residues 394-425.
- The 4 references are 2 Sordariomycetes (P21334, Q8X0Z0), 1 Mucoromycetes (P54982) and 1 Tremellomycetes (O13506), from the lineage of the fetched records. P48537 is the only Dothideomycetes leaf.
- P48537 is appended to `full_length/named_members.faa` (record 14).

### `cart` (0)

There is no file `cart.faa` in this folder. No reviewed fungal entry with EC 1.13.11.59 exists outside the two references and one strain copy (MEASURED). An empty file with a step name was not written.

### Fetched and excluded: the same gene as a reference

| Accession | Organism | Gene | Step | Reference | Difference (MEASURED) | Record evidence |
|---|---|---|---|---|---|---|
| S0EPU6 | Gibberella fujikuroi (strain CBS 195.34 / IMI 58289 / NRRL A-6831) | carB | phytoene_desaturase | Q8X0Z0 (Fusarium fujikuroi carB) | 0 of 570 residues: the sequences are identical | PE 1, score 4.0, ECO:0000269, PMID 10905351, 12172798, 19645721 |
| S0E9A6 | Gibberella fujikuroi (same strain) | carT | cart | A1KQY3 (Fusarium fujikuroi carT) | 7 of 577 residues differ (positions 18, 64, 155, 179, 183, 187, 269), same length | PE 1, score 3.0, ECO:0000269, PMID 17493127 |

Both are the same gene of the same species as a reference, from another strain record. They would add no information to the tree. `build_tree_members.py` checks both difference counts.

The other reviewed fungal entries of the carotenoid oxygenase family are decoys of `cart` (Q7S860, A0A0D1E6L2, A0A384JIP2, Q4W9H1, G2QIL8), a reference of the step `carx` (S0ENQ2), or have another function (Q2HHS4 RCO1, W3VHW6 aromatic dioxygenase). None of them is a `cart` member, and none passes the alignability rule (table below).

## Tree outgroups for `phytoene_desaturase` (7)

File: `outgroups_phytoene_desaturase.faa`. FULL-LENGTH sequences. Q8S4R4 has a chloroplast transit peptide feature without a known end position in the record; nothing was removed. Q64FW2 has a signal peptide feature at 1-25; nothing was removed.

All record values are MEASURED. Every one has "Belongs to the carotenoid/retinoid oxidoreductase family" in the UniProt record and ECO:0000269 on its EC number or catalytic activity line.

| Accession | Organism | Class, phylum | Gene | Type | EC | PE | Score | Function PMID | Length | Best domain score | HMM coverage |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P21685 | Pantoea ananas | Gammaproteobacteria, Pseudomonadota | crtI | bacterial CrtI, lycopene-forming | 1.3.99.31 | 1 | 4.0 | 1400305 | 492 | 288.7 | 0.853 |
| P17054 | Rhodobacter capsulatus (strain SB1003) | Alphaproteobacteria, Pseudomonadota | crtI | bacterial CrtI, neurosporene-forming | 1.3.99.28 | 1 | 5.0 | 8830054 | 524 | 317.0 | 0.874 |
| P54979 | Myxococcus xanthus | Myxococcia, Myxococcota | carA2 | bacterial phytoene desaturase, zeta-carotene-forming | 1.3.99.29 | 1 | 3.0 | 7588751, 17662111 | 517 | 494.4 | 0.862 |
| Q02861 | Myxococcus xanthus | Myxococcia, Myxococcota | carC | bacterial all-trans-zeta-carotene desaturase | 1.3.99.26 | 1 | 3.0 | 17662111 | 529 | 302.7 | 0.865 |
| Q7WT72 | Flavobacterium sp. (strain P99-3) | Flavobacteriia, Bacteroidota | crtD | bacterial 1-hydroxycarotenoid 3,4-desaturase (other function) | 1.3.99.27 | 1 | 4.0 | 15251462 | 488 | 323.1 | 0.851 |
| Q64FW2 | Mus musculus | Mammalia, Chordata | Retsat | animal all-trans-retinol 13,14-reductase (other function, CrtISO subfamily) | 1.3.99.23 | 1 | 5.0 | 15358783, 17253779, 19139408 (also 19940255, 28855500) | 609 | 45.5 | 0.544 |
| Q8S4R4 | Solanum lycopersicum | Magnoliopsida, Streptophyta | CRTISO | plant prolycopene isomerase (other function, CrtISO subfamily) | 5.2.1.13 | 1 | 4.0 | 11884678 (EC evidence also 15557094) | 615 | 41.4 | 0.501 |

For comparison (MEASURED, same search): the 4 references score 986.1 to 1048.9. P54979 scores 494.4, the highest non-fungal score. The 5 bacterial outgroups score far above every decoy (highest decoy: 24.3). They must not be read as decoys.

Two groups are in this file on purpose (PROPOSED): 5 bacterial CrtI-type enzymes (score 288 to 495) and 2 distant CrtISO-subfamily proteins (score 41 to 46). Which group roots the tree is an owner decision. `PC_OUTGROUP_FA` roots on the MRCA of all sequences in the file.

### Candidates that FAIL the alignability rule (not written)

| Accession | Entry | Organism | Best domain score | HMM coverage | Fails on |
|---|---|---|---|---|---|
| Q07356 | PDS_ARATH | Arabidopsis thaliana | 14.3 | 0.093 | coverage |
| P28554 | PDS_SOLLC | Solanum lycopersicum | 17.3 | 0.092 | coverage |
| P80093 | PDS_CAPAN | Capsicum annuum | 17.5 | 0.095 | coverage |
| P26294 | PDS_SYNE7 (pds, crtP) | Synechococcus elongatus PCC 7942 | 10.5 | 0.093 | coverage |
| P29273 | CRTI_SYNY3 (pds) | Synechocystis sp. PCC 6803 | 17.9 | 0.095 | coverage |
| P21134 | CRTI_SYNY4 | Synechocystis sp. PCC 6714 | no hit | no hit | no hit at E <= 1000 |
| Q38893 | ZDS_ARATH | Arabidopsis thaliana | 17.0 | 0.095 | coverage |
| Q9SMJ3 | ZDS_CAPAN | Capsicum annuum | 18.7 | 0.093 | coverage |
| P74306 | ZDS_SYNY3 (crtQ) | Synechocystis sp. PCC 6803 | 17.5 | 0.086 | coverage |
| Q9R6X4 | ZDS_NOSS1 (crtQ) | Nostoc sp. PCC 7120 | 24.2 | 0.088 | coverage |
| Q9M9Y8 | CRTSO_ARATH | Arabidopsis thaliana | 44.7 | 0.496 | coverage (limit 0.5) |
| P72449 | CRTU_STRGR | Streptomyces griseus | 18.6 | 0.069 | coverage |

All 10 reviewed plant and cyanobacterial PDS entries of the family fail with coverage 0.090 to 0.095 (MEASURED; the 4 in the table and Q40406, A2XDA1, Q0DUI8, C3VEP9, P49086, P28553). The hit is HMM positions 4-62 only. Of the 76 reviewed entries of the UniProt family, 54 pass the rule; 49 of the 54 are not fungal.

### Candidates that PASS and were not written

| Accession | Entry | Best domain score | HMM coverage | Why not written |
|---|---|---|---|---|
| P54980 | CRTI_CERS4 (Cereibacter sphaeroides) | 295.9 | 0.850 | Second Alphaproteobacteria CrtI. P17054 was taken. |
| P22871 | CRTI_PSEVU (Pseudescherichia vulneris) | 288.5 | 0.857 | PE 3 (homology). P21685 covers this type. |
| P54978 | CRTI_PARSN (Paracoccus sp.) | 276.1 | 0.858 | PE 3 (homology). |
| P54981 | CRTI_STRGR (Streptomyces griseus) | 261.3 | 0.853 | PE 3 (homology). The record was not read. |
| O07855 | CRTN_STAAE (Staphylococcus aureus crtN) | 325.9 | 0.857 | The EC number in the record is incomplete (1.3.8.-). The `ec` column of `provenance.tsv` needs four numbers. |
| Q4VKV1 | CRTN_METSP (Methylomonas sp crtN) | 349.1 | 0.869 | The EC number has ECO:0000305 (curator inference) only. The FUNCTION text has ECO:0000269. |
| A0A0A1GKA2 | CRTD_HALJT (Haloarcula japonica crtD, Archaea) | 331.9 | 0.857 | A candidate for a later addition: the only archaeal entry that was read. Not needed for 7. |
| P17059 | CRTD_RHOCB (Rhodobacter capsulatus crtD) | 141.5 | 0.850 | The record was not read. |

A CrtN (C30) leaf can be added if the owner accepts an incomplete EC number or curator-inference evidence.

## Tree outgroups for `cart` (5)

File: `outgroups_cart.faa`. FULL-LENGTH sequences. All record values are MEASURED. Every one has "Belongs to the carotenoid oxygenase family" in the UniProt record and ECO:0000269 on its EC number.

| Accession | Organism | Class, phylum | Gene | Type | EC | PE | Score | Function PMID | Length | Best domain score | HMM coverage |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q9HAY6 | Homo sapiens | Mammalia, Chordata | BCO1 | animal beta,beta-carotene 15,15'-dioxygenase (retinal-forming) | 1.13.11.63 | 1 | 5.0 | 11401432, 17951468, 24668807 | 547 | 132.0 | 0.628 |
| Q9JJS6 | Mus musculus | Mammalia, Chordata | Bco1 | animal BCO1 | 1.13.11.63 | 1 | 5.0 | 11092891, 15951442 | 566 | 142.9 | 0.633 |
| Q9I993 | Gallus gallus | Aves, Chordata | BCO1 | animal BCO1 | 1.13.11.63 | 1 | 5.0 | 10799297 | 526 | 133.7 | 0.628 |
| Q16518 | Homo sapiens | Mammalia, Chordata | RPE65 | animal retinoid isomerohydrolase (other function) | 3.1.1.64 | 1 | 5.0 | 25112876, 29659842 (also 16116091, 17848510, 28874556) | 533 | 118.2 | 0.625 |
| Q28175 | Bos taurus | Mammalia, Chordata | RPE65 | animal retinoid isomerohydrolase (other function) | 3.1.1.64 | 1 | 5.0 | 16096063, 19805034, 20100834 (also 15186777) | 533 | 116.2 | 0.625 |

Notes (MEASURED): the best domain of all 5 lies between HMM positions 218 and 610. A second domain of 37 to 48 bits lies between HMM positions 46 and 167. The two references score 1194.4 and 1297.4. The highest decoy scores 19.8. RPE65 is not a carotenoid cleavage enzyme in its record (EC 3.1.1.64; the Q16518 record also has EC 5.3.3.22 on one catalytic activity line). The Q28175 record has ECO:0000250 on two more catalytic activity lines.

### Candidates that FAIL the alignability rule (not written)

"All domains" = the fraction of the HMM covered by all domains of 10 bits or more together. The rule does not use it. It is given so the owner can see why these fail.

| Accession | Entry | Organism | Best domain score | HMM coverage (rule) | All domains | Fails on |
|---|---|---|---|---|---|---|
| Q9BYV7 | BCDO2_HUMAN | Homo sapiens | 79.2 | 0.187 | 0.691 | coverage |
| Q99NF1 | BCDO2_MOUSE | Mus musculus | 80.0 | 0.187 | 0.668 | coverage |
| Q9VFS2 | NINAB_DROME | Drosophila melanogaster | 53.1 | 0.179 | 0.688 | coverage |
| P74334 | ACOX_SYNY3 (ACO) | Synechocystis sp. PCC 6803 | 53.4 | 0.158 | 0.158 | coverage |
| O65572 | CCD1_ARATH | Arabidopsis thaliana | 40.5 | 0.195 | 0.409 | coverage |
| O49675 | CCD4_ARATH | Arabidopsis thaliana | 25.0 | 0.185 | 0.312 | coverage |
| Q7XJM2 | CCD7_ARATH | Arabidopsis thaliana | 38.2 | 0.184 | 0.184 | coverage |
| Q8VY26 | CCD8_ARATH | Arabidopsis thaliana | 77.3 | 0.193 | 0.501 | coverage |
| Q9LRR7 | NCED3_ARATH | Arabidopsis thaliana | 18.7 | 0.147 | 0.359 | coverage |
| O24592 | NCED1_MAIZE (VP14) | Zea mays | 26.9 | 0.187 | 0.187 | coverage |
| Q84KG5 | CCD_CROSA | Crocus sativus | 35.6 | 0.190 | 0.403 | coverage |
| Q94IR2 | CCD1_PHAVU | Phaseolus vulgaris | 35.7 | 0.188 | 0.409 | coverage |
| P9WPR5 | CCO_MYCTU | Mycobacterium tuberculosis | 39.2 | 0.195 | 0.301 | coverage |
| A4XES9 | CARDI_NOVAD | Novosphingobium aromaticivorans | 27.1 | 0.140 | 0.140 | coverage |
| Q53353 | LSDX1_SPHPI | Sphingomonas paucimobilis | 24.3 | 0.192 | 0.192 | coverage |
| S0ENQ2 | CARX_GIBF5 (fungal, reference of `carx`) | Gibberella fujikuroi | 11.2 | 0.085 | 0.085 | coverage |
| W3VHW6 | ADO_MOEAP (fungal) | Moesziomyces aphidis | 12.8 | 0.105 | 0.105 | coverage |
| Q2HHS4 | RCO1_CHAGB (fungal) | Chaetomium globosum | no hit | no hit | no hit | no hit at E <= 1000 |
| Q7S860 | CAO1_NEUCR (fungal, decoy) | Neurospora crassa | 12.8 | 0.101 | 0.101 | coverage |
| A0A0D1E6L2 | RCO1_MYCMD (fungal, decoy) | Mycosarcoma maydis | 13.0 | 0.105 | 0.105 | coverage |
| A0A384JIP2 | RCO1_BOTFB (fungal, decoy) | Botrytis cinerea | 18.8 | 0.193 | 0.193 | coverage |
| Q4W9H1 | RCO1_ASPFU (fungal, decoy) | Aspergillus fumigatus | 12.8 | 0.130 | 0.130 | coverage |
| G2QIL8 | ADO_THET4 (fungal, decoy) | Thermothelomyces thermophilus | 7.7 | 0.081 | 0.000 | score and coverage |
| Q52008 | LSDX3_SPHPI (decoy) | Sphingomonas paucimobilis | 19.8 | 0.145 | 0.145 | coverage |

Whole family (MEASURED): of the 69 reviewed entries of the UniProt carotenoid oxygenase family, 67 have a hit and 19 pass the rule. The 19 are the 2 references, S0E9A6, 4 vertebrate BCO1 and 12 vertebrate RPE65 entries. No plant entry (CCD1, CCD4, CCD7, CCD8, NCED), no bacterial entry and no fungal entry of another function passes. For every entry of the table above that has a domain of 10 bits or more, the best domain lies between HMM positions 485 and 613 (the C-terminal 20% of the HMM).

### Candidates that PASS and were not written

| Accession | Entry | Best domain score | HMM coverage | Why not written |
|---|---|---|---|---|
| Q91XT5 | BCDO1_RAT | 137.6 | 0.634 | The EC number is by similarity to Q9HAY6 (ECO:0000250). PE 2. |
| Q6PBW5 | RP65A_DANRE | 99.3 | 0.630 | The EC number is by similarity to Q28175 (ECO:0000250). PE 2. |
| 9 more RPE65 entries | Q9YGX2, Q9TVB8, Q9XT71, Q91ZQ5, O70276, Q9YI25, A9C3R8, A9C3R9, Q8AXN9 | 99.9 to 123.3 | 0.623 to 0.630 | More vertebrate RPE65 entries. Two were taken. The records were not read. |

### What else could serve as a `cart` outgroup (PROPOSED, not done)

- Accept the 5 vertebrate outgroups for rooting, and accept that no plant or bacterial leaf is in the tree.
- Change the rule to count all domains of the best hit. Then BCO2 (0.67 to 0.69), NinaB (0.69) and CCD8 Q8VY26 (0.501) would pass, and CCD1 (0.41) and ACO (0.158) would still fail (MEASURED). This changes `bin/filter_decoys_for_tree.py` and is an owner decision.
- Build the tree on an alignment that is not limited by the 2-reference HMM. This was not tested.

## mafft check (MEASURED, 2026-11-22)

`mafft --auto` (7.505) of the references with a set of candidates. "Core" = alignment columns with a residue in every reference. The table gives the fraction of core columns where the candidate has a residue, and the fraction where it has the same residue as the first reference.

| Step | Sequence | Rule | Residue in core columns | Identical to first reference |
|---|---|---|---|---|
| phytoene_desaturase (534 core columns of 892; first reference P21334) | P48537 (member) | pass | 0.996 | 0.575 |
| | P21685, P17054, P54979, Q7WT72 | pass | 0.903 to 0.948 | 0.275 to 0.360 |
| | Q64FW2 | pass | 0.919 | 0.215 |
| | Q8S4R4 | pass | 0.897 | 0.191 |
| | Q07356 (plant PDS) | fail | 0.837 | 0.140 |
| | P26294 (cyanobacterial PDS) | fail | 0.822 | 0.146 |
| cart (569 core columns of 820; first reference A1KQY3) | A7UXI1 (reference) | - | 1.000 | 0.503 |
| | Q9HAY6, Q9JJS6, Q9I993 | pass | 0.831 to 0.833 | 0.253 to 0.262 |
| | Q16518, Q28175 | pass | 0.847 | 0.241 to 0.246 |
| | Q9BYV7 (BCO2) | fail | 0.896 | 0.236 |
| | O65572 (CCD1) | fail | 0.831 | 0.211 |
| | P74334 (ACO) | fail | 0.819 | 0.188 |
| | Q9LRR7 (NCED3) | fail | 0.849 | 0.172 |
| | Q7S860 (cao-1, decoy) | fail | 0.814 | 0.158 |

Reading: the kept outgroups are not mostly gaps (0.83 to 0.95 of the core columns). The failing candidates are not mostly gaps either (0.81 to 0.90). Gap content in a mafft alignment does not tell the two groups apart. The Q02861 outgroup was not in this alignment. The alignment files were not kept.

## Who reads this folder (checked again by grep on 2026-11-22)

`grep -rn "tree_members"` over `*.py`, `*.sbatch`, `*.sh`, `*.nf`, `*.yaml`, `*.yml`, `*.toml`, `*.config` under `analysis/pathway_census/`, outside this folder: 2 files.

- `slurm/step5_prepare.sbatch`, 1 line: a comment for the option `PC_MEMBERS_PROV` (the provenance table gives the class of each named member). This is the tree step. It sets no cutoff and makes no absence call.
- `tests/test_tree_members.py`.

No file in `bin/`, `pathway_census/` or `slurm/` other than the step 5 tree code names the folder. No file in `bin/`, `pathway_census/` or `slurm/` has `rglob` or a `**` glob (MEASURED). `tests/test_tree_members.py` now checks every `.py`, `.sbatch` and `.sh` file of these three folders, with `slurm/step5_prepare.sbatch` and `bin/step5_place.py` as the only allowed readers.

The risk stated above stays: `phytoene_desaturase.faa` in this folder has a step name. Never pass this folder as `--refs-dir`.

## Files of this addition

| File | Content |
|---|---|
| `phytoene_desaturase.faa` | 1 full-length named member (P48537). Header: accession, `region=1-621`, then the UniProt header. |
| `outgroups_phytoene_desaturase.faa` | 7 full-length outgroup sequences. |
| `outgroups_cart.faa` | 5 full-length outgroup sequences. |
| `cart.faa` | Not written. There is no named member. |
| `full_length/named_members.faa` | Now 14 records: the 13 above and P48537. The first 13 records are byte-identical to the earlier file (MEASURED). |
| `provenance.tsv` | Now 46 rows. Rows 1-33 are byte-identical to the earlier file (MEASURED). 13 new rows: 1 `named_member` and 7 `tree_outgroup` for `phytoene_desaturase`, 5 `tree_outgroup` for `cart`. The `reason` text of each new row gives the best domain score and the HMM coverage. |
| `build_tree_members.py` | Extended with `second_set()`. One run writes all files of both sets. The measured scores are constants in the script; when `hmmsearch` and the HMMs are present the script measures again and stops on a difference. A second run gave the same bytes (md5, MEASURED). |

The counts "33 rows" and "13 full-length sequences" in the "Files" table above describe the first set only.

Not done here: no tree, no change to any reference, decoy, YAML, HMM, cutoff or tool, and no paper was read.
