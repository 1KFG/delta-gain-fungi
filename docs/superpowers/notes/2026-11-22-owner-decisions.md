# Owner decisions at the Stage 1 review gate (2026-11-22)

Source: the owner's reply in the dgf-pathways session. Packet: `docs/superpowers/plans/2026-11-21-stage1-review-packet.md`.

## Decided

1. **Rule 1: strict reading.** The tables stay as they are. The quoted sentence must name the protein of the link and state an assay or a product. Applies to `carotenoid`, `melanin_dhn` and `melanin_ldopa`. The context reading stays recorded as the alternative (DC1, DD5).
2. **Gene-name rule: accepted.** A reference whose UniProt entry has no gene name may use the protein abbreviation of the paper, or the ORF name. The name is valid only together with its organism. If a name is ambiguous, prefix it with a species prefix. Covers DD1, DL1 and DL2.

3. **`min_score` rule (step 2): accepted.** Multi-reference steps use the lowest leave-one-out score. Single-reference steps use the highest decoy score as a floor. Tree placement (step 5) makes the real calls. Basis: `results/step2/carotenoid/` (MEASURED, 2026-11-22).
4. **Lycopene cyclase decoys: add them.** Other-family cyclases as outgroups. A Fable agent curates them (dispatched 2026-11-22).
5. **Pilot order:** carotenoid first, then the other three pathways. Work in the worktree until step 12. Run without the contamination reference set first (raw counts). Tyrosinase and laccase get separate trees with decoys inside.

## Deferred (current tables stand)

DD3, DD4, DL3, DL4, DL5, DL6, the mycosporine decisions (`refs/mycosporine/UNVERIFIED.md`), and carotenoid items after DC1.

## Names chosen by the owner (2026-11-22, may change with new information)

| Accession | Name to use |
|---|---|
| Q9HFV6 | 4HNR |
| Q2TL94 | Q2TL94_tyrosinase_Psang |
| A0A1S9DK56 | OAory_01107480 |
| Q12634 | MGG_02252_DHN_Pgris |

Pattern for a name that needs disambiguation: `<accession or locus>_<role>_<species code>`.
These are working labels. A later step may replace them with one common name per ortholog group.

## Planned: ortholog-based naming (owner direction, not yet done)

After the pathway census, build phylogenies of candidate genes, match them to named copies, and re-cluster gene sets that fall in one ortholog group under a common name.

## Entries covered by decision 2

| Accession | Organism in the table | Name in the table | Pathway | Why the name is open |
|---|---|---|---|---|
| Q9HFV6 | Pyricularia grisea | 4HNR | melanin_dhn | UniProt has no gene name. 4HNR is the paper's abbreviation (PMID 10956664). |
| Q2TL94 | Trametes sanguinea | tyrosinase | melanin_ldopa | UniProt has no gene name. "tyrosinase" is the protein name. The paper says Pycnoporus sanguineus. |
| A0A1S9DK56 | Aspergillus oryzae | OAory_01107480 | melanin_ldopa | UniProt has the ORF name only. The name melB is not in the entry. melB is also a synonym of tyrP in Aspergillus terreus (Q0CRX0). |
| Q12634 | Pyricularia grisea | MGG_02252 | melanin_dhn | ORF name. The UniProt protein name belongs to the other DHN reductase. |

## Follow-up (not done)

- The verifier reports 6 failures, all for Q9HFV6 and Q2TL94 (`ref_gene` and `link_gene`). Each reads "no match; names in record: (empty)".
- Next step: add a reviewed per-entry name alias in `bin/verify_kb.py` for these two entries.
- Expected result: 0 failures. Not yet run.

## Lycopene cyclase decoys (2026-11-22)
- Q9HNE5 (Halobacterium crtY): tree outgroup only, not a decoy (`refs/carotenoid/outgroups/`). BLAST bit score against the best fungal reference 87.0 over 192 aligned residues, 33% identity (MEASURED).
- Q9KK78 and Q9KK79: owner rule is to move them to outgroups only if they could cross-react in detection. MEASURED: HMM scores 29.3 and 2.4 (cutoff floor is 174.1, the lowest held-out score). BLAST bit scores against the fungal references 34.3 and 21.9 over 39 and 27 residues (E 5e-6 and 0.075). They do not reach the cutoff, so they stay decoys.
- CruA, CruP and Q38932 (epsilon cyclase) are not added: no good evidence found.

## Step 3 and 4 decisions (2026-11-22)
- Chunk size 4,000 genomes, 6 jobs (owner). Evidence: `2026-11-22-step4-results-and-decisions.md`.
- Lowered cutoffs for `astaxanthin_synthase` (999.4) and `ylo1` (480.5) kept for step 5 (owner), after the side evaluation. Evidence and files: same note.

## BUSCO filter for the gain/loss analysis (2026-11-23)
- Owner decision: no BUSCO completeness filter in the main results. Reason: the BFD BUSCO run used `fungi_odb12` for every genome and the budding yeasts score low as a group (Saccharomycetes median 75.4%, 10th to 90th percentile 72.2 to 77.9, 2,589 genomes), which may be a marker-set effect and not poor assemblies (UNVERIFIED). A filter at 80 removes 27 of 43 Saccharomycetes genus tips.
- The runs with cutoffs 65 and 70 stay as a sensitivity check: counts change by 13% or less and the core co-loss result is unchanged (report section 6h; `results/gainloss/carotenoid_busco_compare/`).
- A Saccharomycotina-specific BUSCO set is a note for the next BFD BUSCO run, not run now (`Fungi_BFD_runs/docs/NOTES_TO_FOLLOW.md`, open question Q2).
