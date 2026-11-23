# Decoys of the `lycopene_cyclase` step (carotenoid pathway)

Written 2026-11-22. Curator: Claude (agent). No person reviewed this list. Revised 2026-11-22 after the owner decision on Q9HNE5: 9 decoys and 1 tree outgroup (the first version had 10 decoys).

Labels: MEASURED = read from a fetched record or computed on 2026-11-22. PROPOSED = the choice of the curator. UNVERIFIED = not checked.

## Limits

- The decoy list is PROPOSED. The owner approved the idea of decoys from other cyclase families. The owner did not review these nine entries.
- Owner decision (2026-11-22): Q9HNE5 (*Halobacterium salinarum* crtY) is a true homolog that does the reaction of the step. It is a tree outgroup only. It is not a decoy, because a true homolog in the decoy set would wrongly raise the decoy-calibrated score cutoff. It is in `refs/carotenoid/outgroups/`.
- Seven of the nine decoys have `verified: false`. See "The `verified` flag" below. This is an owner decision.
- Three decoys are unreviewed TrEMBL entries (Q9KK78, Q9KK79, Q93QX4). Their UniProt records state no function. The function comes from the PubMed abstract only. No full text was read.
- For the four PF05834 decoys no paper was read. The function comes from the reviewed UniProt record.
- Q9KK78 and Q9KK79 are the two subunits of one heterodimeric lycopene cyclase that does the reaction of the step (abstract of PMID 10821176). The owner kept them as decoys. They are decoys by lineage only. The four PF05834 decoys (except Q42435) also do the reaction of the step (EC 5.5.1.19) with a protein of another family. Whether the reasoning that moved Q9HNE5 out also applies to these entries needs an owner decision.
- No tree was built. The statement "the tree must separate it" is a requirement for a later step, not a result.
- `refs/carotenoid/UNVERIFIED.md` section F item 2 now describes the 9 decoys and the outgroup.
- `curation_scripts/carotenoid/build.py` holds the 9 decoys and the outgroup. A re-run reproduces the committed files (see "build.py").

## What was done

1. UniProt searches (REST API, through `pathway_census.kbsources`, 2026-11-22). MEASURED counts:
   - `reviewed:true AND xref:ncbifam-TIGR03462`: 23 entries. 19 are fungal fused proteins. 4 are not fungal: Q9HNE5 and B0R753 (*Halobacterium salinarum* crtY, two strains), A1XI29 and A1XI30 (*Dietzia* sp. CQ4 lbtA, lbtBC).
   - `reviewed:true AND xref:pfam-PF18916`: 5 entries (P37295, D1Z4K7, B0R753, Q9HNE5, A1XI30).
   - `reviewed:true AND xref:pfam-PF05834`: 21 entries.
   - `gene_exact:crtYc OR gene_exact:crtYd`: 19 entries, all unreviewed.
   - `gene_exact:crtYe OR crtYf OR crtYg OR crtYh`: 13 entries, all unreviewed.
   - `(gene_exact:cruA OR gene_exact:cruP) AND taxonomy_id:2`: 13 entries, all unreviewed.
2. Each chosen accession was fetched again as a full UniProt JSON record. The script checked: accession, gene name, organism lineage (no "Fungi"), family cross-references, EC numbers, the cited PMID, and sequence length.
3. Four PubMed abstracts were fetched and read: PMID 12003928, 10821176, 17008032, 11432736.
4. Homology to the fungal cyclase domain was measured in three ways (HMMER 3.4, MAFFT 7.505):
   - `hmmsearch --max -E 1000 --domE 1000` with `results/hmms/carotenoid/lycopene_cyclase.PROVISIONAL.hmm` (241 match states, built from the five references) against 17 candidate sequences. E-values use Z = 17. They are not database-size E-values.
   - `phmmer --max -E 1000` with the five reference records as queries against the same 17 sequences. The best score of each candidate is reported.
   - `mafft --localpair --maxiterate 1000` of the five references plus one decoy. Identity = identical residues / columns where both sequences have a residue. Reported for the best reference. Done only for the seven TIGR03462 candidates. At identity below about 30 % the alignment is not reliable; read these values as rough.

## The nine decoys

All values in the table are MEASURED, from the first version of this note. The Q9HNE5 row was removed from the table; its numbers are in the next section. "Step 2 score" is the best domain score in the step 2 run (job 32374977). "No hit" means no hit at `-E 1000`.

| Accession | Organism (UniProt) | Gene | Entry | Family ids in the record | Function, and source | hmmsearch domain score | phmmer best score (reference) | Identity to best reference |
|---|---|---|---|---|---|---|---|---|
| Q9KK78 | Brevibacterium linens | crtYc | TrEMBL | TIGR03462, PF18916 | Subunit of a heterodimeric lycopene cyclase. Abstract of PMID 10821176 | 29.3 | 28.7 (Q7Z859) | 27.3 % (33 of 121, Q7Z859) |
| Q9KK79 | Brevibacterium linens | crtYd | TrEMBL | TIGR03462, PF18916 | Subunit of a heterodimeric lycopene cyclase. Abstract of PMID 10821176 | 2.4 | 3.2 (Q8X0Z1) | 26.6 % (25 of 94, Q9UUQ6) |
| A1XI29 | Dietzia sp. (strain CQ4) | lbtA | Swiss-Prot | TIGR03462 | Putative C(50) carotenoid beta-cyclase subunit A, EC 5.5.1.-. UniProt record; abstract of PMID 17008032 | 10.1 | 15.9 (Q9UUQ6) | 21.6 % (25 of 116, Q9UUQ6) |
| A1XI30 | Dietzia sp. (strain CQ4) | lbtBC | Swiss-Prot | TIGR03462, PF18916, PF01040, NF009608 | C(50) cyclase subunit fused to a lycopene elongase/hydratase. UniProt record; abstract of PMID 17008032 | 1.6 | 0.7 (Q9P854) | 22.0 % (54 of 245, Q7Z859) |
| Q93QX4 | Corynebacterium glutamicum | crtYe | TrEMBL | TIGR03462 | Subunit of the C(50) carotenoid epsilon-cyclase (flavuxanthin to decaprenoxanthin, with CrtYf). Abstract of PMID 11432736 | 12.9 | 19.5 (Q9UUQ6) | 29.2 % (33 of 113, Q9UUQ6) |
| P21687 | Pantoea ananas | crtY | Swiss-Prot | PF05834, TIGR01790, TIGR01789 | Lycopene beta-cyclase, EC 5.5.1.19. UniProt record | no hit | -1.7 (Q8X0Z1) | not computed |
| Q55276 | Synechococcus elongatus (strain ATCC 33912 / PCC 7942 / FACHB-805) | crtL | Swiss-Prot | PF05834, TIGR01790, NF045687 | Lycopene beta cyclase, EC 5.5.1.19. UniProt record | no hit | -0.6 (Q7Z859) | not computed |
| Q38933 | Arabidopsis thaliana | LCY1 | Swiss-Prot | PF05834, TIGR01790 | Lycopene beta cyclase, chloroplastic, EC 5.5.1.19. UniProt record | no hit | -1.3 (P37295) | not computed |
| Q42435 | Capsicum annuum | CCS | Swiss-Prot | PF05834, TIGR01790 | Capsanthin/capsorubin synthase, EC 5.3.99.8. UniProt record | no hit | -0.9 (Q9UUQ6) | not computed |

## Homology: what the measurements show

- **Q9HNE5 is a homolog of the fungal cyclase domain (tree outgroup, not a decoy).** MEASURED: Swiss-Prot, TIGR03462 and PF18916 in the record, domain score 115.6, alignment over HMM positions 35 to 237 of 241, phmmer best score 79.0 (Q9P854), 32.3 % identity to Q7Z859. It does the reaction of the step. The tree must place it outside the fungal proteins.
- **Q9KK78 is a homolog of part of the domain.** MEASURED: domain score 29.3, HMM positions 48 to 91.
- **Q93QX4 and A1XI29 give weak hits** (12.9 and 10.1). **Q9KK79 and A1XI30 give no usable hit with this HMM** (2.4 and 1.6). Their family membership comes from the UniProt cross-reference to TIGR03462, not from this measurement. Homology of these two to the fungal domain is UNVERIFIED by the measurements here.
- **The four PF05834 proteins give no hit** with the step HMM at `-E 1000` and phmmer scores below 0. Homology to the fungal domain is not detected. They test the absence check, not the tree: a method that calls them `lycopene_cyclase` is wrong.
- The abstract of PMID 10821176 (year 2000) says that CrtYc and CrtYd "show no similarity to any of the known lycopene cyclases". The abstract of PMID 12003928 (year 2002) names "a plausible evolutionary connection between heterodimeric lycopene cyclases in bacteria and bifunctional lycopene cyclase-phytoene synthases in fungi". Both sentences are MEASURED (read in the fetched abstracts).

## Step 2 result (job 32375043, 2026-11-22, MEASURED; nine decoys)

`results/step2/carotenoid/step_scores.PROVISIONAL.tsv`, row `lycopene_cyclase`:

| n_refs | ref_min | ref_max | proposed_min_score | n_decoys_at_or_above | max_decoy_score | n_other_step_refs_at_or_above | loo_min |
|---|---|---|---|---|---|---|---|
| 5 | 347.5 | 384.2 | 347.5 | 0 | 29.3 | 0 | 174.1 |

- No decoy scores at or above 347.5 (lowest own-reference score). 0 of 9.
- The highest decoy score is 29.3 (Q9KK78, same value as the hmmsearch column above). The lowest held-out score is 174.1. The gap is 144.8 bits.
- Before the owner decision (job 32374977, ten decoys) `max_decoy_score` was 115.6. That value came from Q9HNE5, which is no longer a decoy. The proposed minimum did not change (347.5), because it is the lowest own-reference score.
- Only Q9KK78 scores above 13 among the nine decoys. The decoy set tests the cutoff weakly. A cutoff anywhere between 29.3 and 174.1 is consistent with the decoys. None of these values is tested on a genome.
- `targets_at_or_above_min.PROVISIONAL.tsv` has a header and no rows.
- The other 13 steps did not change (diff of `step_scores.PROVISIONAL.tsv` and `leave_one_out.PROVISIONAL.tsv` against the previous run: only `max_decoy_score` of `lycopene_cyclase` differs).
- Limit: the HMM is PROVISIONAL and comes from five sequences. The held-out minimum 174.1 comes from five leave-one-out runs.

## The `verified` flag (owner decision)

The header of `pathways/carotenoid.yaml` defines `verified: true` for a decoy: the entry was fetched, organism and gene name matched, the entry carries the family of the step, and its reviewed UniProt record states a function that is not the function of the step.

All nine decoys and the outgroup were fetched and organism and gene name matched (MEASURED). Under the full rule:

| Decoys | verified | Why |
|---|---|---|
| A1XI29, A1XI30 | true | Reviewed. TIGR03462 in the record. No EC 5.5.1.19. The record states a C(50) function. |
| Q9KK78, Q9KK79, Q93QX4 | false | Unreviewed entries. |
| P21687, Q55276, Q38933, Q42435 | false | Not in the family of the step. |
| Q9HNE5 (outgroup, in `outgroups/provenance.tsv`) | false | Reviewed, in the family, but it does the reaction of the step. The flag is kept from the first version. |

The rule was written for decoys of the same family with another function. The owner asked for decoys of other families, which the rule does not cover. The curator did not change the rule and did not set `true` where the rule does not hold. The owner can widen the rule, or accept `false` for these entries.

## Tests

The first version of this note recorded 2 failing tests: `test_unverified_refs_are_listed_in_unverified_md` and `test_committed_packet_is_current`. The revision fixed both by editing `refs/carotenoid/UNVERIFIED.md` and by making the review packet again. No test and no validator was changed.

## build.py

`curation_scripts/carotenoid/build.py` holds the nine decoys (`DECOY["lycopene_cyclase"]`), the per-decoy `retrieved_on` and `verified` values (`DECOY_META`) and the outgroup (`OUTGROUP`). It writes `decoys/*`, `outgroups/*` and the yaml. A re-run reproduced the files byte for byte (`diff -r` against a backup).

Family check for the decoys of this step: the UniProt record must carry NCBIfam TIGR03462 (`DECOY_XREF`). A Pfam-only check does not fit, because A1XI29 and A1XI30 have no Pfam PF18916 in UniProt. The four PF05834 decoys have an explicit allowance (`DECOY_OTHER_FAMILY`): the record must carry PF05834 and must not carry TIGR03462. The check of the 13 other steps did not change.

## Candidates that were fetched and not added

| Accession | Organism | Gene | Why not added | hmmsearch domain score (MEASURED) |
|---|---|---|---|---|
| B0R753 | Halobacterium salinarum (strain R1) | crtY | Second strain of the same protein as Q9HNE5. Seen in the search list only; the full record was not fetched. | not run |
| Q93QX3 | Corynebacterium glutamicum | crtYf | Partner of Q93QX4, unreviewed. One subunit of the pair is enough. | 4.5 |
| Q38932 | Arabidopsis thaliana | LUT2 | Lycopene epsilon cyclase (EC 5.5.1.18), PF05834. Q38933 and Q42435 already cover plant PF05834. A good addition if the owner wants an epsilon cyclase. | no hit |
| Q8NSP8 | Corynebacterium glutamicum (strain ATCC 13032) | crtEb | Lycopene elongase/hydratase (EC 2.5.1.149), PF01040. Not a cyclase. A1XI30 already carries this domain. | -3.4 |
| M0L7V9 | Haloarcula japonica | lyeJ | Lycopene elongase/hydratase (EC 2.5.1.150), PF01040. Not a cyclase. | -3.3 |
| Q8KF75 | Chlorobaculum tepidum (strain TLS) | CT0456 | CruA candidate. Unreviewed, no gene name "cruA" in the record, no family cross-reference, cites the genome paper only. | no hit |
| A0ACS0KHX7 | Picosynechococcus sp. (strain PCC 7002) | cruA | Unreviewed, entry version 1, no family cross-reference, no paper with a PMID. | no hit |
| A0ACS0KCE4 | Picosynechococcus sp. (strain PCC 7002) | cruP | Same as above. | no hit |

CruA and CruP: no reviewed entry exists (MEASURED, search above). No entry ties the sequence to a paper on the function. For this reason no CruA or CruP decoy was added. The owner can ask for one with `verified: false`.

## Files

- `refs/carotenoid/decoys/lycopene_cyclase.faa`: nine full-length sequences, header style of the other decoy files.
- `refs/carotenoid/decoys/provenance.tsv`: nine new rows at the end. The 60 old rows are unchanged.
- `refs/carotenoid/outgroups/lycopene_cyclase.faa` and `outgroups/provenance.tsv`: Q9HNE5 (same columns as the decoy provenance, plus `role` = `tree_outgroup`). The absence check and the decoy-calibrated `min_score` do not read this directory.
- `pathways/carotenoid.yaml`: a `decoys:` block with nine entries in the `lycopene_cyclase` step. No other line changed. The schema has no outgroup key.
- `results/step2/carotenoid/*`: made again by job 32375043.
- The exploration scripts of the first version were run from a scratch directory and are not in the repository.
