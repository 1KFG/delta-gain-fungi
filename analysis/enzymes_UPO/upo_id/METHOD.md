# UPO identification: the UPObase method and our replication

Written 2026-11-07. Tool version 0.1.0.

Labels used on every fact:
- **[paper]** stated in a cited paper (quote or section given).
- **[measured]** measured by me in this folder on 2026-11-07 (file named).
- **[assumption]** a choice I made without data that proves it is right.

## 1. Sources

| Short name | Citation | Access |
|---|---|---|
| Faiza 2019a | Faiza M, Huang S, Lan D, Wang Y (2019). New insights on unspecific peroxygenases: superfamily reclassification and evolution. *BMC Evol Biol* 19:76. doi:10.1186/s12862-019-1394-3 | PMC6417270, full text and supplementary Tables S1, S4 read |
| Faiza 2019b (UPObase) | Faiza M, Lan D, Huang S, Wang Y (2019). UPObase: an online database of unspecific peroxygenases. *Database* 2019:baz122. doi:10.1093/database/baz122 | PMC6902001, full text read |
| Kinner 2021 | Kinner A, Rosenthal K, Lütz S (2021). Identification and expression of new unspecific peroxygenases – recent advances, challenges and opportunities. *Front Bioeng Biotechnol* 9:705630. doi:10.3389/fbioe.2021.705630 | PMC8293615, full text read |
| Structures | PDB 2YOR (AaeUPO), 1CPO (LfuCPO), 5FUJ (MroUPO), 7O1R (HspUPO), 7ZCL (CviUPO) | RCSB PDB files; LINK, SSBOND, DBREF records read |

The UPObase site (upobase.bioinformaticsreview.com) does not resolve from HPCC (REFERENCES.md section 5).
Both papers are online-only articles without page numbers. I cite section names and table numbers instead.

Follow-up papers found by a Europe PMC citation search (6 papers cite Faiza 2019b; 34 cite Faiza 2019a).
None of them publishes a new UPObase release. Two are relevant to classification. I read only their abstracts:
- Zámocký M, Harichová J (2022). Evolution of heme peroxygenases: ancient roots and later evolved branches.
  *Antioxidants* 11:1011. Abstract: HTP genes "are distributed also in the clades of Mucoromycota and
  Chytridiomycota with deep ancient evolutionary origins".
- Zámocký M et al. (2025). Molecular diversity of unspecific heme peroxygenases... *Int J Biol Macromol*,
  doi:10.1016/j.ijbiomac.2025.147823. Abstract: "UPOs from the phyla of Mucoromycota, Glomeromycota and
  Chytridiomycota remain not sufficiently investigated". This verifies the claim listed as unverified in PLAN.md
  section 10 (abstract only; full text not read).

## 2. The UPObase method as published

UPObase was built in two stages. Stage 1 is Faiza 2019a (113 UPOs). Stage 2 is Faiza 2019b (1,948 sequences).

### 2.1 Stage 1: Faiza 2019a (BMC Evol Biol)

1. **Database.** [paper] "Peptide sequences of all fungal genomes constituting 812 different species (or strains)
   were downloaded from the Ensembl fungal genome database" (Methods: Data collection). Genome list = Additional file
   16 (Table S3). CPO outgroup sequences "were downloaded from NCBI" (Table S4, 31 CPOs + MroUPO 5FUJ).
2. **Search.** [paper] "AaeUPO was used as a query to perform a similarity search using PHMMER (version 3.1b2) ...
   with an E-value of 10.0 and an inclusive E-value of 0.01" (Methods: Homology search and identification of UPO
   motifs).
3. **Clustering.** [paper] cd-hit "at 90% similarity cut-off and a word-length of 5 residues", then "MCL ... at an
   inflation value of 1.4 to eliminate the dissimilar sequences corresponding to the seed sequence" (same section).
   The paper does not say which MCL clusters were kept.
4. **Motif filter.** [paper] "The resultant sequences were then searched for motifs providing the final resultant
   UPO-encoding sequences" (same section). Required motif: [paper] "PCP---EGD---R----E" (Results: Conserved motifs;
   Fig. 4). The paper also names "S[IL]G" (between PCP and EGD) and "SXXRXD" (after EGD) as present in "virtually all
   UPOs", and "six amino acids in between the acid-base catalyst pair".
5. **Result.** [paper] "113 putative UPO sequences ... in 35 different fungal species" (Abstract; Table S1 gives
   counts per species). One sequence per species (35, Table 3, with GenBank/Ensembl IDs) went into the phylogeny.
6. **Phylogeny.** [paper] MEGA7, ProtTest3 model LG+G+F, ML, 1000 bootstraps (Methods: Phylogeny analysis).
7. **Classes.** [paper] Five UPO subfamilies (I–V) with signature motifs, and a "Peroxidase-peroxygenase (Pog)
   superfamily" that holds MroUPO and LfuCPO. "Classic CPOs" form a separate outgroup family. Eight NCBI "CPO"
   sequences fall among the UPOs ("8 CPO sequences now classified as a kind of UPOs", Fig. 3 legend), some with
   "EAD/ETD motif ... instead of EGD".
8. **Catalytic motif labels.** [paper] "conserved motif patterns for the catalytic activity of AaeUPO is
   -PCP-EGD-R--E, and for MroUPO and CPO is -PCP-EHD-E" (Introduction). Key residues named: AaeUPO Cys36 (heme
   ligand), Glu196 and Arg189 (acid-base pair), Glu122-Gly123-Asp124 (EGD); LfuCPO Cys29, His105, Glu183.

### 2.2 Stage 2: Faiza 2019b (UPObase, Database)

1. **Seeds.** [paper] "A set of previously found 113 UPO encoding sequences belonging to different subfamilies and a
   superfamily were used" (Methods: Database construction). The 113 sequences are not listed in either paper.
2. **Search.** [paper] Each seed queried by "PHMMER (version 3.1b2) against the generated fungal genome database with
   an E-value and an inclusive E-value set to 10.0 and 0.01". Same 812 Ensembl Fungi genomes.
3. **Clustering and motifs.** [paper] cd-hit 90 % (word 5), MCL inflation 1.4, then "searched for sequence motifs
   corresponding to their subfamily type". Repeated "for each new sequence that appeared in the similarity search".
4. **Dereplication.** [paper] cd-hit at 95 % gave "a total of 1948 clusters"; representatives form the database.
5. **Phylogeny.** [paper] MEGA7, WAG+G+F (ProtTest3), ML, 300 bootstraps (Methods: Phylogenetic analysis).
6. **Classes and counts.** [paper] (Results: Database sequences; Table 1)
   - UPO subfamilies: Ia 54, Ib 16, II 3, IIIa 29, IIIb 5, IV 6, V 10 (total 123).
   - Pog superfamily 1,825: Pog-a 90, Pog-b 47, Pog-c 26, Pog-d 17, Pog-e 8, Pog-f 128, Pog-g 1,509.
   - "Pog-g sequences do not have any signature motif pattern of their own except the Cys ligation to the heme".
7. **Motif strings.** [paper] Table 1 lists one ordered motif string per class (copied verbatim into
   `upo_id/upo_id.py`, `UPOBASE_MOTIFS`). Gap lengths between motifs are not given.
8. **Length rule.** [paper] Only for user submissions: "must be longer than 100 amino acid residues".
9. **Validation.** [paper] No benchmark, no negative set and no error rate are reported in either paper. The only
   check is phylogeny plus motif agreement.

### 2.3 Kinner 2021 Family I / Family II (primary scheme for this project, user decision)

[paper] "Family I, which includes 'short-type' UPOs, and family II containing 'long-type' UPOs (Kellner et al.,
2014; Hofrichter et al., 2015)". Family I: ">1,100 short-type UPOs, occurring in most fungal phyla, with a mean size
of around 30 kDa, but mainly without predictable signal peptides"; CglUPO, MroUPO, MweUPO and LfuCPO are "assigned to
this short-type family". Family II: ">900 long-type UPOs", "mean size of 44.4 kDa", "only ... in Ascomycota and
Basidiomycota", "most sequences possess a C-terminal disulfide bridge and predicted signal peptides". Motifs:
"-EHD-S-E- and -EGD-S-R-E- are found in family I and II, respectively, along with the characteristic PCP motif".

## 3. What the structures say about EHD vs EGD [measured]

PDB LINK records show that the E-X-D motif is a metal site. The Glu side chain, the X backbone O and a Ser side
chain three residues later bind Mg2+ (Mn2+ in 1CPO):

| Enzyme (PDB) | Heme Cys | Metal site | Acid-base Glu | Arg | C-terminal S-S |
|---|---|---|---|---|---|
| AaeUPO (2YOR, mature numbering) | C36 | E122 G123 (S126) | E196 | R189 | C278–C319 |
| LfuCPO (1CPO, mature) | C29 | E104 H105 (S108) | E183 | – (His105) | – (C79–C87 internal) |
| MroUPO (5FUJ) | C17 | E85 H86 (S89) | E157 (by alignment) | – | – (C227 inter-chain) |
| HspUPO (7O1R) | C39 | E109 H110 (S113) | E180 (by alignment) | – | – |
| CviUPO (7ZCL) | C19 | E89 H90 (S93) | E162 (by alignment) | – | – |

Result: three characterized short UPOs (MroUPO, HspUPO, CviUPO) have EHD, like LfuCPO. So "EHD" does not mean
"not a peroxygenase". Kinner's label (EHD = Family I) and Faiza's label (EHD = CPO-type) describe the same split;
they differ only in whether the EHD group is called "UPO". This resolves conflict 3 in REFERENCES.md.

The reference alignment (`models/ref_alignment.afa`, MAFFT L-INS-i) places all these residues in the same columns
as the structures predict (`models/posmap.tsv`). That is 5 of 5 structures, 4 sites each.

## 4. Replication design

### 4.1 Reference set (`ref/`)

`ref/manifest.tsv` lists every sequence with accession, source URL, original header, length, SHA-256 and fetch date.
92 sequences, fetched 2026-11-07 from UniProt REST, NCBI E-utilities, Ensembl REST and RCSB (`scripts/fetch_refs.py`;
one request at a time, back-off on HTTP 429).

| Group | n | Source |
|---|--:|---|
| Characterized enzymes | 13 | REFERENCES.md section 6. HspUPO: UniProt A0A1Y2TH07 is now inactive (demerged), so PDB 7O1R was used. MroUPO, CviUPO, DchUPO, DspUPO, HspUPO from PDB. |
| UPObase stage-1 seeds | 35 | Faiza 2019a Table 3 (one per species of the 113). |
| Faiza CPO outgroup | 31 | Faiza 2019a Table S4. |
| Negatives | 13 | 5 P450s, 7 classical/other peroxidases (CCP, HRP, LiP, MnP, ARP, KatG, DyP), 1 vanadium CPO. |

[measured] The GenBank lengths of the 35 seeds do not match Table 3 lengths (for example Mixia osmundae GAA97831 is
454 aa; Table 3 says 189). Table 3 lengths probably describe the SwissModel models, not the proteins
[assumption]. Species names match for all 35.

### 4.2 Training sets (`ref/train_sets.tsv`)

I grouped the 79 positive references by three measured features: metal-site middle residue, C-terminal Cys pair at
AaeUPO C278/C319, and clade in a FastTree (LG) tree of the reference alignment (no bootstrap; [assumption] that the
clades are stable).

- **FamII (47):** AaeUPO, CraUPO, PabUPO-II, UPO1, all 35 seeds, and the 8 "CPO" sequences that Faiza 2019a moved
  into the UPOs. [measured] All 47 have E[GAT]D. 36 of 47 have both Cys at the AaeUPO C278/C319 columns; the 11
  without are CraUPO (fragment), 9 Ustilaginomycotina/Tremellomycetes seeds and Acidomyces. Placing the 8 CPOs here
  matches Faiza's statement that 8 CPOs belong with the UPOs.
- **FamI_UPO (8):** the characterized short UPOs MroUPO, MweUPO, TteUPO, DcaUPO, HspUPO, CviUPO, DchUPO, DspUPO.
  [measured] All 8 have EHD and lengths 227–279 aa.
- **CPO (16):** LfuCPO and Faiza CPOs in the clade that holds LfuCPO.
- **Holdout (8):** sequences with unstable placement, truncation or a duplicate. Not used for training.

### 4.3 Classification steps (`upo_id/upo_id.py classify`)

| Step | What it does | Faithful to UPObase? |
|---|---|---|
| a1 | `hmmsearch --cut_ga` with PF01328 (Pfam 38.2, PF01328.23) | Deviation: UPObase used phmmer only. PF01328 is the BFD family net. |
| a2 | `hmmsearch` with FamII, FamI_UPO and CPO HMMs (hmmbuild from 4.2). A protein is a candidate if PF01328 GA hits or any family HMM has full-sequence E ≤ 1e-5 | New. E ≤ 1e-5 is [assumption]. |
| a3 | `phmmer -E 10 --incE 0.01` with 38 seeds (35 Table 3 seeds + AaeUPO, MroUPO, LfuCPO); reported, not used in the decision | Same tool and cutoffs as the paper; seed set incomplete (38 of 113+). |
| b | `hmmalign` of candidates to an anchor HMM. The anchor HMM is built with `--hand`; match state k = AaeUPO precursor residue k (mature = k − 43). Residues are read at fixed sites: C36 (PCP), S94-I95-G96 (SIG), E122-X123-D124 + S126 (metal site), R129/D131 (SXXRXD), R189, E196, C278, C319 | Deviation: alignment-based; the paper does not say how motifs were searched. |
| b2 | Acid-base Glu = E at the AaeUPO E196 column, or an R-X(6,7)-E match in the query segment aligned to AaeUPO 137–219 (mature) | Needed because ascomycete long seeds have A at the E196 column but carry R-X(6,7)-E nearby [measured]. Faiza 2019a reports 6 or 7 residues between the pair. |
| c | Length (≥ 330 aa = long), C-terminal Cys pair, signal peptide from an optional table | 330 aa: [assumption], calibrated on references (short ≤ 308, long ≥ 352 except the CraUPO fragment). |
| d | Decision table (4.4) | Kinner scheme, not UPObase. |
| e | UPObase Table 1 motif strings on the raw sequence, motifs in order, any gap length, earliest match; the class with the most motif elements wins; Pog-g if no class matches but "PCP" occurs | Replication of Table 1. Gap rule, tie rule and the PCP test for Pog-g are [assumption]. |
| f | UPObase catalytic label from the alignment: PCP-EGD-R-E or PCP-EHD-E or neither | From Faiza 2019a/b text. |
| g | Nearest characterized enzyme: Biopython local alignment (BLOSUM62, gap −11/−1) to the 13 characterized enzymes; % identity = identities / (identities + mismatches + gaps) | New. |

### 4.4 Decision table

| Order | Condition | final_class |
|---|---|---|
| 1 | no PF01328 GA hit and no family HMM hit at E ≤ 1e-5 | `no_HTP_hit` (rows only with `--all`) |
| 2 | core coverage (AaeUPO 17–237 mature) < 0.5, or no Cys at C36, or no Glu at E122, or no acid-base Glu (b2) | `HTP_fails_motif` (reason lists each failed check) |
| 3 | best family HMM bit score = FamII | `UPO_FamilyII` (Kinner II) |
| 4 | best = FamI_UPO | `UPO_FamilyI` (Kinner I) |
| 5 | best = CPO | `CPO_like_FamilyI` (Kinner I; Kinner puts LfuCPO in Family I) |

Confidence: `high` if metal-site type (H → I; G/A/T/S → II), length class and C-terminal Cys pair all agree with
the HMM family and the HMM margin is ≥ 20 bits; `medium` with one disagreement; `low` with more. The margin of 20
bits is [assumption].

### 4.5 Deviations from the paper

| # | Paper | Here | Why |
|---|---|---|---|
| D1 | 113 seeds (stage 2) | 35 of them (Table 3) + 13 characterized + 31 CPOs | The 113 are not published; UPObase is offline. |
| D2 | phmmer is the detection step | PF01328 GA + family HMMs decide; phmmer is only reported | Per-protein classification of BFD proteomes; BFD already runs PF01328. phmmer with all seeds was not tested at BFD scale. |
| D3 | cd-hit 90 %, MCL 1.4, cd-hit 95 % | Not done | These build a non-redundant database. MCL cluster selection is not described, so it cannot be replicated. Per-protein calls do not need it. |
| D4 | Motif search method not described | Fixed sites by alignment to AaeUPO (Kinner checks); raw-sequence ordered regex (UPObase Table 1) | Alignment avoids substring false hits; regex is the closest replication of Table 1. |
| D5 | Classes from phylogeny + motifs | Classes from HMM competition + motifs; no tree | A tree per protein does not scale. HMM competition is a proxy [assumption]; see 5.4. |
| D6 | Ensembl Fungi 812 genomes (2018) | Any input FASTA (BFD) | Different data. |
| D7 | Primary classes UPO I–V vs Pog | Primary classes Kinner I / II (user decision); UPObase class in its own columns | User decision. |

## 5. Validation [measured]

### 5.1 Characterized enzymes

Self-classification (`tests/refs_selfclass.tsv`) and leave-one-out (`tests/loo/loo_results.tsv`, models rebuilt
without the test protein; AaeUPO cannot be left out because it is the anchor).

| Enzyme | Expected (basis) | LOO result | Confidence |
|---|---|---|---|
| AaeUPO | II (Kinner text) | II (self only) | high |
| PabUPO-II, UPO1 | not stated in sources read | II | high |
| CraUPO | not stated (UniProt fragment, 261 aa) | II | low (no C-terminal region) |
| MroUPO, MweUPO, CviUPO | I (Kinner text; 7ZCL title "short") | I (UPO_FamilyI) | high, high, medium (margin 18.8) |
| TteUPO, DcaUPO, HspUPO, DchUPO, DspUPO | not stated in sources read | I (UPO_FamilyI) | high |
| LfuCPO | I (Kinner text) | CPO_like_FamilyI | medium |

12 of 12 left-out characterized enzymes get the expected Kinner family. Limits: MroUPO and MweUPO are 96 % identical,
so leaving one out still leaves a near copy. For 7 enzymes the expected family is my inference from length and
motif, not a statement in a paper I read.

### 5.2 All training references, leave-one-out

70 left-out references: the best family HMM equals the training family in 70 of 70. The final class matches in 62.
The 8 misses are CPO-group sequences that fail a motif check (acid-base column N/Q, or metal-site Q/A instead of E).
They fail the same way when they are in the training set, so this is a property of the sequences, not of the HMMs.

### 5.3 Negatives

- 13 of 13 negative reference proteins: no PF01328 hit and no family-HMM hit (`no_HTP_hit`).
- Four whole BFD proteomes (63,913 proteins; Jimgerdemannia flammicorona F268574B, Aureobasidium melanogenum
  F83D4B81, Cyclocybe aegerita F85C3D1B, Coprinopsis cinerea F1F5452B; `tests/proteomes/`): 43 candidates.
  42 are the PF01328 proteins in BFD `pfam.parquet` for these genomes (42 of 42 found). One extra protein
  (F85C3D1B_006062-T1, 160 aa, no PF01328 row) is a fragment with 67.5 % identity to a characterized UPO.
  None of the 479 P450 (PF00067), 32 class II peroxidase (PF00141), 14 catalase (PF00199) or 3 DyP (PF20628)
  proteins in these genomes was flagged.
- The Cyclocybe aegerita proteome has 13 UPO_FamilyII calls; the best match to AaeUPO is 85.8 % identity, not
  100 %. Reason not checked (strain or gene-model difference).

### 5.4 UPObase motif classes on the seeds

UPObase Table 1 regex on the 35 seeds: Ia 13, V 9, IV 5, II 3, IIIb 3, IIIa 1, Pog-g 1 (Mixia osmundae). The paper
places all 113 stage-1 sequences in subfamilies I–V, so 34 of 35 agree. AaeUPO → Ia, MroUPO → Pog-f,
LfuCPO → Pog-g; these agree with the paper (AaeUPO in subfamily I; MroUPO and LfuCPO in Pog). PabUPO-II (long,
EGD) → Pog-g: the regex misses it.

## 6. Results on the 807 early-diverging PF01328 proteins [measured]

Input: `../data/early_pf01328_hits.fa.zst`. Output: `results/early807_upo_id.tsv`, `.parquet`,
`results/early807_class_counts.tsv`. Run time 2.5 min on the login node.

| final_class | n | notes |
|---|--:|---|
| CPO_like_FamilyI | 665 | 569 high, 95 medium (70 of them because length ≥ 330 aa), 1 low |
| UPO_FamilyI | 3 | 2 Chytridiomycota (Boothiomyces, 237 aa), 1 Glomeromycota (Geosiphon pyriformis, 637 aa, possible fused gene model [assumption]); HMM margin 3–5 bits |
| UPO_FamilyII | 6 | all Mucoromycota, NCBI-annotated genomes (Jimgerdemannia flammicorona ×2 genomes, Endogone sp.); all **low**: EHD, 284–303 aa, no C-terminal Cys pair |
| HTP_fails_motif | 133 | 51 acid-base only, 35 Cys only, 22 Mg-Glu only, 25 mixed/partial |

Kinner family: 668 Family I (665 + 3) pass all motif checks. No protein has the full Family II signature (EGD +
long + C-terminal Cys pair). The 6 FamII-HMM calls score 64–75 bits higher with FamII than with CPO, and in a
FastTree placement with the references they form a clade at the base of the long-type group. Their motifs and
length are Family I. Treat them as an unresolved lineage, not as Family II UPOs.

UPObase classes on the same 807: Pog-g 693, Pog-b 24, Pog-f 3, none 87 (no PCP substring); 0 in UPO subfamilies
I–V. Catalytic label: PCP-EHD-E 692, neither 115, PCP-EGD-R-E 0.

Comparison with UPObase: UPObase reports 123 of 1,948 (6.3 %) in UPO subfamilies and 1,509 (77.5 %) in Pog-g.
The data sets differ (812 mostly Dikarya genomes from Ensembl in 2018, 95 % cluster representatives; here 807
proteins from early-diverging fungi, no dereplication). The numbers are not comparable as rates. The direction
agrees: Faiza 2019a found UPO subfamilies only in Basidiomycota and Ascomycota, and none appear here.

## 7. Open questions and limits

1. **CPO-like vs UPO Family I.** This split rests on 8 characterized short UPOs against 16 NCBI "CPO" sequences
   with no activity data. Kinner puts LfuCPO in Family I. I have no data that the CPO_like proteins lack
   peroxygenase activity. Read `CPO_like_FamilyI` as "Family I, nearer the CPO references".
2. **HMM competition as a phylogeny proxy.** Bit scores of HMMs of different length and depth are not strictly
   comparable. The FamII HMM (384 states, 47 sequences) can favor proteins that the FamI_UPO HMM (227 states,
   8 sequences) fits less well. A phylogenetic placement step (e.g. EPA-ng on a fixed tree) would be better.
3. **Ascomycete long-type acid-base pair.** No structure was available. The alignment does not put their R-X-E in
   the AaeUPO columns. The window rule (b2) accepts R-X(6,7)-E anywhere in an 83-residue window, so it can accept a
   chance match.
4. **Seeds.** Only 35 of the 113 stage-1 sequences and none of the 1,948 stage-2 sequences are available. Asking
   the authors (Faiza, Wang; yonghw@scut.edu.cn) for the UPObase FASTA is still the best route to a full replication.
5. **CglUPO, CciUPO, MthUPO, HinUPO** accessions remain unverified. Seed EAU83038 (Coprinopsis cinerea, 352 aa) may
   be CciUPO; not verified.
6. **Signal peptide** is read from the BFD SignalP table as given (positives only per BFD docs; not tested).
7. The 330-aa length threshold and the 20-bit margin are calibrated on 79 references only.
