# Design: fungal peroxygenase diversity, substrate prediction, and regulation

Written 2026-11-07. This document supersedes the project framing in `PLAN.md`, which was written for the
earlier question ("are there new UPOs in early-diverging fungi?"). That question has been answered and the
answer was largely negative; see section 1.2. `PLAN.md` remains valid for its census, data paths and the
decisions recorded at its top.

Companion documents, all with measured numbers:
- `CENSUS.md` — genome and PF01328 counts per phylum.
- `SUBSTRATE_CALIBRATION.md` — the published enzyme/substrate data this project calibrates against.
- `POCKET_TRANSFER.md` — the pilot that tests whether the key structural feature is computable at scale.
- `upo_id/METHOD.md` — the existing classification tool.

---

## 1. Intent

### 1.1 The question

Three questions, asked of the whole BFD dataset rather than one clade:

1. **Diversity and novelty.** What peroxygenase subtypes exist across fungi, how are they distributed, and
   is any subtype lineage genuinely new relative to what is already public?
2. **Function.** Which subtypes are worth testing experimentally, judged from sequence, structure and
   genomic organisation?
3. **Regulation.** Do peroxygenase genes share regulatory features that would tell us how fungi control
   them?

### 1.2 Why the framing changed

The original project asked whether BFD's re-annotation of genomes NCBI left unannotated reveals new
unspecific peroxygenases in Mucoromycota, Glomeromycota and Chytridiomycota. It ran, and the answer was no:
of 807 early-diverging PF01328 proteins, the 9 called as UPO all sat in genomes NCBI had already annotated.

A review of that result found it was also **stated wrongly**. Under the user's own decision to treat UPO
Family I and CPO-like as one Family I, the correct statement is that 668 proteins pass all motif checks as
Family I, and 273 of those are in BFD-only genomes. (Verified 2026-11-07 from
`upo_id/results/early807_class_counts.tsv`: `CPO_like_FamilyI` rows sum to 665, of which the `BFD_only` rows
sum to 273; adding the 3 `UPO_FamilyI` rows, all NCBI, gives 668 Family I. This is consistent with
`CENSUS.md` §3's 309 BFD-only of all 807, which counts every PF01328 hit including the 133 that fail motif
checks.) The "no new UPO" headline came from a narrow class
definition the project had already abandoned. Two further problems: the 9 calls had HMM margins of 3–5 bits
against the tool's own 20-bit rule, and the "none in BFD-only" finding is partly a sampling artifact, because
the relevant lineage has exactly one BFD-only genome and it carries no PF01328 protein at all.

So the negative result does not close the subject. It shows the original question was too narrow, and that
the interesting variation is **within** Family I rather than between Family I and Family II.

### 1.3 Why this is worth doing

Peroxygenases oxidise C–H bonds using only hydrogen peroxide, which makes them attractive where a
cytochrome P450 would need a reductase, a cofactor and a membrane. The published enzyme set is small —
on the order of a dozen well characterised proteins — against roughly 97,520 PF01328 proteins in BFD and
thousands of putative genes noted in the literature. The gap between what exists and what has been
characterised is the opportunity, and BFD's uniform re-annotation across ~23,000 genomes is a view of that
space that has not been available before.

---

## 2. Goals, non-goals, anti-goals

### Goals

- **G1.** A classified, lineage-aware catalogue of PF01328 proteins across BFD, with subtype assignments and
  a measure of how novel each subtype lineage is against public sequence databases.
- **G2.** A ranked candidate list for heterologous expression, scored on features that are measurable from
  sequence, structure and genomic context, with the ranking's predictive content tested against published
  enzyme data before it is applied.
- **G3.** A first description of peroxygenase gene regulation in fungi: shared upstream features, and what
  the genes co-express with where RNA-seq allows it.
- **G4.** Everything reproducible from the repository, with measured and assumed numbers labelled
  differently. This is an existing project norm and it holds here.

### Non-goals

- Predicting turnover rate, thermostability or expression yield. Not possible from sequence and structure.
- Characterising enzymes outside the heme-thiolate peroxygenase family.
- A complete fungal phylogeny. The species tree is an input, not a deliverable.
- Wet-lab work inside this project. Expression and assay are expected within about a year and are a
  downstream consumer of G2's output, not part of it.

### Anti-goals — these would make the project a failure even if it ran

- **A ranking that looks confident and is not.** With roughly 19 calibration enzymes and 4 measured
  negatives, a many-parameter model will fit noise and report it as insight. Any scoring scheme must stay
  small enough that its parameters are all traceable to a published observation.
- **Claiming novelty that is only an annotation gap.** "Absent from NCBI's annotation of that assembly" is
  not "absent from public databases". The novelty measure must be against sequence databases, not against
  `prior_annotation_status.tsv`.
- **Mistaking phylogeny for biology.** Both the regulatory analysis and the subtype-novelty measure will
  produce strong-looking signals from shared ancestry if related genomes are counted as independent
  observations.
- **A catalogue nobody can act on.** G2 must end in something an expression partner can pick from, with the
  reason for each candidate attached.

---

## 3. Constraints

- **Compute.** UCR HPCC. `-p epyc,intel --constraint="ryzen|broadwell|cascade|sapphire"`; never `abu_dhabi`
  nodes, which lack AVX2. Work in `$SCRATCH`, copy results to `/bigdata`. Size jobs toward 1–1.5 h each.
- **Storage.** Compress tabular and sequence output; `zstd` for internal intermediates, `.gz` where a tool
  demands it. HMM files stay plain because `hmmsearch` and `hmmalign` cannot read compressed input.
- **Interpreters.** `/usr/bin/python3.12`. The default `python3` on HPCC is 3.9 and will fail.
- **No `BASH_SOURCE`** in anything that runs under SLURM.
- **Scale.** PF01328 covers about 97,520 proteins BFD-wide. Any step that is more expensive than an
  alignment must run on a dereplicated or representative subset, not the full set.
- **Review.** Earlier attempts to get this work reviewed were stopped three times by model safeguards. A
  review did succeed on 2026-11-07. Reviews should state plainly that this is computational annotation of
  industrial biocatalysts.

---

## 4. What is already measured

Facts, with sources. Nothing here is an estimate.

**Census** (`CENSUS.md`): 827 early-diverging genomes, 813 with proteins, 10.17 M proteins. 807 PF01328
proteins in early-diverging genomes, 309 in BFD-only genomes. Zero PF01328 at the gathering threshold in
Mortierellomycota, Kickxellomycota, Blastocladiomycota, Neocallimastigomycota, Entomophthoromycota and
Microsporidia.

**Copy number** (measured 2026-11-07 from `pfam.parquet`; counts are PF01328 *rows*, so an upper bound on
distinct proteins):

| Group | Genomes | With ≥1 hit | Median/genome | Max |
|---|---|---|---|---|
| Eurotiales | 2,489 | 2,402 | 4 | 12 |
| *Aspergillus* | 1,598 | 1,593 | 5 | 12 |
| *Penicillium* | 646 | 609 | 2 | 8 |
| *Monascus* | 25 | **0** | 0 | 0 |
| Agaricomycetes | 1,856 | 1,790 | 7 | 147 |
| Tremellomycetes | 559 | 190 | 0 | 8 |
| Malasseziomycetes | 85 | **1** | 0 | 2 |

The clean zeros are as informative as the expansions. *Monascus*, *Malassezia* and two thirds of
Tremellomycetes have none, which looks like real loss in yeast-like and lipid-dependent lineages. One
Agaricomycete has 147 copies; that is either a large expansion or an assembly artifact and must be checked
before it is quoted anywhere.

**Calibration data** (`SUBSTRATE_CALIBRATION.md`): a 2025 study assayed 19 UPOs across 5 clades in one host
under identical conditions, including naphthalene, with 4 measured negatives. It reports that a narrow
access channel plus a conserved "aliphatic pocket" (HspUPO residues I73, L77, L81, A183, I187, L228) tracks
peroxygenase rather than peroxidase activity. Separately, Karich et al. 2017 give a 14-substrate PAH panel
for two enzymes with substrate dimensions, and Linde et al. 2020 show a single channel residue (Phe88)
moving aromatic activity 5–21 fold.

**Pocket transferability** (`POCKET_TRANSFER.md`, measured here): across 87 PF01328 proteins spanning 29
lineages at ~25 % identity to HspUPO, the six pocket columns are alignable **95.6 %** of the time against a
catalytic-site ceiling of 97.1 %, but only **68 %** of mapped positions are aliphatic. Length class does not
explain the variation. The feature is computable everywhere; it is graded, not binary.

### 4.1 Measured 2026-11-07: the aliphatic pocket does NOT predict aromatic ring activity

Found by an independent builder review of this design and then verified directly. This falsifies the first
version of decision 5.2 #7 and it is the most important finding in this document.

The four UPOs with no detectable naphthalene activity in the 2025 study are **AbrUPO-II, CviUPO, PanUPO and
TteUPO** (verified verbatim from the paper, [DOI](https://doi.org/10.1186/s13068-025-02675-w); the enzyme
names were stripped from the text first retrieved and were confirmed against the published page).

Cross-referencing those calls against the pocket scores measured in `POCKET_TRANSFER.md` §3.1:

| Enzyme | Clade | Pocket /6 | Naphthalene | Pocket residues |
|---|---|---|---|---|
| AaeUPO / PaDa-I | I | 3 | **active** | FAAFL– |
| MroUPO | III | 5 | **active** | VIAFVL |
| HspUPO | V | 6 | **active** | ILLAIL |
| DspUPO-I | V | 6 | **active** | ILLAIL |
| CviUPO | II | 4 | **none** | LLIKTM |
| TteUPO | V | **6** | **none** | VLLAAL |

**Mean pocket score is 5.00 for naphthalene-active enzymes and 5.00 for naphthalene-inactive ones.** Zero
discriminating power. TteUPO carries the full six-residue pocket and has no naphthalene activity at all.

This is mechanistically coherent rather than surprising. The pocket was defined as predicting peroxygenase
activity on **4-propylguaiacol**, a single phenolic ring with an aliphatic propyl side chain, where the
reaction is C–H abstraction on that side chain or on the methoxyl. An aliphatic pocket should favour holding
an aliphatic chain near the heme. Hydroxylating a naphthalene **ring** is different chemistry — epoxidation
of an aromatic system. There is no reason the same feature should predict both, and it does not.

**Consequence:** the aliphatic pocket is a feature for aliphatic-substrate peroxygenation. This project
targets aromatic ring hydroxylation. The pocket is therefore a useful covariate and a possible negative
marker, but it is **not** the predictor for the project's goal. T1 and T2 must be scored separately.

**Weak lead, not a finding.** Counting aromatic residues (F/W/Y) at the same six positions gives 0.75 for
the naphthalene-active enzymes and 0.00 for the inactive ones, and AaeUPO — the strongest aromatic
hydroxylator in the literature — has the most. The direction matches the stacking argument and matches
Linde et al.'s Phe88 result. **n = 6 and this is not evidence.** It is the hypothesis §7.3 must test once
the full calibration matrix exists.

Caveats on the table: AaeUPO here is wild-type B9W4V6 while the paper's Clade I entry is PaDa-I, the
laboratory-evolved variant; and our DspUPO (PDB 9WA6) is assumed to be the paper's DspUPO-I but this is not
confirmed. Both should be checked before the table is used for anything beyond the conclusion above, which
rests on TteUPO and does not depend on either.

---

## 5. Decisions

### 5.1 Made by the user

1. Scope is the **whole BFD dataset**, not only early-diverging fungi.
2. Substrate focus is **aromatic hydroxylation of carbon rings, one to three rings**. Larger PAHs are out of
   scope as a target but may be reported as extrapolation.
3. The deliverable for G2 is a ranked candidate list, not a catalogue alone.
4. Regulation (G3) starts from **sequence**, then brings in *Penicillium* / *Aspergillus* RNA-seq.
5. Expression will be in **Pichia**, of a Eurotiales or Basidiomycete copy.
6. Classification keeps the earlier decisions recorded at the top of `PLAN.md`: Kinner Family I / Family II
   scheme; report all classes inclusively; build our own reference set rather than depending on UPObase;
   treat UPO Family I and CPO-like as one Family I for now; replace HMM score competition with phylogenetic
   placement.

### 5.2 Made by me, and open to reversal

7. **CORRECTED 2026-11-07 — see §4.1. The model has TWO separate target variables, scored and validated
   independently, never combined into one score:**
   - **T1, peroxygenase-versus-peroxidase propensity.** Predicted by the aliphatic pocket and channel shape.
     Supported by the 2025 study.
   - **T2, aromatic ring hydroxylation (the project's actual target).** **Not** predicted by the aliphatic
     pocket. Features for T2 are an open problem; the best lead is aromatic channel lining (§7.1).

   Why not ring-size band: the 1–3 ring substrates all have minimum widths of 3.0–4.0 Å and clear both the
   narrow (7 Å) and wide (11 Å) reference channels. Channel width separates 4–5 ring PAHs from everything
   else and is flat inside the chosen range. A width-keyed model would fit noise.

   The original version of this decision folded T1 and T2 into one score. §4.1 shows that was wrong.
8. **The aliphatic pocket is scored 0–6, graded.** Forced by the measurement in `POCKET_TRANSFER.md`.
9. **HspUPO is the numbering anchor** for pocket and channel positions, recovered by matching the published
   residue numbers against the reference set.
10. **Scoring is a small, transparent rule set with every parameter traceable to a publication — not a fitted
    classifier.** Forced by n ≈ 19 with 4 negatives.
11. **Compute is tiered**: sequence features on everything, structure on clade representatives, docking on a
    shortlist.

### 5.3 Decisions still needed from the user

- **D1. The Eurotiales tension.** Eurotiomycetes scored 2.3/6 on the aliphatic pocket, among the lowest of 29
  lineages, while Ustilaginomycetes scored 6.0/6, the highest. If the pocket is the feature that matters, the
  best candidates are not in *Aspergillus* or *Penicillium*. This does not affect *Pichia* as the host, only
  the source organism. Is the Eurotiales preference firm, or was it a convenience?
- **D2. Scope of the regulatory work.** Sequence-only motif analysis across all of BFD, or restricted to the
  lineages where RNA-seq can follow up?
- **D3. Whether subtype novelty is measured per protein or per clade.** Per clade is the more meaningful
  claim and the more expensive one.
- **D4. RESOLVED (user, 2026-11-07): the 2025 study's Clades I–V define a subtype.** The user's words:
  "2025 study's Clades I–V are better explanation". This is the scheme the aliphatic pocket, the substrate
  panel and all of the §4.1 calibration data are defined on. Kinner Family I/II is retained as a **separate
  column**, not as the primary partition, so user decision 6 (merge CPO-like into Family I) still applies to
  that column only. D3's novelty unit and G1's tree both key on Clades I–V.
  Verified clade memberships from the paper ([DOI 10.1186/s13068-025-02675-w](https://doi.org/10.1186/s13068-025-02675-w),
  PMC12306115): Clade I — CabUPO-III, CmaUPO-I, GmaUPO-II, PaDa-I, LspUPO-II; Clade II — CviUPO, PanUPO;
  Clade III — MorUPO, MroUPO; Clade IV — MphUPO; Clade V — AbrUPO-I, AbrUPO-II, AluUPO, AtuUPO, AacUPO,
  TteUPO, HspUPO, RneUPO, DspUPO-I.
  *Open sub-item for the builder:* BFD proteins must be assigned to Clades I–V by phylogenetic placement
  against these 19 plus the wider reference set. Only 19 sequences define the clades, so placement
  confidence must be reported per protein and low-confidence proteins must stay unassigned rather than be
  forced into the nearest clade.
- **D5. RESOLVED (user, 2026-11-07): AlphaFold3 is acceptable; try Boltz-2 as well.** The user's words:
  "AlphaFold3 if you can run it is okay, we are all non-commercial here but if Boltz-2 is installed I am
  happy to see it tried." The non-commercial licence is therefore not a blocker. **Both are installed and
  verified on this cluster (2026-11-07):** `alphafold/3.0.0`, `alphafold/3.0.2`, `boltz/2.2.1`
  (`/opt/linux/rocky/8.x/x86_64/pkgs/boltz/2.2.1/env/bin/boltz`), and `boltzdesign1/0.4.1`.
  **Builder instruction:** run both on the same small set of clade representatives first and compare, rather
  than committing to one. The comparison must be made on something measurable, not on preference — propose
  per-model confidence at the heme-contacting residues, and agreement between the two models on the channel
  geometry that Tier 2 actually reads. If the two models disagree on channel geometry for the same sequence,
  that disagreement is a result and caps how far any structure-derived feature can be trusted. Caveat: both
  are structure predictors validated on folds, not on cavity geometry; neither has been benchmarked here for
  channel measurement. GPU budget is still unquantified and must be measured on the pilot set before any
  scale-up.
- **D6. G3 scope (was D2), restated with its cost.** All-BFD means extracting and background-matching
  promoters for ~97,000 genes across ~23,000 genomes. Restricting to RNA-seq-followable clades means a few
  hundred Eurotiales promoters. These are different projects, not different settings.

---

## 6. Workstream G1 — diversity and novelty

**Input:** all PF01328 proteins in BFD, with lineage from `samples.parquet`.

**Steps.**

1. **The `upo_id` training-set fix is now largely moot — do not treat it as a blocker.** The original
   finding stands: the 16-sequence CPO training set contains a Mucoromycota sequence (*Absidia repens*) and
   a Lecanoromycetes sequence (*Umbilicaria pustulata*), neither with activity data, so a query from those
   lineages scores better against a profile containing its own lineage. **But** the contamination is in the
   CPO profile, and user decisions 6 and 7 (merge CPO-like into Family I; replace HMM score competition with
   phylogenetic placement) delete the CPO-versus-FamI_UPO comparison that the bias acts on. Under the
   current design the bias has nothing to distort.
   Action: remove the two sequences, re-run the 807-protein set once, and record how far labels move, as a
   measurement of how much the earlier results were affected. Then retire the HMM-competition path. This is
   a half-day of bookkeeping, not a gate.
2. **Dereplicate.** 97,520 proteins is too many for one alignment or tree. Cluster, keep representatives,
   retain the mapping back to every member.
3. **Subtype tree.** Build a reference tree from characterized enzymes plus a broad sample, with an explicit
   outgroup and bootstrap support, then place BFD representatives on it. **Write the decision rule before
   coding**: what placement result counts as which call. Without that rule, placement yields "sister to X"
   and no classification. Note the 2025 study's observation that enzymes in one clade can differ in
   regioselectivity while two at 96 % identity agree — so clade membership alone will not determine function,
   and residue-level features are needed as well.
4. **Novelty, measured properly.** DIAMOND the representatives against nr and UniProt (`nr.dmnd` and
   `swissprot.dmnd` under `$DIAMOND_DB`, `module load db-diamond/20250623`); report best public identity per
   protein and per clade at two identity cutoffs, not one. A clade is novel when it has no close public
   relative. Do **not** use `prior_annotation_status.tsv` as the novelty measure — that records whether NCBI
   annotated an assembly, which is a different question.

   **Known residual confound, must be reported with every novelty call.** A BFD-only protein is absent from
   nr partly *because* nobody annotated that assembly. So protein-database novelty still partly measures
   deposition rather than biology — a weaker form of the anti-goal in §2, not an escape from it. Mitigate by
   also searching nucleotide and WGS space, where an unannotated genome's sequence may still be present, and
   report public-protein status alongside the identity figure so a reader can separate the two.
5. **Lineage dynamics.** Copy number per genome per clade mapped onto the species tree: expansions, losses,
   single conserved copies. The measured zeros in *Monascus*, *Malassezia* and Tremellomycetes are a starting
   point.

**Controls.** Contamination screening is required before any lichen or Glomeromycota result is reported:
neighbour-gene best-hit taxonomy on the same contig, GC and codon usage against the host genome, and a flag
for contigs carrying only the candidate. Lichen assemblies may contain photobiont or bacterial sequence; this
is untested. Genome-quality confounding must also be reported — several zero-hit phyla have poor assemblies
(Kickxellomycota: 123 of 187 genomes below 80 % BUSCO), so their zeros are not yet interpretable as absence.

---

## 7. Workstream G2 — substrate and activity scoring

### 7.1 What is being predicted

**Two independent outputs per protein. They are never summed into one number.** See §4.1 for why.

**T1 — peroxygenase-versus-peroxidase propensity.** Predicted by aliphatic pocket score (0–6) and channel
shape (constrained versus funnel). Basis: the 2025 study. This axis is supported and ready to build.

**T2 — aromatic ring hydroxylation.** This is the project's actual target and **there is currently no
validated predictor for it.** The aliphatic pocket has been measured to have zero discriminating power
(§4.1). Candidate features, in order of support:

1. **Aromatic residues lining the access channel.** Rationale: a ring substrate must stack against the
   channel wall to be held above the heme. Support: Linde et al. 2020, where a single phenylalanine
   substitution (F88L in CviUPO) changed efficiency on naphthalene and other simple aromatics 5- to 21-fold
   ([DOI](https://doi.org/10.1128/AEM.02899-19)); plus the weak n=6 signal in §4.1. **This is a hypothesis.**
2. Channel minimum cross-section, as a permissive filter only — it does not discriminate inside 1–3 rings.
3. Electronic environment near the heme. Unspecified; would need the primary structural literature.

**If T2 has no predictor after §7.3, that is a reportable result and the ranking falls back to T1 plus
distance from characterized space.** Do not ship a T2 score that failed its gate.

Substrate band is reported only as 1-ring / 2–3 ring, with anything larger flagged as extrapolation.

### 7.2 Tiered computation

**Tier 1 — sequence, runs on everything.** Read out at HspUPO-anchored positions: the six aliphatic-pocket
positions (score 0–6); channel-lining residues, with Phe88 (Linde et al.) and Phe302/Pro277 (2025 study) as
worked examples; the catalytic Cys and the EHD metal site as validity checks.

Two implementation constraints, both from the builder review:

- **MAFFT L-INS-i does not scale and must not be used here.** The pilot used it on 97 sequences; it is
  quadratic in pairwise alignments and will not run on ~97,520. Use `hmmalign` against an HspUPO-anchored
  profile built with `--hand`.
- **The 95.6 % alignability figure in `POCKET_TRANSFER.md` was measured under L-INS-i, not under
  `hmmalign`, so it does not transfer.** Re-measure it on the same 87-protein sample under the profile
  method that will actually run, before trusting Tier 1 at scale. If `hmmalign` is materially worse, that
  changes the design.
- **Anchor mismatch to resolve:** `upo_id/METHOD.md` §4.3b is anchored on AaeUPO precursor numbering; this
  design is anchored on HspUPO. Write and test an explicit position map between the two, or re-anchor the
  existing tool. Do not maintain two unmapped coordinate systems.

**Tier 2 — structure, clade representatives only.** Model with the heme in place. Standard AlphaFold2 does
not place cofactors, so this needs an AlphaFold3-class model or template-based placement from 2YOR and 5FUJ.
Then measure the access channel: shape (constrained versus funnel with a narrow base), minimum
cross-section, and lining hydrophobicity. The 2025 study used this combination, which is evidence the
approach works in this family — not evidence that our implementation will reproduce it.

Two practical constraints:

- **CAVER is not installed on HPCC**, and Caver Analyst 2.0 is a GUI desktop application. Batch use needs
  the CAVER 3.0 jar installed and scripted, or a substitute (MOLE, or a custom probe-based measure). Budget
  this as real setup work, not a module load.
- **Model choice: settled by D5 — run both and compare.** Licensing is not a blocker; the user confirmed
  the work is non-commercial. Installed and verified 2026-11-07: `alphafold/3.0.0`, `alphafold/3.0.2`,
  `boltz/2.2.1`, `boltzdesign1/0.4.1`. Run AF3 and Boltz-2 (heme as an explicit ligand) on the same clade
  representatives and compare before committing. Neither tool has been benchmarked here for cavity
  geometry, which is the only thing Tier 2 reads from the model, so the agreement between them is the
  available estimate of how far that geometry can be trusted.

**Tier 3 — docking, shortlist only.** Dock a small substrate ladder, anchored on the 1–3 ring range, plus a
fatty acid as the contrasting case. Filter solutions on the oxidised carbon sitting 3.5–4 Å from the heme,
following the published criterion.

### 7.3 Validation, which comes before any large run

Run Tiers 1–3 on the characterized enzymes first. **Two gates, one per target variable.**

**Gate T1** — does the scheme reproduce the 2025 study's separation of peroxygenase-leaning from
peroxidase-leaning enzymes? Note this gate is **partly circular**: the pocket residues were defined on
Clade V, and the enzymes that score highest are the 59–73 % identity Clade V relatives. Passing it is close
to guaranteed and proves little. Treat it as a smoke test that the pipeline runs, not as validation.

**Gate T2** — does any candidate feature separate the four measured naphthalene negatives (AbrUPO-II,
CviUPO, PanUPO, TteUPO) from the actives, across all 19 enzymes? This is the real gate. It is not circular,
because the features were not defined on this outcome. **The aliphatic pocket has already failed it**
(§4.1). Test aromatic channel lining next.

**Do not use the Karich MroUPO/AaeUPO bulky-PAH contrast as a gate.** It is a 4–5 ring observation, which
§5.2 #7 puts out of scope and which `SUBSTRATE_CALIBRATION.md` §11 shows is uninformative in the 1–3 ring
range. Keep it as an out-of-scope sanity check only.

If T2 fails for every candidate feature, drop the T2 claim and report the negative. That outcome is cheap,
honest and more useful than an unvalidated ranking.

**Statistical honesty requirement.** With 19 enzymes, 4 negatives, and two members at 96 % identity, the
effective sample is smaller than 19. Any accuracy figure computed on this set is **in-sample**. Report it as
a sanity check, never as validation, and never quote a cross-validated-sounding number. Free parameters —
the pocket cutoff, the channel-class boundary, how features combine — must be counted and listed explicitly
in the results, because "a rule set with traceable parameters" is still a fitted model if those are tuned.

**A conflict to resolve first:** the 2025 study found no naphthalene activity for CviUPO, while Linde et al.
2020 report a catalytic efficiency for exactly that reaction. Different hosts and assay conditions are the
likely cause. CviUPO cannot be used as a positive or a negative until this is settled.

### 7.4 Genomic organisation

With a substrate class fixed, gene neighbourhood becomes interpretable. Look for copies near aromatic
catabolism machinery — ring-cleaving dioxygenases, quinone reductases, aryl-alcohol oxidases, laccases — and
for membership in biosynthetic gene clusters, which would supply an endogenous substrate hypothesis. Also
record tandem array versus dispersed, transposon proximity, and synteny conservation, which separate ancient
orthologs from recent duplicates. Note that *Moesziomyces aphidis* carries two UPOs with opposite preferences
(one PAH, one n-alkane), so expansions may partition substrate space rather than duplicate a function.

### 7.5 Expressibility filter

Since the endpoint is *Pichia* expression, carry per-candidate flags: N-glycosylation site count, cysteine
count and pairing, length, and whether the native signal peptide would be replaced. These do not predict
yield. They flag copies likely to give trouble, and they are a negative filter only.

### 7.6 Output

One row per candidate: lineage, clade, annotation group, pocket score, channel class, predicted
peroxygenase:peroxidase propensity, predicted naphthalene activity, genomic context, nearest characterized
enzyme and identity, expressibility flags, QC flags. Ranked by distance from characterized space **with the
reason attached**, so a reader sees why each candidate is unusual rather than a bare score.

---

## 8. Workstream G3 — regulation

This is the most open of the three and the literature supports it as a genuine gap: the 2017 review states
plainly that it is unclear under what conditions UPOs are induced and expressed in fungi.

**Phase A — upstream sequence.** Extract a fixed window upstream of each clean-gene-model PF01328 copy,
**truncated at the neighbouring gene**. Fungal genomes are compact and a fixed 1 kb will often run into the
next coding sequence; the truncation is not optional. Then two analyses in parallel:

- *A priori scan* for known fungal regulator motifs, counted per clade against background genes from the same
  genomes. The strong prior is iron and heme, since this is a heme-thiolate enzyme, so GATA-type sites for
  SreA and HapX are a specific testable hypothesis. Alongside: CreA for carbon catabolite repression, PacC
  for pH, AP-1 type sites for oxidative stress.
- *De novo discovery* against a matched background. Expect noise; fungal promoter discovery usually returns
  low-information motifs.

**The control that decides whether Phase A means anything.** Copies in related genomes share upstream
sequence by descent, not by shared regulation. Without correction, any enrichment measures phylogeny. Take
one copy per species per clade, or test within a clade against same-genome background. Collapse tandem
arrays, which share flanking sequence.

**Phase B — expression.** BFD's own RNA-seq list covers 52 species including four *Aspergillus*/*Penicillium*
and three *Talaromyces*, but it is gene-prediction support — roughly one condition per species — so it can
confirm gene models and cannot give co-expression. Co-expression needs many conditions in one species, which
means public SRA for a well-sampled host such as *A. niger*, *A. nidulans*, *A. fumigatus* or *P. rubens*.
Measured copy numbers support this: Eurotiales genomes carry about 4 copies each.

The outcome that would matter: if UPO copies co-express with iron and heme genes, Phase A's prior is
confirmed. If they co-express with a biosynthetic cluster, that is a substrate hypothesis feeding G2.

---

## 9. Risks

| Risk | Status | Mitigation |
|---|---|---|
| **No validated predictor exists for aromatic ring activity** | **CONFIRMED — pocket failed, §4.1** | Test aromatic channel lining at Gate T2; if it fails, drop the T2 claim and report it |
| Calibration matrix (S1/S2/S4) not yet retrieved | **Confirmed blocker for G2** | Retrieve before writing any scoring rule; without it n drops to ~10 and Gate T2 is not a test |
| Gate T1 is partly circular | **Confirmed** | Treat as a smoke test, not validation (§7.3) |
| Tier 1 alignability measured under a method that will not scale | **Confirmed** | Re-measure the 87-protein sample under `hmmalign` before scaling (§7.2) |
| Two unmapped coordinate systems (AaeUPO vs HspUPO) | **Confirmed** | Write and test an explicit position map (§7.2) |
| Pocket score is a sequence proxy of unknown structural fidelity | **Open** | Model high and low scorers, check channels (§7.2) |
| Calibration set too small for a real model | **Confirmed**, n≈19, 4 negatives, 2 at 96 % id | Count and publish free parameters; report in-sample only (§7.3) |
| Novelty-vs-nr inherits a deposition confound | **Confirmed** | Search nucleotide/WGS too; report public-protein status alongside (§6 step 4) |
| CAVER not installed | **Confirmed** | Budget install, or substitute MOLE / custom probe measure |
| Structure tool choice | **Resolved (D5)** | AF3 and Boltz-2 both installed; run both, compare on channel geometry |
| Lichen and Glomeromycota contamination | **Untested** | Contig-level checks before any result is reported |
| Phylogenetic non-independence inflating regulatory signal | **Expected** | One copy per species per clade; same-genome background |
| Zero-hit phyla confounded by assembly quality | **Confirmed** | Report BUSCO alongside every absence claim |
| G3 Phase B unscoped; may exceed G1+G2 combined | **Confirmed** | Scope SRA curation explicitly before committing; copy number is not evidence co-expression is measurable |
| The 147-copy Agaricomycete | **Unverified**, and `CENSUS.md` §3 gives a Basidiomycota max of 144 | Reconcile the two queries before quoting either |

---

## 10. Sequencing of the work

1. ~~**Retrieve the calibration matrix**~~ — **DONE 2026-11-07. See `CALIBRATION_MATRIX.md`.**
   There is no Table S4; all quantitative data is in figures. The supplement was retrieved, Figure S2 (the
   active-site matrix) is transcribed to `data/calib_dl/active_site_S2.tsv`, and the CviUPO conflict is
   resolved (different assay sensitivity, not a contradiction). **Outcome: Gate T2 cannot be passed on this
   dataset.** Pocket-vs-naphthalene Pearson r = 0.18 across all 19 enzymes, and the only two interpretable
   negatives (CviUPO, PanUPO) are the complete membership of Clade II, so clade and label are perfectly
   confounded. Read `CALIBRATION_MATRIX.md` §5 before planning any G2 work.
2. **G2 RESOLVED (user, 2026-11-07): pivot to target T1, and report the T2 negative.**
   T2 (aromatic ring hydroxylation) is not modellable on the available calibration data — see
   `CALIBRATION_MATRIX.md` §4–5. It is written up as a measured negative result, not worked around.
   G2's modelling effort moves to **T1, peroxygenase-versus-peroxidase propensity**, which is the axis the
   aliphatic pocket was actually defined on. Figure S8 of the 2025 supplement holds per-enzyme T1 data at
   three pH values and still needs digitizing; that digitization is the gate on whether T1 is viable.
   **D6 (G3 scope) is the one decision still open.**
3. Re-measure Tier 1 alignability under `hmmalign` on the existing 87-protein sample; write the
   AaeUPO↔HspUPO position map. Clean up the `upo_id` training set as bookkeeping.
4. **Gate T2 on the full 19-enzyme set: does aromatic channel lining separate the naphthalene negatives?**
   This is the decision point for the whole of G2. A failure here is a result, and it should be reported
   rather than worked around.
5. Dereplicate; Tier 1 features across BFD; subtype tree with its decision rule written first.
6. Novelty by DIAMOND against protein and nucleotide databases.
7. Tier 2 structures and channel measurement on representatives; Tier 3 docking on the shortlist —
   **only if step 4 passed.**
8. G3 Phase A, which can run in parallel from step 5 onward.
9. G3 Phase B once D6 is settled and a host species is chosen.

---

## 11. Left to the builder

- Clustering identity for dereplication, and whether to reuse the DeltaGain linclust tiers.
- Tree method and the specific outgroup sequences.
- Motif discovery tool and background construction for G3 Phase A.
- Job sizing, chunking and the specific tool versions.
- Whether `upo_id` is extended in place or the scoring becomes a separate tool alongside it.
