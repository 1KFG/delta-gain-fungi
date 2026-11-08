# UPO project plan

**Rewritten 2026-11-07.** Supersedes the earlier sketch, kept verbatim as `PLAN.superseded.md`.
That sketch was written before the substrate-calibration work and before either activity gate was
run. Its user decisions 1-10 are still honoured except where measurement has overtaken them; the
mapping is in section 2.

Companion documents, all in this directory:

| File | What it holds |
|---|---|
| `DESIGN.md` | Intent, goals, constraints, decisions D1-D6, workstream definitions |
| `CALIBRATION_MATRIX.md` | The 2025 calibration set: retrieval, validation, both gate results |
| `POCKET_TRANSFER.md` | Pilot: does the aliphatic pocket transfer across fungi |
| `SUBSTRATE_CALIBRATION.md` | Literature survey on UPO substrate range |
| `CENSUS.md`, `REFERENCES.md` | Measured BFD counts; reference-set checks |
| `RESULT_T2_negative.md` | The G2 negative result, written up standalone |
| `upo_sites/` | The software built from all of the above |

---

## 1. Where the project actually stands

Three load-bearing questions have been answered with measurements. Two came back negative. This
plan is built on those answers rather than around them.

**Q1. Can the published active-site positions be found in divergent fungal proteins? Yes.**
`hmmalign` against PF01328 places five of the eight aliphatic-pocket positions, and at those five
it reproduces the paper's structure-based matrix for all 19 calibration enzymes at **95/95 =
100 %**. Across 87 phylogenetically spread BFD proteins the pocket columns align at 94.7 %,
against a catalytic-site ceiling of 95.7 %.

**Q2. Does the aliphatic pocket predict aromatic ring hydroxylation (target T2)? No.**
Pearson r = **+0.176** against naphthalene activity across all 19 calibration enzymes. Only two of
the four measured negatives are interpretable, and those two - CviUPO and PanUPO - are the complete
membership of Clade II, so label and clade are perfectly confounded. **T2 is not modellable on the
available data and will be reported as a negative result.**

**Q3. Does it predict peroxygenase-versus-peroxidase propensity (target T1)? Partly, but not
usefully.** Raw correlation is r = +0.748. But "is Clade V" alone gives r = +0.835, and the
**partial correlation of pocket with T1 given clade is only r = +0.129**. Within the ten
non-Clade-V enzymes it is +0.001. The pocket cannot distinguish GmaUPO-II (Clade I, pocket 3/5,
T1 = 53) from CabUPO-III (Clade I, pocket 3/5, T1 = 0).

**The conclusion that drives this plan: clade carries the measured signal, not the pocket.** With
n = 19 across 5 clades, "this residue matters" and "this clade behaves this way" cannot be
separated. The product is therefore a reliable **clade annotator with transparent site evidence**,
not an activity predictor.

### 1.1 What is deliberately not being built

A per-protein activity score. It would be easy to emit across 97,520 proteins and would rank every
fungal UPO by features measured at r = 0.18 (T2) or partial r = 0.13 (T1). That output would be
worse than no output. This is a standing constraint, enforced in the tool's own header and in its
output file header.

This does change the old plan's decision 5, which ended in "a ranked candidate list fit for
heterologous expression". A ranking is still produced, but it ranks **clade membership and
site completeness**, and it is labelled as such. It is not a predicted-activity ranking, because
nothing in the data supports one.

---

## 2. Decisions in force

| ID | Decision | Status |
|---|---|---|
| D1 | Eurotiales sourcing for the *Pichia* construct | **Moot.** It rested on the pocket score, which no longer predicts anything of interest |
| D3 | Novelty measured per protein and per clade | Mine; both, since per-protein is a by-product |
| D4 | **Clades I-V** (2025 scheme) define a subtype; Kinner I/II kept as a separate column | User, 2026-11-07 |
| D5 | AlphaFold3 acceptable (all non-commercial); Boltz-2 also to be tried | User, 2026-11-07 |
| D6 | G3 scope: all-BFD versus RNA-seq-followable clades | **OPEN - the one blocking decision** |
| G2 | T1 pivot, and report the T2 negative | User, 2026-11-07 |

Carried forward unchanged from the superseded sketch: decision 2 (report every class inclusively),
decision 3 (keep early-diverging lineages and add newly sampled Dikarya), decision 4 (build our own
tool; UPObase is offline), decision 7 (**replace HMM score competition with phylogenetic
placement** - now the top correctness item in section 3.1), decision 8 (do not contact UPObase
authors). Decision 1 is superseded by D4: Kinner Family I/II becomes a secondary column, not the
primary partition. Decision 6 (merge CPO-like into Family I) now applies only to that secondary
column.

Installed and verified 2026-11-07: `alphafold/3.0.0`, `alphafold/3.0.2`, `boltz/2.2.1`,
`boltzdesign1/0.4.1`, `hmmer/3.4`, `mafft/7.505`. **CAVER is not installed** and needs budgeted
setup; it is a desktop GUI application, not a module load.

---

## 3. The software: `upo_sites`

```
upo_sites/bin/upo_sites.py            annotate + assign clade, proteome scale
upo_sites/bin/validate_clade_loo.py   leave-one-out validation of the clade caller
```

**Pipeline.** hmmsearch PF01328 at the gathering threshold, then chunked `hmmalign` of the hits
together with the 19-enzyme calibration backbone, then read residues at fixed HspUPO positions,
then assign clade by nearest neighbour over match-state columns with identity and margin reported.

**Why it scales.** The HspUPO-position to HMM-match-state map is computed once from the profile, so
it is identical in every chunk and at any input size. Chunking is therefore safe and cost is
linear. MAFFT L-INS-i, used in the original pilot, is quadratic and cannot run at BFD scale.

**Measured properties.**

- 21 sites requested, **18 scoreable**. Positions **223, 228 and 232 fall in insert columns** of
  PF01328 and have no cross-sequence meaning. The tool detects and reports this itself, which
  independently confirms the reliability finding in `CALIBRATION_MATRIX.md` section 3.1.
- Pocket score is over the five reliable positions 73, 77, 81, 183, 187.
- On the 87-protein BFD sample it reproduces the independently validated distribution exactly
  (pocket score 0:3, 1:3, 2:9, 3:27, 4:33, 5:12).
- **Leave-one-out clade accuracy: 18/19 = 94.7 %.** The single failure is MphUPO, the sole Clade IV
  member; leave-one-out deletes its entire clade, so it cannot be called correctly by construction.
  On clades with two or more members the result is 18/18.

### 3.1 The honest limitation, and what it forces

**The validation does not cover the regime the tool runs in.** Correct leave-one-out calls span
52.7-96.2 % identity to the nearest reference. BFD proteins sit at roughly 25-35 % identity to
anything characterized. Nearest-neighbour assignment is validated at high identity and
**extrapolated** at BFD identity.

Margin does not rescue this. Correct calls had margins from 9.9 upward; the single wrong call had a
margin of 17.7. **Margin does not separate right from wrong in this set**, so the current
`--min-margin` default is a heuristic, not a calibrated threshold, and is labelled as such. On the
87-protein sample only 40 % of calls cleared the confidence flag.

This forces three things, in priority order:

1. **An identity floor below which the tool returns `unassigned` instead of a forced call.**
   A non-answer is better than a confident-looking wrong one.
2. **Replace nearest neighbour with phylogenetic placement** (`epa-ng` or `pplacer`) against a
   fixed reference tree of the 19 plus the wider `upo_id` reference set. This is what the
   superseded plan's decision 7 already called for, and it is now the single biggest correctness
   improvement available. Placement gives a likelihood-weighted assignment and a real uncertainty;
   nearest-neighbour identity gives neither.
3. **Re-validate by leave-one-out at matched identity**, subsampling references toward ~30 %
   identity rather than testing only the easy high-identity case.

---

## 4. Workstreams

### G1 - diversity, subtype annotation, novelty *(primary deliverable)*

1. Run `upo_sites` across all PF01328 proteins in BFD. **Measured 2026-11-07: 97,520 distinct
   proteins across 15,137 genomes** — not 23,451, which is the number of genomes with gene models;
   only 15,137 of them carry a PF01328 at all.
   **Tool built and tested at n = 87; scale-up is the next action.**
2. Dereplicate, then assign Clades I-V with the placement method of section 3.1, carrying Kinner
   Family I/II as a separate column per D4.
3. Describe clade distribution across the fungal tree. The Agaricomycete copy-number conflict is
   **resolved**: the maximum is **144 distinct PF01328 proteins** in one genome (*Sphaerobolus
   stellatus*, Agaricomycetes), followed by *Mycena* species at 66-106. The 147 figure counted
   `pfam.parquet` rows, which double-counts proteins carrying two PF01328 domains. Any copy-number
   query must use `count(distinct protein_id)`.
4. Novelty per protein and per clade (D3) by DIAMOND against protein **and** nucleotide/WGS
   databases. Protein-only inherits a deposition confound: absence from nr can mean "not deposited"
   rather than "novel".
5. Report lineages holding clade assignments not represented among characterized enzymes. These are
   the genuinely novel subtypes and the honest answer to the original goal 1.

**Known contamination risk, untested:** lichen and Glomeromycota genomes. Contig-level checks
before any lineage claim is reported.

### G2 - activity *(scope reduced by measurement)*

1. **Write up the T2 negative.** The aliphatic pocket does not transfer to aromatic ring substrates,
   measured at n = 19 against a widely cited structural feature. This is a publishable result.
2. **T1 as a clade-level prior only**, with observed in-sample ranges: Clade V 33-82, Clade I 0-53
   and bimodal, Clades II-IV 0-5. Always labelled n = 19, in-sample.
3. **Still worth extracting** from the 2025 supplement: Figures S6/S7 (full product heatmaps at two
   pH values) and S1 (the 19x19 identity matrix, stored as EMF so possibly text-extractable).
   These refine T1; they do not change the clade conclusion.
4. **Structure work (D5), clade representatives only.** Run AlphaFold3 and Boltz-2 on the same set
   and compare on channel geometry at the heme-contacting residues. Disagreement between them is
   the available estimate of how far channel geometry can be trusted. Neither is benchmarked here
   for cavity measurement. Measure GPU cost on the pilot set before any scale-up. Needs CAVER or a
   substitute installed.
5. **Do not build a per-protein activity score.** See section 1.1.

### G3 - regulation *(blocked on D6)*

Phase A, sequence only: extract promoters for clade-assigned UPO genes, find enriched motifs against
a matched background, test whether motif content tracks clade. Phase B: *Penicillium* and
*Aspergillus* RNA-seq to test whether predicted regulation matches observed expression.

**D6 decides the scale.** All-BFD is ~97,000 promoters across ~23,000 genomes. Restricting to clades
where RNA-seq follow-up is possible is a few hundred. These are different projects, not different
settings.

---

## 5. Order of work

1. ~~Retrieve the calibration matrix~~ - **done**, `CALIBRATION_MATRIX.md`.
2. ~~Build the Clade I-V backbone~~ - **done**, `data/clade_backbone/raw.fa`, all 19 recovered.
3. ~~Validate alignment against the published matrix~~ - **done**, 95.5 % overall, 100 % at the five
   pocket positions.
4. ~~Re-measure alignability under a method that scales~~ - **done**, `POCKET_TRANSFER.md` section 7.
5. ~~Run the T1 gate~~ - **done**, `CALIBRATION_MATRIX.md` section 8. T1 is clade, not pocket.
6. **Add the identity floor and `unassigned` state to `upo_sites`.**  <- next
7. **Scale `upo_sites` to full BFD.** **Measured 2026-11-07: 5,000 sequences annotate in 11 s
   wall, 124 MB peak RSS** (SLURM job 31846982). All 97,520 proteins therefore take roughly
   3 minutes in a single job, so the planned array is unnecessary and would be pure scheduling
   overhead. One job, 8 cores, 24 GB. The real cost in this workstream is the **sequence
   extraction** from 15,137 gzipped per-genome files, not the annotation.
8. Phylogenetic placement to replace nearest neighbour; re-validate at matched identity.
9. G1 steps 2-5.
10. G2 write-up of the T2 negative; structure comparison on clade representatives.
11. G3 Phase A once D6 is settled.

## 6. Risks

| Risk | Status | Handling |
|---|---|---|
| Clade caller extrapolates beyond its validated identity range | **Confirmed** | Identity floor; phylogenetic placement; re-validate at matched identity |
| Margin does not separate correct from incorrect calls | **Confirmed**, n = 19 | Threshold labelled heuristic, not calibrated |
| Clade IV has one member, so it cannot be validated or reliably called | **Confirmed** | Report Clade IV calls as provisional, always |
| PaDa-I activity data attached to the wild-type AaeUPO sequence | **Confirmed** | Kept distinct in the backbone and documented |
| Figure S4 values read off bar charts by eye | **Confirmed** | Used only for the T2 negative, which does not depend on precision |
| Figure S8 values are printed data labels | Not a risk | These are the authors' own numbers |
| Novelty-vs-nr inherits a deposition confound | **Confirmed** | Search nucleotide/WGS too |
| Lichen and Glomeromycota contamination | **Untested** | Contig-level checks before reporting |
| CAVER not installed | **Confirmed** | Budget install, or substitute MOLE / custom probe measure |
| Agaricomycete copy-number conflict (147 vs 144) | **RESOLVED 2026-11-07** | 144 distinct proteins is correct; 147 counted `pfam.parquet` ROWS, double-counting multi-domain proteins. Top genome is *Sphaerobolus stellatus* (144), then *Mycena* spp. (66-106). Always `count(distinct protein_id)` |
| 5 UniProt accessions from Table S1 return HTTP 200 with an empty body | **Confirmed** | Scripts must check body length, not status code; NCBI accessions are the reliable route |
