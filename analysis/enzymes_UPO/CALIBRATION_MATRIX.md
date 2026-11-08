# The 2025 calibration set: what the supplementary data actually contains

Retrieved and analysed 2026-11-07. This file replaces the "retrieve supplementary tables S1, S2 and S4"
step in `DESIGN.md` §10. **That step is complete. Its outcome is negative for the planned G2 model.**

## 1. Source

Barros-Reguera M, Lopez-Tavera E, Schröder G, Ayuso-Fernández I, Eijsink VGH, Sørlie M.
*Structure-function relationships in unspecific peroxygenases revealed by a comparative study of their
action on the phenolic lignin monomer 4-propylguaiacol.*
**DOI [10.1186/s13068-025-02675-w](https://doi.org/10.1186/s13068-025-02675-w)**, PMC12306115.
*Biotechnology for Biofuels and Bioproducts*, 2025. NMBU (Ås, Norway).

Supplementary Material 1, `13068_2025_2675_MOESM1_ESM.docx`, 9.3 MB, downloaded from

```
https://static-content.springer.com/esm/art%3A10.1186%2Fs13068-025-02675-w/MediaObjects/13068_2025_2675_MOESM1_ESM.docx
```

Stored at `data/calib_dl/supp1.docx`; extracted figures at `data/calib_dl/figs/imageN.*`.

**Note for anyone repeating this:** the PMC mirror link
(`pmc.ncbi.nlm.nih.gov/articles/instance/12306115/bin/...`) returns an HTML interstitial, not the file.
Use the Springer static-content URL above.

## 2. First finding: there is no Table S4

The supplement contains **two tables only** — S1 (the enzyme list) and S2 (selected docking poses, chosen
by the authors "primarily for illustrative purposes"). Everything quantitative is in **figures S3–S13**,
as bar charts and heatmaps.

This matters for the design. There is no machine-readable per-enzyme calibration matrix to download. Any
numeric use of this dataset requires digitizing bar charts by eye, or contacting the authors. Every number
in section 4 below is a read-off, not a reported value.

## 3. Second finding: the anchor and the pocket definition are now confirmed and were incomplete

Figure S2 is the structure-based active-site alignment for all 19 enzymes. Its caption states verbatim:

> "Residues were selected based on the catalytic center of HspUPO. Residue numbering corresponds to
> HspUPO ... the structural location of each residue for clade V is described, with labels corresponding
> the aliphatic pocket (AP), the entrance (E), the catalytic pocket (CP), and the acid-base pair: Glu180
> for all UPOs, and His110 in short UPOs and Arg173 in long UPOs (Clade I)."

**HspUPO numbering is confirmed.** The inference in `POCKET_TRANSFER.md` §1 — recovered by testing all 90
reference sequences for the I/L/L/A/I/L pattern — was correct.

**But the aliphatic pocket is eight positions, not six.** The paper's own AP labels are:

| Location | Positions (HspUPO numbering) |
|---|---|
| Aliphatic pocket (AP) | **73, 77, 81, 183, 187, 228, 232**, plus **223** labelled AP/E |
| Entrance (E) | 80, 84, 172, 175 (and 223) |
| Catalytic pocket (CP) | 108, 176, 179 |
| Base | 110 (His in short, Gly in long) |
| Acid | 180 (Glu, universal across all 19) |

The pilot used six of these (73, 77, 81, 183, 187, 228) and missed **223 and 232**. The pilot's six were a
correct subset, so `POCKET_TRANSFER.md` is not wrong, but it is incomplete. Section 5 reports what changes.

The transcribed matrix is at **`data/calib_dl/active_site_S2.tsv`** (19 enzymes × 18 positions, plus the
paper's own "similarity relative to HspUPO" column). It was transcribed by eye from the figure image and
carries that warning in its header.

### 3.1 Validation: our alignment versus the published structure-based matrix

**Done 2026-11-07, superseding the five-enzyme spot check first written here.** All 19 calibration
sequences were retrieved (section 7.1), aligned with MAFFT L-INS-i, and read at the 18 Figure-S2 positions
using HspUPO numbering. Results in `data/clade_backbone/observed_S2_positions.tsv`.

**Agreement: 322 of 337 scorable residues = 95.5 %.** A sequence alignment recovers a 3DM structure-based
alignment at these positions. That is real support for the Tier 1 approach.

The 15 disagreements are not scattered. They concentrate almost entirely at the positions the pilot did
**not** use:

| Position | Location label | Mismatches / 19 | Which enzymes |
|---|---|---|---|
| **73, 77, 81, 183, 187** | AP (pilot's core) | **0** | — perfect agreement |
| 108, 110, 173, 176, 179, 180, 84 | CP / Base / Acid / E | 0 | — perfect agreement |
| 80 | E | 2 | AluUPO, AtuUPO (values swapped between them) |
| 172 | E | 1 | MroUPO |
| **175** | E | 3 | CviUPO, MorUPO, MroUPO |
| **223** | AP/E | 2 | AbrUPO-II, CviUPO |
| 228 | AP | 1 | MroUPO |
| **232** | AP | 6 | all 5 Clade I enzymes, plus MroUPO |

**This is the most useful result of the exercise.** The pilot's six positions (73, 77, 81, 183, 187, 228)
are exactly the subset that sequence alignment places reliably — five of them at 100 % agreement, and 228
at 18/19. The two positions I proposed adding in section 3, **223 and 232, are the two least reliable
columns in the whole matrix.**

Three distinct causes, all identified:

1. **Position 232, Clade I — a gap-placement disagreement, not a residue disagreement.** The paper places a
   gap at 228 for all five long Clade I enzymes and a residue at 232. Our alignment places the residue at
   228 and shifts 232. The residues themselves match once the offset is accounted for: CabUPO-III reads
   `228=V` in our alignment and `232=V` in the paper. So this is one ambiguous indel in the C-terminal
   region of long UPOs, not a mis-identified residue. It still means **232 cannot be scored for long UPOs
   by sequence alignment alone.**
2. **MroUPO is genuinely divergent** at 172, 175, 228 and 232 — four of the fifteen mismatches come from
   this one enzyme. The MroUPO 228 discrepancy flagged earlier is real and is part of a local cluster, not
   an isolated slip. MroUPO is only 239 aa, the second shortest in the set. Resolving it needs the
   structure (PDB 7ZBP), not more alignment.
3. **Position 80, AluUPO versus AtuUPO: the two values are transposed.** These rows are adjacent in the
   figure and the two proteins are near-identical in sequence. This is most likely a transcription slip on
   my side rather than an alignment error. I have **not** silently rewritten the transcription to match our
   alignment — that would be circular. `active_site_S2.tsv` records the figure as published, and this table
   records the discrepancy. Settle it by re-reading the figure at higher magnification if either position
   ever becomes load-bearing.

**Consequence for the design: keep the pilot's six positions. Do not add 223 and 232.** Section 3's claim
that the pocket "is eight positions, not six" is correct about the paper's structural definition but wrong
as advice for this project, because two of the eight are not recoverable from sequence. For structure-based
work in Tier 2 all eight are available; for the BFD-scale sequence scan, six is the defensible set.

## 4. Third finding: the aliphatic pocket does not predict naphthalene activity, now at n = 19

Figure S4 gives naphthalene activity for all 19 enzymes as a **continuous** value (percent relative to the
maximum), not as the binary presence/absence the main text implies. Digitized by eye, ±5 units:

| Enzyme | Clade | AP aliphatic /8 | Naphthalene | Best other assay |
|---|---|---|---|---|
| HspUPO | V | 6 | **100** | 54 |
| MroUPO | III | **7** | **100** | 100 |
| MorUPO | III | **7** | 95 | 50 |
| CabUPO-III | I | 4 | 61 | 47 |
| AtuUPO | V | **8** | 52 | 19 |
| DspUPO-I | V | **8** | 33 | 100 |
| CmaUPO-I | I | 3 | 21 | 29 |
| AbrUPO-I | V | 6 | 17 | 22 |
| AluUPO | V | **8** | 15 | 6 |
| MphUPO | IV | 6 | 5 | 3 |
| LspUPO-II | I | 3 | 4 | 0 |
| PaDa-I | I | 4 | 3 | 3 |
| AacUPO | V | **8** | 2 | 2 |
| GmaUPO-II | I | 4 | 2 | 3 |
| RneUPO | V | 6 | 2 | 5 |
| CviUPO | II | 6 | **0** | **100** |
| PanUPO | II | 6 | **0** | 11 |
| AbrUPO-II | V | **7** | 0 | 0 |
| TteUPO | V | **7** | 0 | 1 |

**Pearson r between aliphatic-pocket score and naphthalene activity = 0.176 across all 19.**
(Scored on the paper's 8 AP positions. On the 6 alignment-reliable positions of section 3.1 the picture is
the same: means of 5.00 versus 4.40. The conclusion does not depend on which pocket definition is used.)
Restricted to enzymes with measurable activity on at least one other assay, r = 0.231.

Mean pocket score, naphthalene ≥ 10 versus < 10: **6.33 versus 5.70** out of 8. On the six-position pilot
score: **5.00 versus 4.40**. There is no useful separation.

The pattern that does stand out runs the wrong way for the design. **AacUPO scores a perfect 8/8 and has
essentially no naphthalene activity (2). AluUPO scores 8/8 and reaches only 15. Meanwhile MroUPO and
MorUPO, both Clade III, are the joint highest on naphthalene — and Figure S13 states explicitly that Clade
III "lacks a constrained aliphatic pocket."** The feature and the target are decoupled in both directions.

This confirms, with the full calibration set rather than the n = 6 subset used in `DESIGN.md` §4.1, that
**the aliphatic pocket is a 4-propylguaiacol feature, not an aromatic-ring-hydroxylation feature.** The
earlier conclusion stands and is now better supported.

### 4.1 Two of the four "negatives" are uninformative

`DESIGN.md` §4.1 treated AbrUPO-II, CviUPO, PanUPO and TteUPO as four measured naphthalene negatives. The
supplementary data shows this is too generous. **AbrUPO-II and TteUPO are at or near zero on every one of
the five assays** (ABTS, NBD, 2,6-DMP, naphthalene, 4-NA). Figure S4's caption notes the enzymes were
"only partially purified" and that the graph "primarily illustrates differences in substrate specificity."
An enzyme that is flat across all five assays is most simply explained as poorly expressed or inactive
protein, not as an enzyme that is selectively unable to oxidise naphthalene. Figure S3's colony screen is
consistent with this for AbrUPO-II.

**That leaves exactly two interpretable negatives: CviUPO and PanUPO.** Both are demonstrably active
enzymes — CviUPO is the single most active on ABTS in the whole panel (100) — and both have zero
naphthalene.

### 4.2 And those two negatives are completely confounded with clade

CviUPO and PanUPO are **the entire membership of Clade II**. There are no other Clade II enzymes in the
calibration set, and no Clade II enzyme with naphthalene activity.

So the usable training signal for aromatic ring hydroxylation is: *n = 2 negatives, perfectly confounded
with a clade label.* Any classifier trained on this will learn "is it Clade II", not "does it hydroxylate
aromatic rings". It would look excellent in cross-validation and mean nothing.

## 5. What this does to the design

**Gate T2 cannot be passed on this dataset.** `DESIGN.md` §7.3 set Gate T2 as: does aromatic channel
lining separate the naphthalene negatives? It can be run, but it cannot produce a trustworthy answer,
because the only two interpretable negatives are one complete clade. This is a property of the data, not of
the method, and no amount of feature engineering fixes it.

The honest options, in the order I would take them:

1. **Report the negative and keep the descriptive work.** G1 (diversity, subtype assignment under Clades
   I–V, novelty) and G3 (regulation) do not depend on the predictive model at all. They remain fully
   deliverable. G2 becomes a *measured, documented* statement that the published structural feature does
   not transfer to aromatic substrates, which is a publishable result and is more defensible than a model.
2. **Widen the calibration set beyond this one paper.** Karich 2017 (14 PAHs, AaeUPO and MroUPO),
   Schmitz 2023 (AbrUPO on substituted benzenes), Linde 2020 (CviUPO kinetics). These use different hosts,
   assays and units, so merging them needs an explicit, written normalisation rule, and the combined set
   will still be small. This is the only route to a model, and it is a real project of its own.
3. **Do not ship a score.** A 0–8 pocket score computed across 97,520 BFD proteins would be easy to
   produce and would rank every fungal UPO by a feature now measured to have r ≈ 0.18 against the property
   of interest. That output would be worse than no output.

### 5.1 The CviUPO conflict is resolved

`DESIGN.md` listed an unresolved conflict: the 2025 study reports no CviUPO naphthalene activity, while
Linde 2020 reports a catalytic efficiency. Both can be true. The 2025 value is a relative colorimetric
slope on partially purified enzyme at pH 4.5, normalised to the best enzyme in the panel; Linde 2020 is
steady-state kinetics on purified protein. A low-but-nonzero k_cat/K_M can read as zero on a relative
colorimetric panel. **The 2025 naphthalene column measures "high enough to see against the best in the
panel", not "catalytically incapable".** This weakens the negatives further and reinforces section 4.2.

## 6. What is still worth extracting

Not yet digitized, and genuinely useful for G1/G3 rather than for the dead G2 model:

- **Figure S8** — peroxygenase fraction versus undetected product, per enzyme, at three pH values. This is
  the closest thing in the paper to the peroxygenase-vs-peroxidase axis (target T1 in `DESIGN.md` §5.2).
  Unlike naphthalene, this axis *is* what the aliphatic pocket was defined on, so T1 remains live.
- **Figures S6 and S7** — full product heatmaps at pH 3.5 and 5.5, with and without ascorbic acid.
- **Figure S11** — CAVER access-channel surfaces for Clade I, computed by the authors. Useful as a check
  target when CAVER is set up here.
- **Figure S1** — the 19 × 19 percent identity matrix (stored as `image1.emf`, a vector file, so it may be
  text-extractable rather than needing digitization).

## 7. Table S1, verbatim

Nineteen enzymes, selected from 2,044 UPOs in the 3DM database, biased toward solved crystal structures and
prior recombinant expression; MphUPO included as the one uncharacterized member. All expressed in
*Komagataella phaffii* (*Pichia*).

| Name | UniProt | NCBI | Organism | Long/short | MW (kDa) | Clade |
|---|---|---|---|---|---|---|
| AacUPO | A0A1L9X7H3 | XP_020060613.1 | *Aspergillus aculeatus* | Short | 28.4 | V |
| AbrUPO-I | A0A1L9U8A9 | OJJ67899.1 | *Aspergillus brasiliensis* | Short | 27.8 | V |
| AbrUPO-II | A0A1L9UNI9 | OJJ73116.1 | *Aspergillus brasiliensis* | Short | 27.3 | V |
| AluUPO | A0A7R7WEQ5 | BCS01637.1 | *Aspergillus kawachii* | Short | 27.7 | V |
| AtuUPO | A0A100IPW8 | GAQ45152.1 | *Aspergillus niger* | Short | 27.8 | V |
| CabUPO-III | A0A4Q2DAL3 | RXW15716.1 | *Candolleomyces aberdarensis* | Long | 39.4 | I |
| CmaUPO-I | A0A5C3KVD0 | TFK24496.1 | *Coprinopsis marcescibilis* | Long | 39.1 | I |
| CviUPO | UPI002016A77D | — | *Collariella virescens* | Short | 25.9 | II |
| DspUPO-I | A0A1Y2X6J4 | OTB17553.1 | *Daldinia* sp. EC12 | Short | 27.3 | V |
| GmaUPO-II | A0A067SWH8 | KDR72024.1 | *Galerina marginata* | Long | 38.3 | I |
| HspUPO | A0A1Y2TH07 | OTA57433.1 | *Hypoxylon* sp. EC38 | Short | 26.4 | V |
| LspUPO-II | A0A137QW42 | KXN91485.1 | *Leucoagaricus* sp. SymC.cos | Long | 38.7 | I |
| MorUPO | A0A9P7RSA8 | XP_043005295.1 | *Marasmius oreades* | Short | 25.0 | III |
| MphUPO | K2QYP1 | EKG15041.1 | *Macrophomina phaseolina* | Short | 27.0 | IV |
| MroUPO | UPI002016A75C | — | *Marasmius rotula* | Short | 25.9 | III |
| PaDa-I | — | — | *Agrocybe aegerita* | Long | 35.9 | I |
| PanUPO | B2B6C6 | XP_001911526.1 | *Podospora anserina* | Short | 28.3 | II |
| RneUPO | A0A1W2TUZ6 | GAP92448.1 | *Rosellinia necatrix* | Short | 26.4 | V |
| TteUPO | G2R8M6 | XP_003654578.1 | *Thielavia terrestris* | Short | 26.6 | V |

Correction to `DESIGN.md` §5.3 as first written: the Clade I member is **CabUPO-III**, not CabUPO-II.

These accessions make the clade assignment operational — all 19 can be pulled from UniProt/NCBI and used
directly as the placement backbone for assigning BFD proteins to Clades I–V.

### 7.1 Sequences retrieved — the Clade I–V placement backbone

Done 2026-11-07. All 19 are in `data/clade_backbone/raw.fa`, aligned in `backbone.afa`.

- 13 came from UniProtKB, 2 from UniParc (CviUPO, MroUPO — these have no UniProtKB entry).
- **5 UniProtKB accessions are now dead**: A0A100IPW8 (AtuUPO), A0A1Y2X6J4 (DspUPO-I), A0A1Y2TH07
  (HspUPO), A0A137QW42 (LspUPO-II) and K2QYP1 (MphUPO). The REST API returns HTTP 200 with
  `content-length: 0` in release 2026_04, which reads as success — any script must check for an empty body,
  not just the status code. All five were recovered from the NCBI accessions in Table S1.
- **PaDa-I has no accession in Table S1** because it is an engineered laboratory variant of the
  *Agrocybe aegerita* wild-type enzyme. The backbone carries the wild-type AaeUPO (UniProt B9W4V6) from
  `upo_id/ref/refs_all.fa`, labelled `PaDa-I_parent_AaeUPO`. **Its activity data in section 4 belongs to
  PaDa-I, not to this sequence.** Any model that uses PaDa-I's numbers against this sequence inherits the
  difference between a directed-evolution variant and its parent. Keep them distinct.

Sanity checks that passed: HspUPO is 261 aa, matching PDB 7O1R and `POCKET_TRANSFER.md` §1; all five long
Clade I enzymes are 371–387 aa and all fourteen short enzymes are 226–274 aa, consistent with the Long/Short
column and the MW column of Table S1.

Also retrieved for the PMC mirror problem noted in section 1: the NCBI accessions in Table S1 are the
reliable route for these sequences, not the UniProt ones.

---

## 8. The T1 gate, run 2026-11-07

Target T1 is peroxygenase-versus-peroxidase propensity — the axis the aliphatic pocket was actually
defined on. Data: Figure S8, transcribed to `data/calib_dl/t1_peroxygenase_S8.tsv`. Unlike Figure S4 these
are **printed data labels**, not bar heights measured by eye, so the numbers are the authors' own.

Feature: aliphatic pocket score over the five alignment-reliable positions (73, 77, 81, 183, 187), read
from the `hmmalign` backbone, which matches the published figure exactly at those positions.

### 8.1 Headline: the feature works

| Comparison | Pearson r | Spearman |
|---|---|---|
| Pocket (0–5) vs T1, max over pH | **+0.748** | +0.775 |
| Pocket (0–5) vs T1 at pH 5.5 | +0.719 | +0.746 |

Compare T2 (naphthalene): r = +0.176. **The pivot to T1 was correct** — the pocket predicts the property
it was defined on, and does not predict the one it was not.

### 8.2 But the feature is nearly collinear with clade, and clade is the better predictor

| Comparison | Pearson r |
|---|---|
| "is Clade V" vs T1 | **+0.835** |
| "is Clade V" vs pocket score | **+0.851** |
| **Partial r(pocket, T1 \| is Clade V)** | **+0.129** |

**Knowing the clade beats knowing the pocket, and the pocket adds almost nothing on top of it.**

Two checks confirm this is not an artifact of the partial correlation:

- **Within Clade V (n = 9) the pocket has almost no variance** — eight of nine score 5/5, one scores 4.
  The apparent within-clade r of +0.479 rests entirely on that single point. It is leverage, not a
  relationship.
- **Within the ten non-Clade-V enzymes, r = +0.001.** The pocket is flat against T1 there.

The case that settles it is **GmaUPO-II**: Clade I, pocket 3/5, and T1 = 53, the only high-peroxygenase
enzyme outside Clade V. CabUPO-III and PaDa-I are also Clade I with pocket 3/5 and score 0 and 3. The
pocket gives all three the same value and cannot explain the one that matters.

### 8.3 What this means for the build

T1 is predictable. **It is predictable from clade membership, not from the aliphatic pocket.** Any
per-protein T1 score built on pocket residues would be reporting clade through a more complicated route,
and would carry a spurious implication that the residues themselves were measured to matter.

With n = 19 spread over 5 clades, and clade explaining most of the variance, **it is not possible to
separate "this residue feature matters" from "this clade behaves this way".** That is a limit of the
calibration data and no modelling choice overcomes it.

The defensible deliverable therefore shifts:

1. **Primary product: reliable clade assignment.** Clade is what carries the measured signal, so a tool
   that places any fungal PF01328 protein into Clades I–V with a stated confidence is the useful thing to
   build, and it is exactly what G1 needs anyway.
2. **Site annotation reported as evidence, not as a prediction.** Emit the active-site residues, the
   pocket score and the catalytic-machinery flags per protein, because they are informative and cheap.
   Do not combine them into an activity number.
3. **T1 reported as a clade-level prior** with its observed range, explicitly labelled in-sample, n = 19.
   Clade V mean T1 = 53.8; Clade I = 11.8 and bimodal; Clades II–IV ≈ 0–5.
