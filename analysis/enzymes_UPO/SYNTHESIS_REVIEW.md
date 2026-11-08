# Synthesis review: all plans and reviews against the current work

Written 2026-11-08 from the laptop session `novelty-enzymes-50`, at the user's request, after
committing the 2026-11-07 HPCC work (`c090229`) and the previously-unshared review (`4dcd560`).

**Note for the session working in this directory:** I did not edit any of your files. This document
and the two commits are all I added. I could not reach you by message — you run on `skylark` and
are not addressable from the laptop's session list — so this file is the channel.

Documents read in full: `PLAN.md`, `PLAN.superseded.md`, `DESIGN.md`, `CALIBRATION_MATRIX.md`,
`POCKET_TRANSFER.md`, `RESULT_T2_negative.md`, `G1_BFD_SUMMARY.md`, `SUBSTRATE_CALIBRATION.md`,
`upo_sites/README.md`, `upo_reg/README.md`, `REVIEW_PLAN_superseded.md`.
Numbers labelled **measured here** were computed by me against `data/bfd_full/out/*.tsv.zst`
during this review; everything else is quoted from the documents.

---

## 1. Headline assessment

**The project is in better shape than either the superseded plan or the review that attacked it.**
Two activity gates were run and both came back negative, and the response in each case was to
narrow the claim rather than to rescue the model. `PLAN.md` §1.1 ("what is deliberately not being
built") and the refusal written into `upo_sites/README.md` are the strongest things in the
repository. A per-protein activity score across 97,520 proteins would have been easy, publishable-
looking and wrong, and it was declined on measured grounds.

The product is now a **clade annotator with site evidence**, which is what the data supports.

Three things now hold it back, in order:

1. **The 19-sequence backbone is too small for the job**, and the fix the review identified a day
   ago — a published 174-sequence curated set and a live database — is still untouched.
2. **`DESIGN.md`, the longest document and the one a newcomer reads first, is stale** and tells the
   superseded story.
3. **The confidence tiering, which is the project's main defence, has a load-bearing soft spot at
   Clade IV** that nothing currently flags in the summary tables.

---

## 2. Scorecard: the review against the current work

`REVIEW_PLAN_superseded.md` was written against the old sketch and never reached this work. Most of
it has been overtaken. Mapping each point:

| # | Review's point | Status now |
|---|---|---|
| 1 | EGD/EHD is not a UPO/CPO discriminator; use the acid-base pair | **Resolved by a better route.** Motif-based UPO/CPO calling was abandoned entirely. The acid-base positions are now scored directly (Glu180 universal; His110 short / Arg173 long, HspUPO numbering) |
| 2 | No engagement with Zamocky 2022 prior art | **STILL OPEN — and now the biggest gap.** See §3.1 |
| 3 | RedoxiBase is live and is the curated source | **STILL OPEN.** See §3.1 |
| 4 | Tier 1 as defined is nearly empty; SP gating is a Dikarya trait | **Moot.** Scope moved to all-BFD; tiering replaced by identity-confidence tiers |
| 5 | Absence claims used the wrong denominator | **Moot by rescope, and denominator discipline has improved markedly** — 97,520 proteins / 15,137 genomes / 23,451 with gene models are now kept distinct, and the 147-vs-144 copy-number error was found and fixed (`count(distinct protein_id)`) |
| 6 | jackhmmer is the wrong instrument; calibrate thresholds against decoys | **Half resolved.** jackhmmer dropped for `hmmalign`. But the underlying prescription — set thresholds by ROC against characterized positives and a decoy set — is **unimplemented**, and both the 30 % identity floor and `--min-margin` are explicitly labelled heuristics. See §3.2 |
| 7 | Structure prediction budgeted but never specified | **Partly.** D5 now names AlphaFold3 and Boltz-2 with a real comparison criterion (agreement on channel geometry). Still no acceptance criteria, and CAVER is not installed |
| 8 | Novelty metric: use global alignment and patristic distance, not local DIAMOND identity | **STILL OPEN.** G1 step 4 still specifies DIAMOND. The review's argument is now *confirmed by the project's own data* — see §3.3 |
| 9 | RNA-seq support unmeasured; promote miniprot to a gene-model check | **STILL OPEN, and the rescope made it worse.** See §3.4 |
| 10 | Glomeromycota short length is biology, not truncation | Resolved; consistent with the current data |
| 11 | Chytridiomycota/oomycete question needs an outgroup | **STILL OPEN.** No oomycete or non-fungal sequence is in the backbone |
| 12 | Reproducibility: missing extraction query, inconsistent denominators | **Largely resolved.** The tool is documented, the LOO test is runnable in one command, extraction is now `bin/extract_pf01328.py`, denominators are reconciled |

Six of twelve resolved or moot, two partly, four still open.

---

## 3. The still-open items that matter

### 3.1 Prior art and the reference set — one action fixes both

**Measured here:** the strings `Zamocky`, `RedoxiBase`, `Peroxibase` and `oomycete` appear in
**none** of the current documents. `Faiza` appears only in `upo_id/METHOD.md`.

This is no longer just a citation problem. The project's single worst limitation is that clade
assignment rests on **19 reference sequences**, against which BFD proteins sit at 25–42 % identity,
so 87 % of calls are extrapolations. Zamocky et al. 2022 (Antioxidants 11:1011, PMC9138132)
published **174 curated full-length sequences** of this superfamily, and RedoxiBase
(`peroxibase.toulouse.inra.fr`) is the live curated database behind it.

Folding those in attacks four problems at once:

- **It raises the backbone from 19 to ~174**, which is the only cheap route to lifting proteins out
  of the `extrapolated` bin. Nothing else on the roadmap does this.
- **It supplies the outgroup** the review asked for in point 11, including the oomycete sequences
  needed to test the Chytridiomycota placement.
- **It is the prior art a reviewer reaches for first.** Zamocky published a Mucoromycota-basal
  clade, a distinct Glomeromycota short subfamily, and placed *Spizellomyces* chytrid enzymes inside
  the oomycete family. `G1_BFD_SUMMARY.md` §4 now reports Mucoromycota (370), Glomeromycota (124)
  and Chytridiomycota (188) distributions with no reference to any of it.
- **It converts a liability into a positive control.** Recovering Zamocky's published clades from
  97,520 proteins would be real validation. Failing to recover them would be a finding. Either beats
  silence.

Note the Glomeromycota distribution is striking on its face: **124 proteins, 71 Clade V, 40 Clade
III, 0 Clade I.** A phylum with no Clade I at all is either a genuine result or an artifact of
reference coverage, and the published Glomeromycota subfamily is the obvious thing to test it
against.

**This is the highest-value unexploited action in the project.** It is also cheap: a download and a
backbone rebuild.

### 3.2 The thresholds are still uncalibrated

`PLAN.md` §3.1 and `upo_sites/README.md` are both candid that the 30 % floor and `--min-margin` are
heuristics, and that margin did not separate correct from incorrect calls at n = 19. That honesty is
right, but it leaves the tiering resting on two unvalidated numbers that partition 97,520 proteins.

The review's prescription still applies and is still the fix: calibrate against decoys. A sampled
set of PF01328 proteins that are confidently *not* UPOs gives a negative distribution, and the floor
becomes an ROC operating point instead of a round number. The project already has a natural decoy
source in its own data — see §4.2.

### 3.3 The novelty metric will not work as specified

G1 step 4 still specifies DIAMOND identity. **Measured here**, per-clade mean identity to the
nearest characterized enzyme is 37.5–42.8 % for Clades I, II, III and V, and 24.9 % for unassigned.
That is exactly the flat band the review predicted. Local identity will not rank anything
meaningfully inside it.

The planned move to phylogenetic placement (`epa-ng`/`pplacer`) already produces the right
measure as a by-product — **patristic distance** to the nearest characterized enzyme. Do the novelty
metric in the same pass rather than as a separate DIAMOND step.

### 3.4 Gene-model artifact control was dropped in the rescope, and the exposure grew

The superseded plan had a serious artifact-control section — partial models, missing start/stop,
internal stops, fused models, intron-position concordance, cross-strain comparison. The review
singled it out as the best part of the document. **It does not survive into `DESIGN.md` or
`PLAN.md`.** The current risk table carries only "lichen and Glomeromycota contamination, untested".

Meanwhile the exposure went up by an order of magnitude: the analysis now spans **15,137 genomes**
of self-generated gene models instead of ~800.

This matters directly for the headline numbers. `G1_BFD_SUMMARY.md` §5 reports the full catalytic
set intact in **35.4 %** of proteins and the proximal Cys present in **91.4 %**. Some unknown part
of the missing two-thirds will be truncated or mis-called gene models rather than genuine
divergence, and nothing currently separates the two. Until it does, "65 % lack the full catalytic
set" is not safely a biological statement.

The cheap control is already paid for: `gene_transcripts` in BFD carries `is_partial`,
`has_start_codon`, `has_stop_codon`. Joining those to the output is a single query and would let
every catalytic-machinery number be reported QC-clean versus flagged.

---

## 4. Findings from this synthesis, not raised in any review

### 4.1 `DESIGN.md` is stale and contradicts `PLAN.md`

`DESIGN.md` §10 step 2 reads "**G2 RESOLVED: pivot to target T1**", and calls the Figure S8
digitization "the gate on whether T1 is viable".

That gate was then run, and recorded in `CALIBRATION_MATRIX.md` §8: T1 **also** fails as a
per-protein feature. Raw r = +0.748, but partial r given clade = **+0.129**, and +0.001 within the
ten non-Clade-V enzymes. `PLAN.md` §1 Q3 states this correctly. `DESIGN.md` does not.

So the project's longest document (36 KB), explicitly framed as the intent document, still presents
a pivot that measurement has since overtaken. `PLAN.superseded.md` got an honest supersession
header; `DESIGN.md` needs the same treatment or a corrected §10.

### 4.2 Clade IV supplies 39 % of all validated-range calls, and is the one clade that cannot be validated

This is the sharpest issue I found.

**Measured here**, per-clade identity to the nearest characterized enzyme:

| clade | n | mean identity | n at identity >= 52.7 |
|---|---|---|---|
| IV | 6,922 | **62.1** | **4,968** |
| V | 24,847 | 42.8 | 2,637 |
| II | 3,113 | 39.4 | 403 |
| I | 33,424 | 37.9 | 3,731 |
| III | 10,154 | 37.5 | 1,079 |
| unassigned | 19,060 | 24.9 | 0 |

Clade IV's mean identity is ~20 points above every other clade, and **4,968 of the 12,818
validated-range calls in the entire dataset — 38.8 % — are Clade IV**.

Clade IV has exactly one characterized member, MphUPO. `PLAN.md` §6 and the LOO test both state
that it cannot be validated by construction, because leave-one-out deletes its entire clade. So the
largest single contributor to the "validated" bin is the clade with **zero** validation support.

The `+provisional_cladeIV` suffix exists and is good. But `G1_BFD_SUMMARY.md` §6's catalytic-core
table credits Clade IV with **4,565 validated-range** proteins — the largest of any clade, larger
than Clade V's 1,951 — with no note that every one is provisional. Anyone reading that table for a
candidate list would start exactly where the evidence is weakest.

I did not determine *why* Clade IV identity runs so high. Mean length (271.9 aa) is close to Clade V
(271.2), so a short-overlap inflation is not the obvious explanation. **Measured here**, MphUPO is
the nearest reference for 7,978 proteins. Two readings: MphUPO genuinely sits at the centre of a
large conserved fungal subfamily, or the identity calculation treats it differently from the others.
Worth one check before any of those 4,565 proteins are used.

### 4.3 The real working denominator is ~11,000, not 97,520

Worth stating once, plainly, because three different denominators circulate.

From `G1_BFD_SUMMARY.md` §6, the catalytic-core set is 34,561 proteins, of which the validated-range
subset by clade sums to **11,242**. That — not 97,520, and not 34,561 — is the set with both intact
catalytic machinery and a clade call inside the validated identity range. Subtract the provisional
Clade IV contribution (4,565, §4.2) and the well-supported core is **~6,700**.

6,700 of 97,520 is 6.9 %. That is a perfectly good result and a much more defensible headline than
any number ten times larger. It should be the stated denominator for downstream work.

### 4.4 Two internal numbers look contradictory and need one sentence of reconciliation

`upo_sites/README.md`: "On an 87-protein BFD sample, **only 1 of 87** reached 52.7 %" — 1.1 %.
**Measured here** on the full run: **12,818 of 97,620 = 13.1 %**.

A twelve-fold difference. The likely explanation is sampling design — the 87 were drawn 3 per
lineage across 29 lineages, deliberately spread, while BFD is dominated by Sordariomycetes and
Dothideomycetes which sit closer to the characterized references. If so it is not an error, but both
numbers now live in the same repository and look like a contradiction. One sentence fixes it.

### 4.5 Benign, but worth knowing

**Measured here:** the concatenated output is 97,620 lines, of which **100 are comment headers**
(5 `#` lines per chunk × 20 chunks). Real proteins = 97,520, matching the summary. No data is lost,
but any consumer must strip `#` lines as well as repeated `sample` headers, or every per-chunk
header becomes a phantom record. Worth a line in the README.

### 4.6 The Clade I caveat deserves promotion

`G1_BFD_SUMMARY.md` §6 contains the most important caveat in the document: Clade I calls carry the
universal Glu180 at only **16.8 %**, while all five characterized Clade I enzymes have it. The
document's own reading — that the low-identity Clade I bin is absorbing proteins that are not UPOs,
because the long characterized enzymes are the nearest available match for divergent sequences — is
almost certainly right, and is supported by Clade I's mean length of 437.8 aa against 271–309 for
every other clade.

But Clade I is **33,424 proteins, 34 % of everything**, and it is reported without that caveat in
the §4 lineage table, where Dothideomycetes' 15,054 Clade I calls are the largest single number in
the document. The caveat should travel with the number.

This also supplies the decoy set §3.2 needs: low-identity Clade I proteins lacking Glu180 are a
ready-made negative distribution for calibrating the identity floor.

---

## 5. Recommended order

1. **Mark `DESIGN.md` superseded or correct its §10.** Minutes. Until then the repository tells two
   different stories about T1.
2. **Bring in Zamocky 2022 Table S1 and RedoxiBase; rebuild the backbone at ~174 sequences.**
   Highest value by a wide margin: it attacks the extrapolation problem, supplies the outgroup, and
   closes the prior-art gap in one action.
3. **Check Clade IV** (§4.2) before anything downstream consumes the validated-range set, and add
   the provisional flag to the summary tables.
4. **Join `gene_transcripts` QC flags** to the output and re-report the catalytic-machinery numbers
   QC-clean versus flagged (§3.4).
5. **Carry the Clade I Glu180 caveat into §4**, and state the ~11,000 / ~6,700 working denominator
   (§4.3, §4.6).
6. Then phylogenetic placement, with novelty as patristic distance in the same pass (§3.3), and
   decoy-calibrated thresholds using the Clade I negatives (§3.2, §4.6).

Items 1, 3, 5 and 4.4–4.5 are documentation and checking, and could be done today. Item 2 is the
one that changes the science.

---

## 6. What is genuinely good, and should not be lost

- **Both gates were run and both were reported negative.** T2 at r = +0.176 and T1 at partial
  r = +0.129 are the two results that most of this field would have quietly buried.
- **The refusal to emit a per-protein activity score**, written into the tool's own header and
  output file, with the reason stated as a measurement rather than as caution.
- **The tool detects its own unscoreable positions.** 223, 228 and 232 fall in PF01328 insert
  columns and the tool reports this at startup — independently confirming the reliability finding in
  `CALIBRATION_MATRIX.md` §3.1, which was derived a different way.
- **`clade_T1_insample` is labelled as a property of the calibration set, not a prediction about the
  row it sits on.** That is a distinction that routinely gets lost in output tables.
- **The confidence tiering exists at all**, and `below_identity_floor` returns `unassigned` rather
  than forcing a call. §4.2 is a gap in it, not a refutation of it.
- **The 147-vs-144 copy-number error was caught** and traced to counting `pfam.parquet` rows instead
  of distinct proteins.
- **Chunked `hmmalign` with the position map computed once from the profile.** This is what makes
  the whole thing linear and chunk-safe, and it is the reason the full run costs minutes.

---

## 7. Verification notes

Computed by me during this review, against `data/bfd_full/out/*.tsv.zst` as of 2026-11-08 09:40:

- per-clade n, mean/max identity, and count at identity >= 52.7 (§3.3, §4.2 table)
- nearest-reference usage counts (MphUPO 7,978; LspUPO-II 19,452)
- per-clade mean protein length (§4.2, §4.6)
- total rows 97,620; comment headers 100; real proteins 97,520 (§4.5)
- absence of `Zamocky`/`RedoxiBase`/`Peroxibase`/`oomycete` across all current documents (§3.1)

Not verified by me: everything quoted from the 2025 paper, the r values in
`CALIBRATION_MATRIX.md`, the LOO result, and the Figure S2 transcription. The Clade IV
identity anomaly (§4.2) is reported as an observation with two candidate explanations; I did not
determine which holds.
