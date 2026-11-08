# Pilot: does the Clade V "aliphatic pocket" transfer across the fungal tree?

Measured 2026-11-07. All numbers here were computed in this run, not quoted. Data and alignment are in
`data/pocket_test/` (zstd). Method is reproducible from the commands recorded in section 5.

## 0. Why this was run

The 2025 comparative study (`SUBSTRATE_CALIBRATION.md` section 12,
[DOI](https://doi.org/10.1186/s13068-025-02675-w)) reports that Clade V UPOs carry a conserved "aliphatic
pocket" next to the heme — residues **I73, L77, L81, A183, I187, L228** — and that this feature tracks
peroxygenase rather than peroxidase activity. That paper's enzymes are ascomycetes with high mutual identity.

Open risk: BFD proteins sit near 25 % identity to any characterized enzyme. If those six positions cannot be
located in a divergent sequence, the whole structural scoring plan fails at step one. This pilot tests that
and nothing else.

## 1. Identifying the anchor sequence

The paper gives residue numbers but the species name was stripped from the retrieved text. The anchor was
recovered by testing all 90 sequences in `upo_id/ref/refs_all.fa` for the pattern I73, L77, L81, A183, I187,
L228.

**Exactly one sequence matches all six: HspUPO** (*Hypoxylon* sp. EC38, PDB 7O1R, 261 aa). No other sequence
scored above 3 of 6. A 6/6 match by chance is implausible, so the numbering is HspUPO numbering.

Control positions taken from the same sequence: **Cys39** (the proximal Cys of the PCP motif) and
**Glu109/His110/Asp111** (the EHD metal site). These are the most conserved positions in the family, so they
set the ceiling for what alignment can achieve.

## 2. Sample

- Source: `pfam.parquet`, all PF01328 proteins, joined to `samples.parquet` for lineage.
- Grouping: phylum for non-Dikarya, class for Ascomycota and Basidiomycota.
- **3 proteins per lineage, hash-ordered with a fixed salt** (reproducible, not cherry-picked), from every
  group with at least 3 available.
- Result: **87 proteins across 29 lineages**, length 123–746 aa (median 396). All 87 sequences were recovered
  from `seqdb`; none missing.
- Lineages covered span the tree: Aphelidiomycota, Chytridiomycota, Monoblepharomycota, Zoopagomycota,
  Basidiobolomycota, Mucoromycota, Glomeromycota, 10 Basidiomycota classes, 12 Ascomycota classes.

Alignment: MAFFT L-INS-i (`--localpair --maxiterate 1000`) on the 87 BFD proteins plus 10 characterized
references. HspUPO ungapped positions were mapped to alignment columns, then read across every sequence.

"Aliphatic" is scored as **I, L, V, A or M**. That is my operationalization of "aliphatic pocket"; the paper
does not define the class membership this way, and section 4 treats this as a limitation.

## 3. Results

### 3.1 Reference controls — the method behaves as expected

| Reference | % id to HspUPO | pocket cols aligned | aliphatic | catalytic cols | pocket residues |
|---|---|---|---|---|---|
| DspUPO | 72.7 | 6/6 | 6/6 | 4/4 | ILLAIL |
| DchUPO | 67.9 | 6/6 | 6/6 | 4/4 | LLLAIL |
| DcaUPO | 66.7 | 6/6 | 5/6 | 4/4 | LLLCIL |
| TteUPO | 59.1 | 6/6 | 6/6 | 4/4 | VLLAAL |
| MroUPO | 36.2 | 6/6 | 5/6 | 4/4 | VIAFVL |
| MweUPO | 35.5 | 6/6 | 5/6 | 4/4 | VIAFVL |
| CviUPO | 33.0 | 6/6 | 4/6 | 4/4 | LLIKTM |
| LfuCPO | 28.0 | 6/6 | 5/6 | 4/4 | IVAFLV |
| AaeUPO | 27.6 | 5/6 | 3/6 | 4/4 | FAAFL– |

This reproduces the paper's claim independently: the close Clade V relatives (DspUPO, DchUPO, DcaUPO,
TteUPO, all 59–73 % identity) keep the pocket nearly intact, while AaeUPO — the long-type Family II enzyme
that is mainly a peroxidase on phenolic substrates — has only 3 of 6 aliphatic and loses one column entirely.
The measure separates the enzymes in the published direction. That is the pilot's internal validation.

### 3.2 BFD sample — the headline numbers

**n = 87, mean identity to HspUPO about 25 %.**

| Quantity | Result |
|---|---|
| Pocket columns alignable (non-gap) | **95.6 %** (499 / 522) |
| Catalytic-site columns alignable (ceiling) | **97.1 %** |
| Pocket positions that are aliphatic | **68.0 %** |
| Proteins retaining the catalytic Cys | **84 / 87** (96.5 %) |

**The positions transfer. The chemistry does not.**

Alignability at the pocket columns (95.6 %) is within 1.5 points of the catalytic-site ceiling (97.1 %), so
MAFFT places these columns about as reliably as it places the most conserved residues in the family. The
feature is therefore computable at BFD scale.

But only 68 % of those mapped positions are aliphatic. So "has the aliphatic pocket" cannot be a binary call
across fungi. It must be a graded score from 0 to 6, and the typical BFD protein will be partial.

### 3.3 Per lineage (n = 3 each — indicative only)

| Lineage | mean % id | pocket aligned /6 | aliphatic /6 |
|---|---|---|---|
| Ustilaginomycetes | 23.8 | 6.0 | **6.0** |
| Aphelidiomycota | 28.6 | 6.0 | 5.7 |
| Sordariomycetes | 27.3 | 6.0 | 5.3 |
| Dacrymycetes | 24.9 | 6.0 | 5.0 |
| Tremellomycetes | 24.5 | 6.0 | 5.0 |
| Zoopagomycota | 23.4 | 6.0 | 5.0 |
| Agaricomycetes | 27.7 | 6.0 | 4.7 |
| Chytridiomycota | 26.6 | 5.0 | 4.7 |
| Coniocybomycetes | 25.2 | 6.0 | 4.7 |
| Dothideomycetes | 26.8 | 6.0 | 4.7 |
| Glomeromycota | 27.2 | 6.0 | 4.7 |
| Microbotryomycetes | 23.0 | 6.0 | 4.7 |
| Taphrinomycetes | 24.7 | 6.0 | 4.7 |
| Basidiobolomycota | 27.6 | 6.0 | 4.3 |
| Orbiliomycetes | 27.7 | 5.7 | 4.3 |
| Mucoromycota | 26.4 | 6.0 | 4.0 |
| Monoblepharomycota | 27.3 | 5.0 | 4.0 |
| Saitoellomycetes | 28.8 | 6.0 | 4.0 |
| Leotiomycetes | 26.4 | 6.0 | 3.7 |
| Exobasidiomycetes | 26.3 | 5.0 | 3.7 |
| Wallemiomycetes | 24.5 | 6.0 | 3.7 |
| Geoglossomycetes | 24.9 | 6.0 | 3.3 |
| Pezizomycetes | 25.4 | 6.0 | 3.3 |
| Lecanoromycetes | 26.6 | 5.0 | 3.0 |
| Xylonomycetes | 24.4 | 6.0 | 3.0 |
| Arthoniomycetes | 22.9 | 4.7 | 2.3 |
| Candelariomycetes | 29.4 | 5.0 | 2.3 |
| Eurotiomycetes | 21.3 | 5.0 | 2.3 |
| Pucciniomycetes | 19.5 | 6.0 | 2.3 |

### 3.4 Two things this rules out, and one it raises

**Ruled out: length class does not explain the variation.** Short-type (<320 aa, n=37) mean 4.03/6 versus
long-type (>=320 aa, n=50) mean 4.12/6. Essentially identical. The low-scoring lineages are not simply
Family II proteins.

**Ruled out: it is not an alignment failure.** The columns align at the catalytic-site ceiling. Where the
pocket score is low, the residues really are non-aliphatic.

**Raised, not established: Ustilaginomycetes scores 6.0/6**, the highest of all 29 lineages, at only 23.8 %
identity to HspUPO. Independently, *Moesziomyces aphidis* (Ustilaginomycetes) is reported to use an
extracellular UPO with a stated preference for PAHs
([DOI](https://doi.org/10.1016/j.jhazmat.2024.134437), `SUBSTRATE_CALIBRATION.md` section 6). Those two facts
are consistent, but **n = 3 and this is a single observation — it is a hypothesis to test, not a finding.**

**Note for the Pichia plan: Eurotiomycetes scores 2.3/6**, among the lowest. If the aliphatic pocket turns
out to be the feature that matters, *Aspergillus* and *Penicillium* copies are not where the best candidates
will be. This does not affect *Pichia* as an expression host; it affects where the genes come from.

## 4. Limitations — read before using these numbers

1. **n = 3 per lineage.** The per-lineage column in 3.3 is indicative only. No confidence interval is
   meaningful at this size. The overall figures (n = 87) are more solid than any single row.
2. **"Aliphatic = ILVAM" is my definition**, not the paper's. The paper describes a pocket geometry observed
   in structures. A position can hold a non-aliphatic residue and still form part of a hydrophobic cavity,
   and the converse. This score is a proxy for the published feature, not the feature.
3. **Alignment columns are not structural superposition.** A residue can occupy the correct alignment column
   and still point away from the cavity in three dimensions. Confirming that needs structures, which is
   exactly what Tier 2 of the design is for. This pilot establishes only that the positions are locatable.
4. **One alignment, one method.** No sensitivity test against a different aligner, no bootstrap, no
   structure-guided alignment. MAFFT L-INS-i on 97 sequences was used as-is.
5. **3 of 87 proteins lack the catalytic Cys** and are probably not functional peroxygenases. They were left
   in the totals. Removing them would raise the pocket figures slightly.
6. **The sample is PF01328 at the gathering threshold**, so it inherits that net's blind spots — notably
   fragments and highly divergent members that GA excludes.

## 5. Reproducing this

1. Sample: duckdb over `pfam.parquet` filtered to `PF01328`, joined to `samples.parquet`, grouped as in
   section 2, `row_number() over (partition by grp order by hash(protein_id||'s7')) <= 3`.
2. Sequences: read `seqdb/proteins/<sha1(LOCUSTAG)[0:2]>/<LOCUSTAG>.proteins.fa.gz` directly with the stdlib
   `gzip` module (bgzip is gzip-compatible; no samtools needed).
3. Align: `module load mafft/7.505; mafft --localpair --maxiterate 1000 --thread 8`.
4. Map HspUPO ungapped positions 39, 73, 77, 81, 109, 110, 111, 183, 187, 228 to alignment columns and read
   every sequence at those columns.

Stored inputs and outputs: `data/pocket_test/{sample,bfd_meta}.tsv.zst`, `bfd_sample.fa.zst`,
`refs_sel.fa.zst`, `all.afa.zst`.

## 6. What this means for the design

- The structural-feature approach is **viable**: the positions can be found in any fungal PF01328 protein.
- The pocket must be scored **0–6, graded**, never as a binary "has pocket / lacks pocket".
- A score of 6/6 at 24 % identity (Ustilaginomycetes) shows the feature genuinely can be conserved across
  deep divergence, so a high score in a distant lineage is not automatically an artifact.
- The **next validation step is structural**, not sequence-based: build models for a handful of high-scoring
  and low-scoring BFD proteins and check with CAVER whether the high scorers really do present a constrained
  channel with a hydrophobic cavity. Until that is done, the 0–6 score is a sequence proxy of unknown
  structural fidelity.

---

## 7. ADDENDUM 2026-11-07 — re-measured under `hmmalign`

Section 3 used MAFFT L-INS-i. That is O(n²) and cannot run on 97,520 proteins, so the headline numbers had
to be re-measured with a method that scales before they could be trusted for the full BFD scan. Done here.

**Method.** `hmmalign --amino` (HMMER 3.4) of the same 87 BFD proteins plus the 19 calibration enzymes
(`data/clade_backbone/raw.fa`, see `CALIBRATION_MATRIX.md` §7.1) against `upo_id/models/PF01328.hmm`.
HspUPO aligns at its full 261 residues. Only **match-state** columns are comparable across sequences;
insert columns are not homologous and were excluded. Inputs and output in `data/hmmalign_test/`.

### 7.1 Position 228 is not usable at BFD scale

**HspUPO residue 228 falls in an insert column of PF01328.** It is not a match state, so it has no
cross-sequence meaning under a profile alignment. MAFFT could place it; `hmmalign` cannot, and neither can
any profile-based scan of the full dataset.

The other five pocket positions (73, 77, 81, 183, 187) and all four catalytic control positions (Cys39,
Glu109, His110, Asp111) are match states, and HspUPO reads `I L L A I` and `C E H D` at them as expected.

Combined with `CALIBRATION_MATRIX.md` §3.1, which found 223 and 232 to be the least reliably aligned
columns in the matrix, the usable sequence-level pocket is now **five positions: 73, 77, 81, 183, 187.**
The paper's structural definition has eight. Three of those eight are not recoverable from sequence.

### 7.2 The profile alignment is exact on the calibration set

At those five positions across all 19 calibration enzymes, `hmmalign` agrees with the published
structure-based Figure S2 matrix at **95 of 95 residues — 100 %.** No mismatches.

This is a stronger validation than section 3.1's reference-control table. A profile HMM alignment
reproduces a 3DM structure-based alignment exactly at the positions the project depends on.

### 7.3 Re-measured BFD numbers

| Quantity | MAFFT (§3.2, 6 positions) | `hmmalign` (5 positions) |
|---|---|---|
| Pocket columns alignable | 95.6 % (499/522) | **94.7 %** (412/435) |
| Catalytic columns alignable (ceiling) | 97.1 % | **95.7 %** (333/348) |
| Mapped pocket positions aliphatic | 68.0 % | **71.4 %** (294/412) |
| Proteins retaining catalytic Cys | 84/87 (96.5 %) | **83/87** (95.4 %) |

**The section 3 conclusions survive the method change.** Alignability stays within one point of the
catalytic ceiling, so the positions are locatable at scale. The chemistry still does not transfer: about
three in ten mapped positions are non-aliphatic, so the pocket remains a graded score, never a binary call.
The small differences between the two columns are within what a different aligner and a different position
count would produce; neither method is being presented as more correct.

Pocket score distribution across the 87 BFD proteins, out of 5: mean 3.38 —
0:3, 1:3, 2:9, 3:27, 4:33, 5:12. Only 12 of 87 are complete.

### 7.4 What this does not change

The aliphatic pocket still has no measured power to predict aromatic ring hydroxylation. See
`CALIBRATION_MATRIX.md` §4: Pearson r = 0.176 against naphthalene activity across all 19 calibration
enzymes. Making the feature cheaper to compute does not make it predictive. It is retained for target T1
(peroxygenase versus peroxidase propensity), which is the property it was originally defined on.
