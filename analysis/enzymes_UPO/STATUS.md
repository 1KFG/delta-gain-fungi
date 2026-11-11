# UPO project status

> **2026-11-11:** the clade-calling numbers below (78.9 % unassigned, precision 38 % below 55 %) are
> for the 19-reference caller. The production clade call is now placement on the 187-leaf Zamocky
> 2022 reference (`upo_place/`, column `place_our_clade`): 42 % of BFD unassigned, and at matched
> identity 0 errors on the 19 characterized enzymes down to < 35 % identity (13/13 at 45-55 %).
> Caveats and the comparison: `RECONCILIATION.md` §1.

**Model card: `MODEL_CARD.md`** — the model stated as three separate tasks with measured
operating characteristics. Two are reliable (site finding, clade-V-like ranking). Predicting
activity outside clade V is not, and §I of `MODEL_VALIDATION.md` shows why reanalysis cannot
fix it.

**2026-11-08.** One page. Full reasoning is in the documents linked from each line.

## What is finished and defensible

| | Result | Where |
|---|---|---|
| **G1 diversity** | 97,520 PF01328 proteins in 15,137 genomes; **34,561 (35.4 %) are UPO-like** by catalytic core | `G1_FINDINGS.md`, `G1_BFD_SUMMARY.md` |
| **G1 novelty** | All 7,055 representatives have an `nr` hit. Apparent 18 % novelty is **deposition**: GCF 2.0 % vs GCA-only 18.7 %, z = 7.4 | `G1_FINDINGS.md` §5b |
| **G1 contamination** | Checked. UPO genes sit on contigs **more** gene-rich than their genome's typical contig | `G1_FINDINGS.md` §6 |
| **G2 T2** | Aliphatic pocket does **not** predict aromatic ring activity: r = **+0.176**, n = 19 | `RESULT_T2_negative.md` |
| **G2 T1** | Predictable, but by **clade, not pocket**: partial r = **+0.129** given clade | `CALIBRATION_MATRIX.md` §8 |
| **G3 Phase A** | **Clade IV is the only clade whose motifs replicate across independent lineages and survive both controls**, including a **CCAAT box** naming to the CCAAT-binding complex twice (HAP3, HAP5) | `G3_FINDINGS.md` |
| **Site identification** | **Validated**: 95/95 agreement with the published structure-based matrix; catalytic positions rank 6/12/15 of 186 match states by conservation | `MODEL_VALIDATION.md` §A |
| **G3 Phase B** | **Negative, and interpretable.** 15 Z7 RNA-seq libraries, median mapping 59.3 %. No treatment shifts the 12 UPOs as a class (p 0.14-0.79). UPOs correlate LESS with the CCAAT complex (+0.011) than with random genes (+0.097) | `G3_PHASE_B_RESULTS.md` |
| **Software** | `upo_scan.py`: one command, any proteome -> sites + score + applicability domain. Regression-tested 87/87 against the full BFD pipeline | `upo_sites/README.md` |

## What is bounded, with the bound measured

**The identification model is NOT reliable for finding active enzymes, and this has not been
relabelled.** Out-of-clade AUC is 0.556 on n = 1 active. A per-protein confidence measure was
attempted via applicability distance and **failed its test** (Spearman with model skill -0.281,
p = 0.88); the failure is recorded in `MODEL_VALIDATION.md` §F rather than dropped. What the
model IS reliable for is stated and measured: finding PF01328 sites (95/95), and ranking
clade-V-like pockets (LOO AUC 0.956). 69.6 % of BFD lies outside the calibration domain.

**Clade assignment is reliable above ~55 % identity and nowhere else.** Precision 100 % / 88 % /
38 % at 70–100 / 55–70 / 45–55 % identity. BFD median is **43.8 %**. Of 3,052 assigned
representatives, **2,088 clear the threshold** and the tool says which. (`MODEL_VALIDATION.md` §C)

Three attempts to break that ceiling, each measured rather than assumed, all negative:

- **Reference expansion** with high-confidence BFD members: precision **80 % → 66 %**. Inferred
  labels propagate into topology. (`PLACEMENT_FINDINGS.md`)
- **Phylogenetic placement**: better calibrated than nearest neighbour (57 % unassigned vs 19.5 %,
  92.6 % agreement) but does not extend the reliable range.
- **Clade IV's rank**: unresolved at 20–28 % bootstrap. The *group* is real (652/789 confirmed);
  its status as a separate clade is not. (`PLACEMENT_FINDINGS.md`)

**Root cause, quantified:** three of five clades have one or two reference enzymes, and 71 of 114
matched-identity test cases were impossible because removing a reference's close relatives removed
its whole clade.

## What is open right now

Nothing. Every job in the standing goal has finished.

| Job | State |
|---|---|
| Full-BFD phylogenetic placement (31868164) | DONE, 20/20 chunks, 97,520 queries |
| Boltz-2 modelling of all 19 (31870133) | DONE, and it FAILED validation; the 14 models are discarded |
| Final merged per-protein table | DONE: `data/bfd_full/upo_sites_bfd_final.tsv.zst`, `FINAL_TABLE_REPORT.md` |

### The deliverable

`data/bfd_full/upo_sites_bfd_scored.tsv.zst` — one row per PF01328 protein in BFD,
97,520 rows. Filter on `clade_trust`, rank on `t1_pocket_matches`.

**Shortlist: 5,115 proteins across 1,378 species** carry both the catalytic core and
≥ 3/5 pocket matches. 1,246 of them sit outside clade V, including 29 in early-diverging
phyla. This is a candidate list for assay, not a prediction of activity: the score reaches
LOO AUC 0.956 overall but 0.556 with clade V removed. See `T1_SCORE_REPORT.md`.

| tier | n | % |
|---|---|---|
| `unassigned` (placement declined) | 76,913 | 78.9 |
| `extrapolated` | 10,152 | 10.4 |
| `reliable` | 10,455 | 10.7 |

Site residues and `catalytic_core` are validated for every row. The clade call is
bounded, which is why the tier column exists. No activity score is emitted.

**The Boltz-2 run has a built-in kill switch.** 5 of the 19 have crystal structures. If predicted
channel widths do not reproduce the experimental ones (AaeUPO 4.23, MroUPO 5.19, HspUPO 4.87,
CviUPO 3.71, TteUPO 5.40 Å), the other 14 are discarded and the answer is that this cannot be done
with predicted structures.

## The channel-width question, honestly

It has moved twice today and both moves are recorded:

1. CAVER **3.0.2** gave ~1 Å bottlenecks and I concluded the ring rule was not transferable. That
   was **wrong** — 3.0.2 has a documented bug planting the start vertex "in too narrow cavity".
2. CAVER **3.0.3** reproduces the published ordering (MroUPO wider than AaeUPO) and gave
   r = **+0.933** against naphthalene on 4 enzymes.
3. Adding the fifth structure (**TteUPO, 9IZ8**) — the **widest channel measured, zero activity** —
   dropped it to r = **+0.452**.

TteUPO can be set aside on grounds recorded *before* that measurement (flat on all five assays, so
most likely poorly expressed). But that is a post-hoc exclusion rescuing a correlation, and n = 4
is not evidence. **This is what the Boltz-2 run is for: n = 19 instead of 5.**
(`CHANNEL_FINDINGS.md`)

## Decisions in force

D4 Clades I–V (user) · D5 AF3 and Boltz-2 both acceptable (user) · D6 all-BFD core, clade-stratified
(user) · G2 T1 pivot + report T2 negative (user) · Phase B RNA-seq: run on public SRA data (user) — **complete, negative**.

## The cheapest experiment that would settle the model question

Six characterized UPOs already exist as expressed, purified enzymes in other laboratories and
have **no published 4-propylguaiacol number**. Pre-registered predictions, recorded before any
measurement (`EXTERNAL_VALIDATION.md`):

- **active** (> 10 % conversion): DchUPO (`LLLAI`, 4/5), DcaUPO (`LLLCI`, 3/5)
- **inactive** (< 10 %): CraUPO (0/5), MweUPO (0/5), PabUPO-II (0/5), UPO1_Codo (1/5)

Six existing enzymes on one substrate is far smaller than the 96-well panel and attacks the
exact weakness. A failure on DcaUPO or DchUPO is the informative outcome: it would confirm the
rule overcalls outside clade V, which is what the 0.556 out-of-clade AUC warns about.

Why no computational route remains: structural features from the 19 Boltz-2 models were tested
and rejected (out-of-clade AUC 0.333 for channel, against 0.556 for the residue score; adding
them makes the model worse) - `STRUCTURE_VS_RESIDUES.md`; a second panel on this substrate does not exist; an external
set built from other papers' labels is invalid, because PaDa-I (the AaeUPO variant, the field's
benchmark peroxygenase) converts 4-PG at 1.7 % and T1 correctly calls it inactive; and
applicability distance failed as a confidence measure (p = 0.88).

## The one thing that would change everything

**19 characterized enzymes.** Every limit above traces to it, and it was tested rather than
assumed: expansion with inferred labels makes the caller worse, and the structural route is
calibrated on four. More characterized enzymes — Zamocky 2022 / RedoxiBase, or new biochemistry —
is the only route past it.

## Two tooling traps, recorded so they do not recur

- **`~/.local` numpy 2.5.2 shadows any module-provided Python environment** on this account. Boltz-2
  then fails at import with a message blaming numba. `PYTHONNOUSERSITE=1` fixes it; the environment
  was never broken.
- **CAVER writes `whole_clustering.csv` into the working directory.** Run it from a scratch
  directory, not the repo.
