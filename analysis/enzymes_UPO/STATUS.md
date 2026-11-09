# UPO project status

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

## What is bounded, with the bound measured

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
