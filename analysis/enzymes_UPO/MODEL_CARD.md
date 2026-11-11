# Model card: UPO active-site identification and T1 pocket scoring

> **2026-11-11:** clade-call reliability statements here ("precision 38 % below 55 %", "no clade
> call for 78.9 %") refer to the 19-reference caller. The 187-leaf placement caller (`upo_place/`)
> made no errors on the characterized enzymes at matched identity down to < 35 % and is now the
> production clade call; see `RECONCILIATION.md` §1. Tasks A-C below are otherwise unchanged:
> clade calls do not make Task C (activity outside clade V) reliable.

Version 1.0, 2026-11-08. Software: `upo_sites/bin/upo_scan.py`. Branch `upo-review-and-coordination`.

"Reliable" is not a property a model has in general. It is a property of a model **on a stated
task, measured**. This card states three tasks separately, because the evidence differs sharply
between them. Two are reliable. One is not, and no amount of restating makes it so.

---

## Task A — Find PF01328 heme-thiolate peroxidase proteins and read their active-site residues

**Status: RELIABLE. Use it.**

| property | value |
|---|---|
| Validation | 95/95 agreement with the published structure-based assignment matrix |
| Independent check | recovers AaeUPO's pocket `FAAFL` from sequence alone, matching published Fig. S2 for its variant PaDa-I |
| Conservation check | catalytic positions rank 6th, 12th and 15th of 186 match states |
| Coverage | all 97,520 PF01328 proteins in BFD, 15,137 genomes; set-equal to BFD's own annotation, 0 missing, 0 extra |
| Scoreable positions | 18 of 21; positions 223, 228, 232 fall in insert columns and are not comparable across sequences |
| Throughput | 87 proteins in ~20 s on 2 CPU; full BFD in 20 chunked SLURM jobs |

Evidence: `MODEL_VALIDATION.md` §A, `data/bfd_full/COVERAGE_CHECK.md`, `EXTERNAL_VALIDATION.md`.

## Task B — Rank proteins by similarity to the clade-V aliphatic pocket

**Status: RELIABLE for ranking, on the stated substrate. Use it to prioritize assays.**

| property | value |
|---|---|
| Rule | matches (0-5) to consensus `ILLAI` at HspUPO p73, p77, p81, p183, p187 |
| Cross-validation | leave-one-out, every parameter refit on the 18 training enzymes |
| LOO R² | +0.496 (clade-only model: +0.434; null: −0.114) |
| LOO AUC | 0.956 |
| Permutation p | 0.0005, 10,000 shuffles |
| Calibration set | n = 19, 10 active at > 10 % conversion |
| Alternatives tested and rejected | channel bottleneck (R² +0.108), corridor profile (−0.168), pocket+channel (+0.490) |

**Scope:** 4-propylguaiacol peroxygenase conversion, one substrate, one assay, one publication.
It says nothing about which rings an enzyme cleaves — that was tested and is a separate negative
(`RESULT_T2_negative.md`, r = +0.176).

## Task C — Predict whether an arbitrary protein is an active peroxygenase

**Status: NOT RELIABLE. Do not use it for this. This is measured, not cautious.**

| property | value |
|---|---|
| Out-of-clade AUC (clade V removed) | **0.556 — chance** |
| Out-of-clade actives available | **1** (GmaUPO-II, clade I, 49.3 % active, scored 1/5 — the rule misses it) |
| Per-protein confidence measure | attempted via applicability distance, **failed** (Spearman with skill −0.281, p = 0.88) |
| Structural rescue | attempted with 19 Boltz-2 models, **failed** (out-of-clade AUC 0.333) |
| Proteins outside the calibration domain | 69.6 % of BFD; 40.1 % of those scoring ≥ 3/5 |

**A high score is a clade-V hypothesis worth assaying. A LOW SCORE IS NOT EVIDENCE OF
INACTIVITY.** GmaUPO-II is the proof: 49.3 % active, 1 of 5 matches.

### Why this cannot be fixed by reanalysis

The out-of-clade activity distribution has a gap. GmaUPO-II sits at 49.3 %; the next non-clade-V
enzyme is MphUPO at 2.7 %. **No activity threshold between 3 % and 49 % changes the out-of-clade
positive count from one.** Validating out-of-clade detection requires out-of-clade positives, and
all published 4-propylguaiacol data contains exactly one. The constraint is arithmetic.

Routes closed by measurement, not by argument:

| route | outcome |
|---|---|
| second panel on the same substrate | does not exist (searched 2026-11-08) |
| external labels from other papers | invalid — labels are substrate-specific; PaDa-I, the field's benchmark peroxygenase, converts 4-PG at 1.7 % |
| applicability distance as confidence | failed, p = 0.88 |
| structural features, 19 Boltz-2 models | failed, out-of-clade 0.333 |
| a different activity threshold | cannot help; the gap is in the data |
| **assay 6 already-expressed enzymes on 4-PG** | **open — predictions pre-registered** |

---

## Intended use

Prioritizing candidates for experimental assay, and annotating active-site residues across large
fungal proteomes. Every output row carries `calib_domain`, so a user can see per protein whether
it falls inside the range the 19 enzymes span.

## Out-of-scope use

Calling an enzyme inactive. Claiming substrate scope beyond 4-propylguaiacol. Treating a clade
call below ~55 % identity to a reference as trustworthy — precision there is 38 %, and BFD's
median identity is 43.8 %. Treating `calib_domain` as a probability.

## Known failure modes

1. Misses out-of-clade actives (GmaUPO-II).
2. No clade call for 78.9 % of BFD; placement declines rather than guesses.
3. `calib_domain` is coverage, not confidence — it did not predict model error.
4. Degrades silently on a proteome unlike the calibration set; nothing warns beyond `calib_domain`.

## What would change this card

Assay results for the six pre-registered enzymes in `EXTERNAL_VALIDATION.md` (DchUPO and DcaUPO
predicted active; CraUPO, MweUPO, PabUPO-II, UPO1_Codo predicted inactive), or the 96-protein
panel in `ASSAY_PANEL.md`. A failure among the predicted actives would be the informative result:
it would confirm the rule overcalls outside clade V, which is exactly what the 0.556 warns about.

Until then, Task C stays marked not reliable. Relabelling it without new measurement would make
this document the opposite of what it is for.
