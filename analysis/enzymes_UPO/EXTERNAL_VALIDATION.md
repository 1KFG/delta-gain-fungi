# External validation of the T1 score: why it is blocked, and the cheapest experiment that unblocks it

The identification model is not reliable for predicting activity outside clade V. Out-of-clade
AUC is 0.556 on n = 1. This document records what was tried to fix that computationally, why each
route is closed, and what the smallest piece of bench work is that would settle it.

Nothing here is an argument that the model is better than `MODEL_VALIDATION.md` says. It is an
account of where the ceiling comes from.

## Route 1: another panel assayed the same way — does not exist

The 19 calibration enzymes come from Barros-Reguera et al. 2025
([10.1186/s13068-025-02675-w](https://doi.org/10.1186/s13068-025-02675-w), PMC12306115), which
assayed 19 UPOs from five clades on 4-propylguaiacol. A literature search (2026-11-08) for a
second multi-enzyme panel measured on the same substrate with comparable units found none. That
paper is the only panel of its kind. Other 2024-2025 UPO characterizations exist — MfiUPO
(*Marasmius fiardii*), CmaUPO (*Coprinopsis marcescibilis*), bacterial UPOs — but each reports
its own substrates and units.

## Route 2: an external test set built from other papers' activity labels — INVALID, and here is the proof

The obvious workaround is to collect characterized UPOs from other publications and use a binary
"is it a peroxygenase?" label. **That would measure a different property, and the calibration set
itself proves it.**

**PaDa-I is the directed-evolution variant of AaeUPO (*Agrocybe aegerita*), the most studied UPO
in the field and the standard benchmark peroxygenase.** Its measured 4-propylguaiacol
peroxygenase conversion in the calibration data is **2 %, 0 % and 3 %** across the three pH
values. Mean 1.7 %. It is labelled **inactive**.

A binary label drawn from the general literature would call PaDa-I active, because on
ethylbenzene, naphthalene and many other substrates it is. T1 calls it inactive, correctly, for
4-propylguaiacol. Any "independent test set" assembled this way would score the T1 rule against
a label it was never meant to predict, and would report failure where there is none.

So the barrier is **label incommensurability**, not data scarcity. This is a stronger and more
specific statement than the earlier note that commercial panels withhold sequences. It cannot be
engineered around. Activity is substrate-specific, and T1 is one substrate.

## What the rule says about the characterized enzymes we do hold sequences for

`upo_scan.py` was run over the 13 named characterized references in `upo_id/models/`.

| enzyme | pocket | matches | call | T1 label available? | agrees? |
|---|---|---|---|---|---|
| HspUPO | ILLAI | 5 | clade_V_like | active, 44.0 % | yes |
| DspUPO | ILLAI | 5 | clade_V_like | active, 46.3 % (DspUPO-I) | yes |
| TteUPO | VLLAA | 3 | clade_V_like | active, 33.3 % | yes |
| AaeUPO | FAAFL | 0 | below_cut | via PaDa-I, inactive 1.7 % | yes |
| CviUPO | LLIKT | 1 | below_cut | inactive, 0.0 % | yes |
| MroUPO | VIAFV | 0 | below_cut | inactive, 1.0 % | yes |
| DcaUPO | LLLCI | 3 | clade_V_like | **none** | prediction |
| DchUPO | LLLAI | 4 | clade_V_like | **none** | prediction |
| CraUPO | FGVFF | 0 | below_cut | **none** | prediction |
| MweUPO | VIAFV | 0 | below_cut | **none** | prediction |
| PabUPO-II | FAAFL | 0 | below_cut | **none** | prediction |
| UPO1_Codo | VLQFL | 1 | below_cut | **none** | prediction |
| LfuCPO | IVAFL | 1 | below_cut | chloroperoxidase, not a UPO | n/a |

The six rows with a T1 label agree, but **five of those six are in the calibration set, so they
are in-sample and prove nothing about generalization.** They are shown as a consistency check on
the tool, not as validation of the score.

One of them is worth more than a consistency check. **AaeUPO is not in the calibration set; its
variant PaDa-I is.** Running the AaeUPO sequence through `upo_sites.py` independently recovers
the pocket `FAAFL`, which is exactly the pocket the published Figure S2 lists for PaDa-I. The
site-extraction pipeline reproduced a published active-site assignment from sequence alone,
through a different route than the paper used. That supports the site layer, which
`MODEL_VALIDATION.md` §A already measured at 95/95.

## Route 3: the cheapest experiment that would actually settle it

The 96-protein assay panel (`ASSAY_PANEL.md`) is the thorough answer and needs real bench work on
proteins that have never been expressed.

**There is a much cheaper first step.** Six characterized UPOs above already exist as expressed,
purified, published enzymes in other laboratories. None has a published 4-propylguaiacol number.
The rule makes a falsifiable, pre-registered prediction for each:

**Predicted ACTIVE on 4-propylguaiacol (> 10 % conversion):**
- **DchUPO** — pocket `LLLAI`, 4/5 matches
- **DcaUPO** — pocket `LLLCI`, 3/5 matches

**Predicted INACTIVE (< 10 %):**
- **CraUPO** (`FGVFF`, 0/5), **MweUPO** (`VIAFV`, 0/5), **PabUPO-II** (`FAAFL`, 0/5),
  **UPO1_Codo** (`VLQFL`, 1/5)

Assaying six existing enzymes on one substrate is a far smaller undertaking than a 96-well panel
of novel proteins, and it attacks the exact weakness. **The informative outcome is a failure.**
If DcaUPO or DchUPO comes back inactive, the rule is overcalling outside its calibration clade,
which is precisely what the 0.556 out-of-clade AUC warns about and what no computation here can
resolve.

These predictions are recorded before any such measurement, so they cannot be adjusted afterwards.

## Summary

| route | status |
|---|---|
| second panel, same assay | does not exist |
| external labels from other papers | invalid — PaDa-I proves labels are substrate-specific |
| applicability distance as per-protein confidence | tested, failed (p = 0.88), `MODEL_VALIDATION.md` §F |
| more positions / richer model | n = 19 with 10 actives; adding parameters overfits |
| **assay 6 existing enzymes on 4-PG** | **open, cheap, predictions recorded above** |
| assay the 96-protein panel | open, thorough, expensive |

The model remains what section D measured: reliable for finding PF01328 sites (95/95) and for
ranking clade-V-like pockets (LOO AUC 0.956), and **not** demonstrated to find active enzymes
outside clade V. That has not been relabelled, because no measurement taken since supports
relabelling it.
