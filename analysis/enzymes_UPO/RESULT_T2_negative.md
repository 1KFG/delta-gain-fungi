# Result: the UPO "aliphatic pocket" does not predict aromatic ring hydroxylation

**Status: negative result, measured. 2026-11-07.**
Supporting data and derivations: `CALIBRATION_MATRIX.md` sections 3-5 and 8.

## Summary

A conserved "aliphatic pocket" next to the heme has been proposed as the structural feature that
separates high-performing unspecific peroxygenases (UPOs) from poor ones. We tested whether that
feature predicts activity on **aromatic ring hydroxylation**, the reaction of interest for
hydrocarbon degradation. **It does not.**

Across the 19-enzyme comparative panel of Barros-Reguera et al. 2025
([DOI 10.1186/s13068-025-02675-w](https://doi.org/10.1186/s13068-025-02675-w)), the correlation
between aliphatic-pocket completeness and naphthalene activity is **Pearson r = +0.176**.

The feature is not useless; it is specific to a different reaction. On
peroxygenase-versus-peroxidase propensity toward 4-propylguaiacol — the property it was originally
defined on — the same score gives r = +0.748. **The pocket predicts what it was defined on, and
does not transfer to aromatic ring substrates.**

## Why this was worth testing

The pocket was defined on 4-propylguaiacol: a phenolic ring carrying an aliphatic propyl chain,
where the reaction is C-H abstraction **on that chain**. An aliphatic pocket holding an aliphatic
chain against the heme is a direct mechanistic story. Hydroxylating a bare aromatic **ring** is
different chemistry, and the published work does not claim otherwise. The question was whether the
feature generalises anyway. It does not.

## Evidence

### 1. The correlation is absent

Pocket score (0-8, the paper's own AP positions) against naphthalene activity, read from Figure S4:

| r | n |
|---|---|
| **+0.176** | 19 |

Restricted to the ten enzymes with measurable activity on at least one other assay, r = +0.231.
Scoring on the five alignment-reliable positions instead of eight gives the same picture
(group means 5.00 vs 4.40). **The conclusion does not depend on the pocket definition used.**

### 2. The feature fails in both directions

| Enzyme | Clade | Pocket /8 | Naphthalene |
|---|---|---|---|
| AacUPO | V | **8** (perfect) | **2** |
| AluUPO | V | **8** (perfect) | 15 |
| AbrUPO-II | V | 7 | **0** |
| TteUPO | V | 7 | **0** |
| MroUPO | III | 7 | **100** |
| MorUPO | III | 7 | **95** |
| HspUPO | V | 6 | **100** |

Perfect pockets with no activity, and top activity in Clade III — which the same paper's Figure S13
states explicitly **"lacks a constrained aliphatic pocket."** The feature and the target are
decoupled in both directions.

### 3. The negatives will not support a model even if one wanted to build it

The paper reports four enzymes with no naphthalene activity: AbrUPO-II, CviUPO, PanUPO, TteUPO.
Two of these are uninformative. **AbrUPO-II and TteUPO are at or near zero on all five assays**
(ABTS, NBD, 2,6-DMP, naphthalene, 4-NA). The caption notes the enzymes were "only partially
purified". An enzyme flat across every assay is most simply poorly expressed protein, not an enzyme
selectively unable to oxidise naphthalene.

That leaves **two interpretable negatives, CviUPO and PanUPO** — and those two are **the complete
membership of Clade II**. Label and clade are perfectly confounded. Any classifier trained here
learns "is it Clade II", would cross-validate well, and would mean nothing.

### 4. A related conflict, resolved

An apparent contradiction — the 2025 study reports no CviUPO naphthalene activity while Linde et al.
2020 report a catalytic efficiency — is not a contradiction. The 2025 value is a **relative
colorimetric slope on partially purified enzyme, normalised to the best enzyme in the panel**;
Linde 2020 is steady-state kinetics on purified protein. A low but non-zero k_cat/K_M reads as zero
on a relative panel. The 2025 naphthalene column measures "high enough to see against the best in
the panel", not "catalytically incapable". This weakens the negatives further.

## What we are not claiming

- Not that the aliphatic pocket is wrong or unimportant. It correlates with the property it was
  defined on (r = +0.748 against peroxygenase propensity).
- Not that no sequence or structural feature can predict aromatic hydroxylation. We claim that
  **this** feature does not, and that **this** dataset cannot settle the general question, because
  it has two interpretable negatives and they are one clade.
- Not a conclusion about purified-enzyme kinetics. Figure S4 is a relative colorimetric assay on
  partially purified protein.

## Limitations

1. **n = 19**, one laboratory, one expression host (*Komagataella phaffii*), one substrate panel.
2. **Naphthalene values were digitized by eye** from a bar chart (±5 units); the supplement contains
   no table of them. The conclusion is insensitive to this: r = 0.18 does not become useful under
   any plausible read-off error, and the qualitative failures in section 2 are large.
3. Pocket scores use "aliphatic = ILVAM", our operationalization of a geometric feature. A residue
   can sit in the right alignment column and point away from the cavity.
4. Three of the paper's eight AP positions (223, 228, 232) fall in insert columns of PF01328 and
   cannot be scored across divergent sequences at all.

## Consequence for the project

The planned per-protein activity predictor for aromatic hydroxylation **is not buildable on the
available calibration data, and was cancelled rather than shipped on a feature with r = 0.18.**
Work moved to target T1 (peroxygenase propensity), where a further test showed the pocket is nearly
collinear with clade and adds partial r = +0.129 beyond it. The deliverable is therefore a clade
annotator with transparent site evidence, not an activity score. See `PLAN.md` sections 1 and 1.1.
