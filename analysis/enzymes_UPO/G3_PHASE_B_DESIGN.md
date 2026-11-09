# G3 Phase B: design, and its power limits, recorded BEFORE the results

Written 2026-11-09, before quantification finished, so the analysis cannot be reshaped to fit
whatever comes out.

## Data

Public RNA-seq, *Alternaria alternata* strain Z7, BioProject **PRJNA656145**, 15 runs,
3 biological replicates per condition:

| condition | n | tests |
|---|---|---|
| wild-type, untreated | 3 | control |
| wild-type + H2O2 | 3 | **oxidative stress** |
| wild-type + FeSO4 | 3 | **iron** |
| wild-type + NaCl | 3 | specificity: generic osmotic stress |
| wild-type + CuSO4 | 3 | specificity: another metal |

NaCl and CuSO4 are in the design on purpose. Without them, any H2O2 response could be a generic
stress response, and any FeSO4 response could be a generic metal response.

Reference: BFD `Alternaria_alternata_Y9` (LOCUSTAG `F75F20B0`), chosen because its gene ids carry
our clade calls. Quantification is salmon against its CDS transcripts.

## The power limit, stated up front

This is the weakest part of the design and it is not fixable with more reads.

Y9 contains **13 PF01328 proteins**. By phylogenetic placement, exactly **one** is clade IV
(`F75F20B0_001676`, reliable). A second (`F75F20B0_012248`) was called clade IV by the
nearest-neighbour caller but is `unassigned` by placement, which is the better-calibrated method.

So the clade-IV-specific hypothesis rests on **1 to 2 genes in one genome**. That is an
observation, not a test. It cannot produce a p-value worth reporting, and I will not manufacture
one by treating replicate libraries as if they were independent genes.

## What will actually be tested

1. **Class level, adequately powered.** Do the 13 Y9 UPOs as a group change expression under
   H2O2 or FeSO4 relative to control, compared against all other genes? This has real power.
2. **Clade IV, descriptive only.** Report the single clade IV gene's response alongside the other
   12, with no significance claim about the clade.
3. **Co-regulation.** Correlate UPO expression with CCAAT-complex (HAP) target genes across all
   15 samples. Correlation across 15 samples is computable; whether it means co-regulation is a
   separate question.

## Cross-strain caveat

Reads are strain Z7; transcripts are strain Y9. The salmon mapping rate per run is recorded and
is the QC. If mapping is poor the comparison is void, and that outcome will be reported rather
than worked around.

## Pre-registered interpretation rules

- A UPO class response present under H2O2 **and** NaCl is a generic stress response, not an
  oxidative-stress-specific one.
- A response under FeSO4 **and** CuSO4 is a generic metal response, not iron-specific.
- The single clade IV gene behaving differently from the other 12 is an anecdote to follow up,
  never a finding about clade IV.
- Phase A found the clade IV CCAAT enrichment is partly confounded with ancestry
  (p = 0.0005, see `G3_FINDINGS.md` caveat 5). Nothing here can undo that.
