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

---

## Reference correction, 2026-11-08 — recorded because it nearly caused a false negative

The first quantification used `Alternaria_alternata_Y9` (BFD `F75F20B0`) and mapped at
**3.9 %**. Under the pre-registered QC rule that would have been reported as "cross-strain
quantification has failed, the comparison is void."

That verdict would have been wrong. The same 100,000 reads map at:

| reference | strain | mapping rate |
|---|---|---|
| Alternaria_alternata_Y9 | Y9 | **3.9 %** |
| Alternaria_alternata_JS-1623 | JS-1623 | 44.4 % |
| Alternaria_alternata_XJaa | XJaa | 50.6 % |
| **Alternaria_alternata_Z7** | **Z7** | **50.7 %** |

The Y9 annotation is the outlier, not the reads and not the biology. The reads are clean 150 bp
pairs at 53.3 % GC with a normal transcript-abundance profile, and the Y9 CDS file looks
superficially fine (12,614 transcripts, 17.2 Mb), so nothing but an explicit comparison would
have caught this.

**The reference is now `Alternaria_alternata_Z7` (BFD `F40104A9`), which is the SAME STRAIN the
RNA-seq came from.** The cross-strain caveat in the original design is therefore removed, not
merely reduced.

Consequences for the design above:
- Z7 carries **12** PF01328 proteins, and again exactly **one** clade IV by placement
  (`F40104A9_012798`, reliable). The power limit is unchanged: clade IV remains an observation on
  one gene, never a test.
- The mapping floor is set to 40 %. About 50 % is the realistic ceiling for a CDS-only reference,
  because RNA-seq also contains UTRs, intronic pre-mRNA, rRNA and ncRNA that a CDS file cannot
  represent. A floor of 50 % would have rejected a perfectly good run.
- The CCAAT-complex gene set is re-derived in Z7: `F40104A9_003727` (CBFB_NFYA) plus eight
  PF00808 histone-fold candidates.

The general lesson, which applies beyond this analysis: a QC rule that can only say "pass" or
"fail" will convert a broken reference into a biological conclusion. The floor needed a second
reference to be interpretable.
