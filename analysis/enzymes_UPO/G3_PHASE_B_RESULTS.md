# G3 Phase B results

## Mapping-rate QC (strain Z7 reads against the strain-matched Z7 transcripts)

| run | condition | % mapped |
|---|---|---|
| SRR12147525 | control | 51.1 |
| SRR12147531 | H2O2 | 59.2 |
| SRR12147532 | H2O2 | 59.7 |
| SRR12147533 | H2O2 | 59.3 |
| SRR12147536 | control | 54.3 |
| SRR12147537 | control | 53.5 |
| SRR12629930 | NaCl | 63.8 |
| SRR12629931 | NaCl | 62.5 |
| SRR12629932 | NaCl | 62.2 |
| SRR12629933 | FeSO4 | 57.8 |
| SRR12629934 | FeSO4 | 57.5 |
| SRR12629935 | FeSO4 | 57.3 |
| SRR12629936 | CuSO4 | 67.2 |
| SRR12629937 | CuSO4 | 64.8 |
| SRR12629938 | CuSO4 | 65.7 |

Median mapping rate 59.3 %. 15 of 15 runs pass the 40 % floor.

12954 transcripts quantified across 15 libraries.

12 of 12 UPO genes present in the quantification.

## Test 1: are the UPOs as a class shifted, versus all other genes?

Mann-Whitney U on per-gene log2 fold change, UPO set against every other expressed gene.

| treatment | UPO median log2FC | genome median log2FC | p (two-sided) |
|---|---|---|---|
| CuSO4 | +0.193 | +0.515 | 0.4667 |
| FeSO4 | -0.087 | +0.049 | 0.2907 |
| H2O2 | +0.000 | +0.121 | 0.1444 |
| NaCl | +0.377 | +0.156 | 0.79 |

**Pre-registered rule:** a response under H2O2 that is also present under NaCl is a generic stress response. A response under FeSO4 also present under CuSO4 is a generic metal response.

## Test 2: the individual UPO genes (descriptive)

| gene | clade | trust | pocket | CuSO4 | FeSO4 | H2O2 | NaCl |
|---|---|---|---|---|---|---|---|
| F40104A9_002483 | unassigned | unassigned | 2 | +0.09 | -0.73 | +0.81 | -1.56 |
| F40104A9_004496 | V | extrapolated | 3 | +0.14 | +0.00 | +0.22 | +0.16 |
| F40104A9_012798 | IV | reliable+provisional_cladeIV | 2 | +3.03 | +1.87 | +0.19 | +1.04 |
| F40104A9_000984 | II | extrapolated | 1 | +2.19 | +0.20 | +0.11 | +1.59 |
| F40104A9_010964 | unassigned | unassigned | 2 | +0.41 | +0.05 | +0.06 | +0.18 |
| F40104A9_000876 | unassigned | unassigned | 2 | -0.23 | -0.33 | +0.02 | -0.15 |
| F40104A9_008941 | unassigned | unassigned | 2 | +0.24 | +0.05 | -0.02 | +0.58 |
| F40104A9_003903 | unassigned | unassigned | 2 | -0.14 | -0.53 | -0.08 | -0.43 |
| F40104A9_010939 | unassigned | unassigned | 1 | +0.55 | -0.17 | -0.18 | +0.93 |
| F40104A9_007204 | unassigned | unassigned | 2 | +0.11 | -1.08 | -0.21 | +1.79 |
| F40104A9_000498 | unassigned | unassigned | 2 | -1.40 | -3.78 | -0.47 | -1.24 |
| F40104A9_007621 | V | extrapolated | 2 | +0.59 | +4.04 | -0.55 | +0.92 |

The clade IV row is ONE gene. It is shown, not tested. No claim about clade IV as a clade can rest on it.

## Test 3: co-regulation with the CCAAT-binding complex

9 CCAAT-complex genes. Spearman across 15 libraries.

| UPO gene | clade | mean rho vs CCAAT | mean rho vs random genes |
|---|---|---|---|
| F40104A9_000498 | unassigned | +0.054 | -0.087 |
| F40104A9_000876 | unassigned | +0.109 | +0.023 |
| F40104A9_000984 | II | -0.029 | +0.225 |
| F40104A9_002483 | unassigned | +0.229 | +0.034 |
| F40104A9_003903 | unassigned | +0.121 | +0.072 |
| F40104A9_004496 | V | +0.071 | +0.174 |
| F40104A9_007204 | unassigned | -0.117 | +0.056 |
| F40104A9_007621 | V | -0.121 | +0.025 |
| F40104A9_008941 | unassigned | -0.034 | +0.177 |
| F40104A9_010939 | unassigned | -0.105 | +0.155 |
| F40104A9_010964 | unassigned | -0.031 | +0.130 |
| F40104A9_012798 | IV | -0.013 | +0.186 |

UPO mean rho vs CCAAT +0.011, vs random +0.097.

Correlation across 15 libraries that differ by treatment is produced by any shared response to those treatments. It is not evidence of a regulatory link, and the difference from the random-gene baseline is the only part worth reading.


---

## Conclusion (written after seeing the results; the tests above were pre-registered)

**Phase B is a negative result, and the negative is interpretable.**

The QC no longer limits the conclusion. Median mapping 59.3 %, range 51.1-67.2 %, all 15 runs
pass. Reads and reference are the same strain, Z7. The earlier 3.99 % rate came from a broken Y9
annotation, not from biology.

**Test 1. No treatment moves the UPOs as a class.** All four p values are between 0.14 and 0.79.
Under H2O2, the treatment with the strongest prior, the UPO median log2FC is +0.000 while the
genome median is +0.121: the UPOs move less than the average gene, not more. Under CuSO4 and
FeSO4 the UPO median is below the genome median. The only treatment where the UPO median exceeds
the genome median is NaCl, with p = 0.79, which is no effect.

So the Phase A promoter enrichment does not predict a coordinated transcriptional response to
oxidative stress or to iron in this genome. The hypothesis that UPOs are a stress- or
metal-regulated class is not supported here.

**Test 2 is one gene and stays an anecdote.** The clade IV gene F40104A9_012798 has the highest
CuSO4 fold change of the 12 (+3.03) and the second highest FeSO4 (+1.87). That is a single gene
with no control for multiplicity across 12 genes and 4 treatments. It is written down as
something to look at, not as a result about clade IV. The pre-registered rule stands.

**Test 3 argues against the CCAAT link rather than for it.** Mean UPO-to-CCAAT Spearman rho is
+0.011. The random-gene baseline is +0.097. The UPOs are LESS correlated with the CCAAT-binding
complex than with 300 arbitrary genes. Nine of the twelve UPOs have a lower rho against CCAAT
than against the baseline.

### What this does and does not retire

It does NOT show the CCAAT motif is non-functional. One genome, one laboratory, four treatments,
15 libraries. The motif could be used under conditions not sampled here, in a species not
sampled here, or at a timepoint not sampled here.

It DOES remove the easiest positive explanation. Combined with the promoter-phylogeny control in
G3_FINDINGS.md caveat 5, which showed clade IV CCAAT carriers are more closely related than
chance (z = +3.44, p = 0.0005), the regulatory reading of the Phase A enrichment now has one
confound confirmed and no expression support. The enrichment is still real as a sequence
observation. It should not be described as evidence of shared regulation.
