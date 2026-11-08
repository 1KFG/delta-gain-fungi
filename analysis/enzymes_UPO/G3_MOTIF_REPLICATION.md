# G3 Phase A: cross-lineage motif replication

A motif is reported as **clade-linked** only if it is significant on STREME's held-out test set (test p < 0.05) in two or more INDEPENDENT lineages for the same clade, and the two are matched by Tomtom. Single-lineage motifs are lineage elements.

**Significance filter matters here.** STREME run with `--nmotifs 8` always reports 8 motifs per contrast whether or not they are significant. Across the 17 contrasts it reported 136 motifs, of which only those passing the held-out test are considered below.

STREME p-values are not used for ranking; see the module docstring for why.


## Clade I

Tested in 2 lineages: Agaricomycetes, Lecanoromycetes

- Agaricomycetes vs Lecanoromycetes: no matching motifs

**Clade I: 0 motif(s) replicated across lineages, from 1 lineage pair(s).**

## Clade II

Tested in 2 lineages: Agaricomycetes, Dothideomycetes

- Agaricomycetes vs Dothideomycetes: no matching motifs

**Clade II: 0 motif(s) replicated across lineages, from 1 lineage pair(s).**

## Clade III

Tested in 3 lineages: Agaricomycetes, Lecanoromycetes, Sordariomycetes

- Agaricomycetes vs Lecanoromycetes: no matching motifs
- Agaricomycetes vs Sordariomycetes: no matching motifs
- Lecanoromycetes vs Sordariomycetes: no matching motifs

**Clade III: 0 motif(s) replicated across lineages, from 3 lineage pair(s).**

## Clade IV

Tested in 4 lineages: Dothideomycetes, Eurotiomycetes, Leotiomycetes, Sordariomycetes

- **Dothideomycetes vs Eurotiomycetes: 3 matching motif pair(s)**
  - `1-CGARTACGCGTAY` (Dothideomycetes) matches `5-YTCGTACTCGTACY` (Eurotiomycetes), q = 0.00884448
  - `2-AGWGCGAGWRCGAR` (Dothideomycetes) matches `5-YTCGTACTCGTACY` (Eurotiomycetes), q = 8.65371e-05
  - `3-AGCCAATGARA` (Dothideomycetes) matches `2-AGCTCATTGGGTCGS` (Eurotiomycetes), q = 0.0091604
- **Dothideomycetes vs Leotiomycetes: 3 matching motif pair(s)**
  - `1-CGARTACGCGTAY` (Dothideomycetes) matches `1-CGYAYTCGCAY` (Leotiomycetes), q = 0.0111997
  - `2-AGWGCGAGWRCGAR` (Dothideomycetes) matches `3-TCGTACTCGCANYG` (Leotiomycetes), q = 1.37927e-06
  - `2-AGWGCGAGWRCGAR` (Dothideomycetes) matches `1-CGYAYTCGCAY` (Leotiomycetes), q = 9.79699e-05
- **Dothideomycetes vs Sordariomycetes: 7 matching motif pair(s)**
  - `1-CGARTACGCGTAY` (Dothideomycetes) matches `1-RATACGCGTATY` (Sordariomycetes), q = 9.82089e-05
  - `1-CGARTACGCGTAY` (Dothideomycetes) matches `2-WATACGMGKM` (Sordariomycetes), q = 0.000332401
  - `1-CGARTACGCGTAY` (Dothideomycetes) matches `5-ASCGAGACCGAGAMC` (Sordariomycetes), q = 0.0900042
  - `2-AGWGCGAGWRCGAR` (Dothideomycetes) matches `5-ASCGAGACCGAGAMC` (Sordariomycetes), q = 0.00855438
  - `3-AGCCAATGARA` (Dothideomycetes) matches `3-CCAATRAGAG` (Sordariomycetes), q = 0.00423513
- **Eurotiomycetes vs Leotiomycetes: 2 matching motif pair(s)**
  - `5-YTCGTACTCGTACY` (Eurotiomycetes) matches `3-TCGTACTCGCANYG` (Leotiomycetes), q = 4.05814e-05
  - `5-YTCGTACTCGTACY` (Eurotiomycetes) matches `1-CGYAYTCGCAY` (Leotiomycetes), q = 0.000268595
- **Eurotiomycetes vs Sordariomycetes: 1 matching motif pair(s)**
  - `5-YTCGTACTCGTACY` (Eurotiomycetes) matches `2-WATACGMGKM` (Sordariomycetes), q = 0.0192822
- **Leotiomycetes vs Sordariomycetes: 2 matching motif pair(s)**
  - `3-TCGTACTCGCANYG` (Leotiomycetes) matches `5-ASCGAGACCGAGAMC` (Sordariomycetes), q = 0.0311484
  - `3-TCGTACTCGCANYG` (Leotiomycetes) matches `2-WATACGMGKM` (Sordariomycetes), q = 0.0449285

**Clade IV: 5 motif(s) replicated across lineages, from 6 lineage pair(s).**

## Clade V

Tested in 6 lineages: Agaricomycetes, Dothideomycetes, Eurotiomycetes, Lecanoromycetes, Leotiomycetes, Sordariomycetes

- **Agaricomycetes vs Dothideomycetes: 1 matching motif pair(s)**
  - `8-AAGCTGTGGBTH` (Agaricomycetes) matches `4-AGCCGAGGCTY` (Dothideomycetes), q = 0.012158
- Agaricomycetes vs Eurotiomycetes: no matching motifs
- Agaricomycetes vs Lecanoromycetes: no matching motifs
- Agaricomycetes vs Leotiomycetes: no matching motifs
- **Agaricomycetes vs Sordariomycetes: 1 matching motif pair(s)**
  - `3-CACAGCCTC` (Agaricomycetes) matches `3-GAGGCTGTGC` (Sordariomycetes), q = 2.7638e-05
- Dothideomycetes vs Eurotiomycetes: no matching motifs
- Dothideomycetes vs Lecanoromycetes: no matching motifs
- Dothideomycetes vs Leotiomycetes: no matching motifs
- **Dothideomycetes vs Sordariomycetes: 2 matching motif pair(s)**
  - `1-CGACCGAGGCABCKT` (Dothideomycetes) matches `1-GCCGCGTCAYCRN` (Sordariomycetes), q = 0.0191466
  - `4-AGCCGAGGCTY` (Dothideomycetes) matches `5-WGCCGAGGMTCWCGG` (Sordariomycetes), q = 3.26993e-05
- Eurotiomycetes vs Lecanoromycetes: no matching motifs
- Eurotiomycetes vs Leotiomycetes: no matching motifs
- Eurotiomycetes vs Sordariomycetes: no matching motifs
- Lecanoromycetes vs Leotiomycetes: no matching motifs
- Lecanoromycetes vs Sordariomycetes: no matching motifs
- Leotiomycetes vs Sordariomycetes: no matching motifs

**Clade V: 4 motif(s) replicated across lineages, from 15 lineage pair(s).**


## Summary

| clade | lineages tested | lineage pairs | replicated motifs |
|---|---|---|---|
| I | 2 | 1 | **0** |
| II | 2 | 1 | **0** |
| III | 3 | 3 | **0** |
| IV | 4 | 6 | **5** |
| V | 6 | 15 | **4** |

**Total replicated motifs across all clades: 9.**
