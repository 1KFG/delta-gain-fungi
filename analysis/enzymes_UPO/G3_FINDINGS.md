# G3 Phase A findings: is there clade-linked regulation of fungal UPOs?

**Measured 2026-11-08.** Design and parameters are in `PLAN.md` §4.2, all chosen by the user.
Primary output: `G3_MOTIF_REPLICATION.md`. Scripts: `upo_reg/bin/`.

## Design in one paragraph

7,055 independent units (one promoter per 90 % protein cluster, so *Fusarium*'s 5,103
near-identical copies cannot manufacture significance), 1,000 bp upstream truncated at the
neighbouring gene, against two controls: 5 genome-matched **and length-matched** non-UPO promoters
per genome, plus an order-2 shuffle. **The primary test is within-lineage clade contrasts** — a
motif counts only if it is significant on held-out data in two or more *independent lineages* for
the same clade. This avoids repeating the clade-versus-lineage confound already found once in this
project.

Extraction: **6,760 target and 15,659 background promoters**, zero unreadable genomes, 260 targets
dropped for under 300 bp of usable promoter.

## Significance filtering, which turned out to matter

STREME run with `--nmotifs 8` **reports 8 motifs per contrast whether or not they are
significant.** Across the 17 contrasts it reported 136 motifs; only **57 pass STREME's own
held-out test set at p < 0.05**. Several contrasts returned nothing significant at all
(Eurotiomycetes Clade V, Lecanoromycetes Clades III and V, Leotiomycetes Clade V). The replication
analysis counts only motifs that pass the held-out test in both lineages being compared.

Without that filter the analysis would have reported "8 motifs everywhere" and been meaningless.

## Result

| Clade | lineages tested | lineage pairs | **replicated motifs** |
|---|---|---|---|
| **IV** | 4 | 6 | **5** |
| **V** | 6 | 15 | **4** |
| III | 3 | 3 | 0 |
| I | 2 | 1 | 0 |
| II | 2 | 1 | 0 |

**Clades I, II and III show no cross-lineage replication.** For them this is a negative result.

### Clade IV replicates in all six lineage pairs

Clade IV (Dothideomycetes, Eurotiomycetes, Leotiomycetes, Sordariomycetes) is the only clade whose
motifs match in **every** lineage pair tested. Two distinct motif families:

**Family 1 — a `TACGCGTA` / `TCGTACTCG` element.** Appears as `CGARTACGCGTAY` (Dothideomycetes),
`RATACGCGTATY` and `WATACGMGKM` (Sordariomycetes), `YTCGTACTCGTACY` (Eurotiomycetes),
`CGYAYTCGCAY` and `TCGTACTCGCANYG` (Leotiomycetes). Best cross-lineage match q = 1.4e-06. **No
JASPAR match.** Per the standing caveat, unnamed means "no close match in a thin database", not
novel.

**Family 2 — a CCAAT box.** `AGCCAATGARA` (Dothideomycetes), `CCAATRAGAG` (Sordariomycetes),
`AGCTCATTGGGTCGS` (Eurotiomycetes, reverse strand). Cross-lineage q = 0.0042. **This one is
named**: Tomtom matches it to **JASPAR MA0314.3 = HAP3**, q = 0.023.

HAP3 is a subunit of the CCAAT-binding complex (Hap2/3/4/5 in *Saccharomyces*; HapB/HapC/HapE,
the AnCF complex, in filamentous fungi), which regulates heme biosynthesis, iron homeostasis and
the oxidative stress response. **A CCAAT-binding-complex element upstream of a heme-thiolate
peroxygenase family is mechanistically coherent** rather than arbitrary.

### Clade V replicates more weakly

Clade V replicates 4 motifs over 15 lineage pairs — a much lower rate than Clade IV's 5 over 6.
Its replicated motifs are GC-rich (`AGCCGAGGCTY`, `WGCCGAGGMTCWCGG`, `CACAGCCTC`) and none matched
JASPAR.

## What this does NOT establish

1. **CCAAT is the most common element in eukaryotic promoters.** Finding one is unsurprising on its
   own. The claim here rests entirely on it being *enriched against genome-matched, length-matched
   promoters from the same genomes* and passing a held-out test — not on its presence. Treat the
   enrichment, not the motif, as the finding.
2. **CCAAT-like motifs also appear in Clade V** (Dothideomycetes, Eurotiomycetes), they simply did
   not replicate across lineages there. So this is not a Clade IV-exclusive element.
3. **The HAP3 assignment is one Tomtom hit at q = 0.023 against a 178-motif, heavily
   *Saccharomyces* database.** It names a *family* of CCAAT-binding factors; it does not identify
   which fungal complex subunit, and it is not evidence that HAP3 itself binds these promoters.
4. **Everything about Clade IV inherits Clade IV's provisional status.** It has one characterized
   member (MphUPO), so "Clade IV" here means "proteins at high identity to one uncharacterized
   enzyme". See `G1_FINDINGS.md` Finding 3.
5. **Shared ancestry is not fully excluded.** Clade IV proteins sit at 66.5 % median identity, so
   some promoter similarity could be descent rather than convergent regulation. Cross-lineage
   replication argues against this — Sordariomycetes and Dothideomycetes are different classes, and
   promoters diverge fast — but no explicit test of promoter-level phylogenetic signal was run.
6. **No expression data.** Phase B (RNA-seq) is deferred at the user's direction. Nothing here
   shows these motifs are used, only that they are enriched.

## What would make this solid

- **Phase B**: test whether Clade IV UPOs are co-regulated with known CCAAT-complex targets in
  *Aspergillus* / *Penicillium* RNA-seq, and whether they respond to oxidative stress or iron.
- A promoter-phylogeny control for point 5.
- A wider fungal TF motif database than JASPAR CORE fungi, which would also let Family 1 be named
  or confidently called unknown.
