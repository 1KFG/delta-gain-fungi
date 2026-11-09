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
| **IV** | 4 | 6 | **5** — survives both controls |
| **V** | 6 | 15 | **4**, but see the control comparison below — **likely a composition artifact** |
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

## The two controls disagree, and the disagreement is the useful part

The pooled runs finished after the within-lineage analysis. Both report 10 significant motifs, but
they are **almost entirely different motifs**, and that is diagnostic.

| vs order-2 shuffle | vs genome-matched real background |
|---|---|
| `AAAAAAAAAAAAAAG` | `DATACGCGTATH` |
| `ATCATCATCATC` | `CCAATRRGARG` |
| `CTACCTACCTACY` | `WACGAGTACGAR` |
| `CCCYYCYYCYCCT` | `GTATAWAAGVYS` |
| `GCCGAGGCAS` | `AACGGCGCCGACYTT` |
| `CTCGTWCTCGTWCT` | `GATCTRCA` |

**The shuffle control is dominated by low-complexity and repeat tracts** — poly-A, `ATCATC`,
`CTACCTAC`, pyrimidine runs. A dinucleotide shuffle preserves dinucleotide frequency but destroys
repeats, so any repeat looks enriched against it. **None of these appear against the real
background.** This is why the genome-matched background was made the primary control, and it is
the reason the user's choice of two controls was worth the extra compute.

### This independently supports Clade IV

Both Clade IV motif families appear in the pooled run **against the real background**:

- `2-DATACGCGTATH` (p = 5.4e-07) and `5-WACGAGTACGAR` (p = 2.8e-05) — Family 1
- `3-CCAATRRGARG` (p = 1.5e-05) — Family 2, the CCAAT box

And the CCAAT motif names to the CCAAT-binding complex a **second** time, from a different run and
a different control: Tomtom matches `3-CCAATRRGARG` to **JASPAR MA0316.2 = HAP5**, q = 0.046.
HAP5 and HAP3 are different subunits of the same Hap2/3/4/5 complex. Two independent analyses
landing on two subunits of one complex is more convincing than either alone.

### And it weakens Clade V

Clade V's replicated motifs were GC-rich (`AGCCGAGGCTY`, `WGCCGAGGMTCWCGG`, `CACAGCCTC`). The
pooled shuffle run finds `3-GCCGAGGCAS` — the same element — **but the pooled real-background run
does not find it at all.**

**So the Clade V result is most likely a base-composition artifact rather than a regulatory
signal.** GC-rich promoter content differs between UPO genes and shuffled sequence without implying
a shared binding site. The Clade V row in the table above should be read as weak-to-negative, not
as a positive finding. Clade IV is the only clade whose motifs survive both controls.

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
5. **Shared ancestry is NOT excluded. It is now tested, and it is partly supported.**
   *Updated 2026-11-08; the earlier text said no explicit test had been run.*
   `upo_reg/bin/promoter_phylo_control.py` scans the 760 clade IV promoters with the actual
   STREME CCAAT PWM (`3-CCAATRAGAG`) under FIMO, then permutes carrier labels 10,000 times
   against the protein identity matrix.

   | quantity | value |
   |---|---|
   | carriers | 259 of 760 (34.1 %) |
   | mean pairwise identity, carriers | 57.49 % |
   | mean pairwise identity, all clade IV | 54.36 % |
   | difference | **+3.13 points** |
   | null mean (10,000 permutations) | −0.012 points (sd 0.914) |
   | one-sided p | **0.0005**, z = +3.44 |

   Carriers are more closely related to each other than a random subset of the same size.
   Within clade IV, which genes carry the motif tracks the phylogeny. **So the enrichment is
   at least partly confounded with ancestry and must not be read as convergent regulation.**

   This does not erase the cross-lineage replication: the motif still recurs in Sordariomycetes,
   Dothideomycetes and Eurotiomycetes, which are different classes. Both can hold if the motif is
   ancestral in part of clade IV and retained. The effect is also modest, 3 points on a 54 %
   baseline. But the direction is the one descent predicts, and it is significant.

   A methodological note worth keeping: a bare `CCAAT` 5-mer appears in 93.6 % of these promoters
   by chance and cannot discriminate anything. An earlier run of this control using the 5-mer gave
   a near-meaningless 711-vs-49 split. The test is only informative with the PWM.
6. **No expression data.** Phase B (RNA-seq) is deferred at the user's direction. Nothing here
   shows these motifs are used, only that they are enriched.

## What would make this solid

- **Phase B**: test whether Clade IV UPOs are co-regulated with known CCAAT-complex targets in
  *Aspergillus* / *Penicillium* RNA-seq, and whether they respond to oxidative stress or iron.
- A promoter-phylogeny control for point 5.
- A wider fungal TF motif database than JASPAR CORE fungi, which would also let Family 1 be named
  or confidently called unknown.
