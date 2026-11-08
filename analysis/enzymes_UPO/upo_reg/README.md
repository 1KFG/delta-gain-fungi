# G3: regulation of fungal UPOs

**Status: blocked on decision D6.** The argument contract for Phase A step 1 is settled
(`bin/extract_promoters.py`); the implementation is deliberately not written, because D6 changes
the scale by two orders of magnitude and therefore changes the design, not just a parameter.

## What D6 decides

| Option | Scale | What becomes possible |
|---|---|---|
| All-BFD | ~97,520 promoters across 15,137 genomes | A cross-fungal statement about UPO regulation. No expression validation. |
| RNA-seq-followable only | a few hundred *Penicillium* / *Aspergillus* promoters | A testable claim: predicted regulation checked against observed expression. |

These are different projects. The first is descriptive and large; the second is small and
falsifiable.

## Phase A design (sequence only)

1. **Promoters.** `--upstream` bp 5' of each CDS start, on the gene's own strand, truncated at the
   contig edge and at the 3' end of the nearest upstream gene so promoters never contain a
   neighbouring CDS. Genes with a usable promoter shorter than `--min-len` are reported and
   dropped, never padded.
2. **Background, matched by genome.** Promoters of randomly chosen non-UPO genes from the *same*
   genomes, fixed seed. This is the load-bearing control: across 15,000 fungal assemblies, GC
   content and gene-prediction style vary so much that an unmatched background would manufacture
   enrichment.
3. **Motif discovery** against that background, then
4. **the actual question: does motif content track clade?** Clade is the variable that carries the
   measured functional signal (`CALIBRATION_MATRIX.md` section 8), so "are UPO subtypes regulated
   differently" is the form of the question worth asking.

## Phase B

*Penicillium* and *Aspergillus* RNA-seq, to test whether predicted regulation matches observed
expression. Only meaningful under the second D6 option, or as a follow-up to it.

## Known risks before any result is reported

- **Promoter length is not uniform** across fungi; compact genomes have shorter intergenic regions.
  Any motif enrichment must be checked against promoter length as a confounder.
- **Clade assignment at BFD identity is an extrapolation** (see `../upo_sites/README.md`). Grouping
  promoters by clade inherits that uncertainty, so a clade-linked motif result is only as good as
  the clade calls underneath it.
- Gene models are funannotate predictions, not curated. A mispredicted start codon moves the
  promoter.
