---
title: DeltaGain
---

# DeltaGain

**Quantifying the annotation/novelty yield of large-scale fungal genome sequencing.**

How much new protein, domain and fold diversity, taxonomic breadth and gene-family
expansion the Fungi_BFD genome set contributes beyond what is already in UniProt and
RefSeq — and, alongside that, the analyses built on the same genome set.

## [Reports →](reports/)

Rendered analyses. Published so far:

- [**Delta Gain: findings synthesis**](reports/delta_gain/) — the novelty-yield measurement
  itself: whole-proteome funnel, what the survivors are, family-level diversity against an
  order-matched null, and species accumulation across 23,567 genomes.
- [**The enriched-Pfam signal is mostly Y1000+ deposition**](reports/y1000_deposition/) — a
  companion test showing that the funnel's largest family-level effect is mostly an artefact
  of where 1,033 genomes were deposited.
- [**Unspecific peroxygenases across 15,137 fungal genomes**](reports/upo_clades/) — a
  census of 97,520 PF01328 proteins, phylogenetic clade calling, candidates for
  expression, and a measured account of what the model cannot predict.
- [**Candidate effectors of two *Mucor* strains isolated from *Marchantia***](reports/mucor_secretome.html)
  — secretome, predicted folds, and a structure-based interaction screen.

The reports index lists these and what is deliberately not published yet.

## Code

- Analyses and pipeline code: [github.com/1KFG/delta-gain-fungi](https://github.com/1KFG/delta-gain-fungi)
- Original pipeline: [github.com/1KFG/delta-gain](https://github.com/1KFG/delta-gain)

Every result on this site is computational. Nothing has been tested at the bench, and
each report states its own limitations.
