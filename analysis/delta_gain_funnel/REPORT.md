# The enriched-Pfam signal is mostly Y1000+ deposition

2026-11-15. Code in `analysis/delta_gain_funnel/bin/y1000_join.py`, driver
`slurm/y1000_join.sh` (job 32102558, 1 min 20 s), tests in `test/` (17 passing). Per-family
results in `out/y1000_enriched.tsv` and `out/y1000_control.tsv`.

## The question

The delta-gain funnel's largest family-level result is that Pfam families chosen for being
enriched in new genomes carry **+14.6 to +14.9 points** more novel clusters than an
order-matched baseline null, against **+2.2** for size-matched random Pfam controls. The
funnel's own notes observe that the enriched list "is a Saccharomycotina yeast list" and
read this as a lineage effect, supported by the fact that it replicates when the list is
selected on one half of the species and evaluated on the other.

The objection the split-sample test does not answer: **replication across species does not
distinguish new biology from a shared deposition history.** "Baseline" in the funnel means
"had a prior NCBI annotation". That is a boundary around a database, not around knowledge.
A sequencing project that deposits assemblies to GenBank unannotated and releases its own
proteomes elsewhere puts every one of its genomes on the new side, however long its
proteins have been public — and it does so for all of its species at once, so a split-sample
test replicates it perfectly.

## The join

Y1000+ is such a project, and it is in BFD at scale.

| BioProject | NCBI title | genomes in BFD |
|---|---|---|
| `PRJNA736342` | Y1000+ Saccharomycotina yeasts genome sequencing and assembly (UW-Madison) | 914, across 883 species |
| `PRJNA429441` | Sequencing of 196 yeast species from the subphylum Saccharomycotina (UW-Madison) | 116 |

**1,033 genomes in total, every one of them on the funnel's new side.** Opulente et al.
(2024) state it plainly in their data-availability section: assemblies and raw reads to
GenBank, everything else to Figshare. GenBank therefore holds these assemblies without an
NCBI annotation, BFD re-annotated them, and their proteins are absent from nr — while the
project's own annotations have been public since 2024 and InterProScan was run across all
1,154 of them.

The test re-asks the funnel's own question with those genomes moved across the line. A
cluster is novel when no member comes from a baseline genome; the second column asks how
many clusters have no member from a baseline **or Y1000+** genome. Same rule, larger
known set.

## Result

At c70, the tier the funnel reports:

| family set | clusters | novel now | novel if Y1000+ counted as known | share of novelty lost |
|---|---|---|---|---|
| **Enriched Pfam (204)** | 108,774 | 77,206 (71.0 %) | 19,010 (**17.5 %**) | **75.4 %** |
| Random-control Pfam (197) | 51,479 | 23,796 (46.2 %) | 20,879 (**40.6 %**) | 12.3 % |

Every tier behaves the same way:

| tier | enriched: novel → if known (lost) | control: novel → if known (lost) |
|---|---|---|
| c90 | 71.4 % → 26.5 % (62.9 %) | 51.1 % → 46.0 % (9.9 %) |
| c70 | 71.0 % → 17.5 % (75.4 %) | 46.2 % → 40.6 % (12.3 %) |
| c50 | 65.6 % → 17.0 % (74.1 %) | 43.9 % → 39.3 % (10.4 %) |
| c30 | 57.8 % → 21.4 % (62.9 %) | 42.9 % → 39.2 % (8.6 %) |

Per family at c70, the two sets barely overlap:

| | enriched (204) | control (197) |
|---|---|---|
| median share of novelty lost | **86.6 %** | **0.8 %** |
| families losing > 50 % | 194 | 28 |
| families losing > 90 % | 23 | 0 |
| families losing < 10 % | 2 | 130 |

**The control is what makes this an argument.** Y1000+ contributes ~883 distinct species at
roughly one genome each, so it will add singleton clusters to any family; if that were the
whole story the control families would shrink too. They lose 12.3 % against the enriched
list's 75.4 %, and their median family loses 0.8 %.

**The excess does not merely shrink — it inverts.** Before the join the enriched families
are 71.0 % novel against the control's 46.2 %, a 24.8-point excess. After it they are
17.5 % against 40.6 %: **23.1 points below the controls.** Whatever "enriched in new
genomes" selected for, with Y1000+ counted as public it is the opposite of a novelty-rich
set. 1,033 genomes, 5.8 % of the 17,827 new-side genomes in these families, carry three
quarters of their novelty.

## What this does and does not establish

- **It is not a recomputation of the +14.6-point figure.** That statistic is a comparison
  against an order-matched baseline null with its own resampling; this is the raw
  novel-cluster share under two definitions of "known". The direction and size are not in
  doubt, but the published number should be recomputed properly by re-running the
  order-matched null with Y1000+ in the baseline, not adjusted by hand from this table.
- **It covers 401 Pfam families, not the corpus.** The funnel's genome-wide figures
  (50.2 / 29.2 / 22.4 % novel at 95/70/50 %) are untouched by this run. Y1000+ is 5.8 % of
  the new-side genomes here, and its species are nearly all distinct, so some effect on the
  genome-wide numbers is likely — it has not been measured.
- **It does not show the enriched families are uninteresting.** Yeast-specific families that
  Y1000+ sampled deeply are real biology; the claim it undercuts is that the funnel's
  database comparison *discovered* them.
- **Y1000+ is the largest such project, not the only one.** Any project depositing
  assemblies without annotation has the same effect. The same join against other large
  unannotated submissions would say whether this is one project or a general pattern, and
  that has not been run.

## Consequence for the funnel

The funnel's limitations already state that baseline is "a sampling, not a knowledge
boundary" and that the enriched list "is selected on the outcome; it is a candidate list,
not a genome-wide estimate." Both are correct and both are now quantified: for this result
the sampling boundary accounts for three quarters of it.

The honest summary of the enriched-Pfam result is that it measures **where Saccharomycotina
genomes were deposited**, and that the lineage reading is not supported once deposition is
controlled for. It should not be published as a novelty finding in its present form.

## Upstream issue found

`SUBPHYLUM` is **NULL for all 23,680 rows** of both `samples.parquet` and `species.parquet`
in `Fungi_BFD_runs/tables/`. `CLASS` (23,462 populated) and `ORDER` (23,380) are fine, so
subphylum-level grouping has to be reconstructed from the class list — this analysis used
the nine Saccharomycotina classes present in BFD rather than the subphylum field. Worth
fixing where the table is built; it is the second such column found after
`is_gpi_anchored`, which is TRUE for every protein.
