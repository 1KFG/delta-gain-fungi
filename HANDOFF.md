# DeltaGain_Fungi — Handoff / Status (2026-09-07)

## Where things actually stand

**Stage 1 (genome classification, `DeltaGain/genome_classify.nf`):** code complete,
validated end-to-end against real data at small scale (`--n_test 2`, real
genomes from `Fungi_BFD_runs/samples.csv`) — not yet run at full 22k-genome
scale. Three real bugs found and fixed during that validation (see "Bugs
found" below); the pipeline had never previously been exercised against
anything but tiny synthetic stub data, so this was the first real test.

**Stage 3 (protein novelty, `DeltaGain/protein_novelty.nf`):** the UniProt
fungi reference-set half is **done and real**:
- `results/protein_novelty/refdb/uniprot_fungi_reviewed.fasta` — small REST
  smoke-test output (Swiss-Prot only, 20 records; not meant to be used for
  anything beyond confirming connectivity)
- `results/protein_novelty/uniprot_fungi_nr.fasta` — **15,089,757
  non-redundant representative sequences**, MMseqs2 linclust at 95%
  identity / 90% coverage, from the full comprehensive UniProtKB fungi set
  (Swiss-Prot 38,451 + TrEMBL 17,534,069 = 17,572,520 — matches the live
  REST total-count check exactly)
- `results/protein_novelty/uniprot_fungi_clusters.tsv` — cluster membership
  (representative accession → member accession), joins back to the
  InterPro/Pfam/GO xrefs table produced alongside the FASTA conversion

The BFD-side of Stage 3 (extracting BFD's own ~176M proteins, dereplicating
them, and running the actual novelty search against both this UniProt-nr set
and the already-staged `nr_cluster_seq.dmnd`) is **not started** — see
"Next steps" below, item 1, which is the real blocker.

## What's next (priority order)

1. **Locate/confirm the BFD annotated-protein source.** Every downstream
   novelty-screen step depends on this and nothing else can start until it's
   settled. Likely `Fungi_BFD/nextflow/input/pep/` per naming convention seen
   elsewhere, but not confirmed for the full ~22k-genome set — some genomes
   may not be annotated yet at all (which is itself a relevant fact: no
   annotation yet = can't screen its proteins for novelty yet, a separate
   bookkeeping concern from the analysis itself). **This is the one thing
   most worth resolving first**, and it's mostly a lookup/confirmation task,
   not new code.

2. **Run Stage 1 at meaningfully larger scale before the full 22k.** Try
   `--n_test 50` or `100` next, not straight to the full set — confirms
   `MASH_PREFILTER_GLOBAL`'s all-vs-all step and memory scaling hold up
   before committing a `highmem`/`epyc` allocation to the real ~22k run.
   Real thresholds (95%/90% ANI, 80% mash prefilter floor) are still the
   simplifying defaults flagged in the design doc — worth a scientific
   sanity-check on a real, larger result before trusting them at full scale.

3. **Wire `DIAMOND_NOVELTY_SCREEN` into a workflow** — module already
   written and correct (confirmed `--taxonlist` approach works), just not
   connected to any `.nf` file yet. Needs: BFD protein extraction (item 1) →
   MMseqs2 linclust dereplication of BFD's own proteins (same tool, same
   linclust settings already proven at UniProt-fungi scale) → DIAMOND search
   against `nr_cluster_seq.dmnd` and/or `uniprot_fungi_nr.fasta`.

4. **E-value/coverage binning** (`bin/bin_novelty_hits.py`, not yet written)
   — turns raw DIAMOND hits into strong-hit / weak-hit / no-hit per the
   design doc's framing (E-value + coverage, not raw %identity).

5. **Stage 3b controls** (artefact/contamination/fast-evolving-family
   filters) and **Stage 5 domain/fold sorting** (Pfam/InterPro + Foldseek,
   sorting each candidate-novel protein into one of the five novelty axes)
   — both still design-doc-only, no code yet.

6. **Reports website.** `docs/`/`reports/` are still a placeholder landing
   page — once Stage 1 has a real (even partial) result, the first real
   report page (genome classification counts by `cluster_class`) is
   low-effort to add and would make the site actually useful rather than
   just a connectivity test.

## What you can do in parallel right now (doesn't block on my end)

- **Track down the BFD protein FASTA location** (item 1 above) — this is
  the single highest-leverage thing you could resolve independently, since
  it's a "where does X live" question about `Fungi_BFD`/`Fungi_BFD_runs`
  that you likely know or can find faster than I can by searching.
- **Sanity-check the ANI/mash thresholds** (95% cluster / 90% outlier /
  80% mash prefilter floor) against your own knowledge of fungal species
  boundaries — these are flagged as simplifications in the design doc and
  are a scientific judgment call, not something to leave to a default.
- **Decide the E-value/coverage cutoffs** for the Stage 3 novelty screen
  (item 4) — same kind of call, better made deliberately than defaulted.
- **Skim `results/protein_novelty/uniprot_fungi_clusters.tsv`** if you want
  a sense of what the real UniProt-fungi redundancy structure looks like
  before we build on top of it (17.5M → 15.1M is a smaller collapse than
  the design doc originally guessed — worth knowing why: TrEMBL fungi may
  be less redundant per-protein than assumed, or the 95%/90% thresholds are
  conservative).

## Bugs found during this session's validation (for future reference)

- A `-c`-supplied Nextflow config's `params{}` block can't see a `params.x`
  value set by the auto-loaded root `nextflow.config` via string
  interpolation (general Nextflow limitation, reproduced in isolation) — fix
  was having each project conf resolve shared values (like
  `singularity_cache`) independently rather than depending on cross-file
  propagation.
- Calling `.countLines()` on a process input `Path` inside a resource
  directive (`cpus`/`memory`/`time`) throws `ProviderMismatchException` —
  fix was passing the count as an explicit `val()` computed in Groovy at
  channel-construction time (same convention `Fungi_BFD/nextflow`'s
  `SKANI_COMPARE` already uses).
- `highmem`'s nodes are mixed hardware — 6 of 7 predate AVX2 (Ivy Bridge),
  only `h07` is modern. Any AVX2-requiring container should request `epyc`
  instead (20 nodes, uniformly modern AMD Milan, ~950GB RAM each) rather
  than `highmem`, unless a job genuinely needs more than epyc's per-node
  memory ceiling.
- `MMSEQS_LINCLUST` was missing a `publishDir` — real output existed only in
  `work/` until caught after the run completed. Fixed; worth double-checking
  any other module for the same gap before the next real run.
