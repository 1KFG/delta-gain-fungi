# DeltaGain_Fungi

**Project instance of [DeltaGain](https://github.com/1KFG/delta-gain), applied to the Fungi_BFD dataset.**

This repo holds the Fungi_BFD-specific run configuration, reports, and
(eventually) results for the DeltaGain analysis — quantifying how much new
protein/domain/fold diversity, taxonomic breadth, and gene-family expansion
the ~22k-genome Fungi_BFD set contributes beyond UniProt/RefSeq. It does not
contain the pipeline code itself.

## Split from DeltaGain (2026-09-07)

DeltaGain started as a single repo/folder combining pipeline code with
project-specific results and a reports website. It split in two, mirroring
how an nf-core-style pipeline separates from the project instances that run
it:

- **[`DeltaGain`](https://github.com/1KFG/delta-gain)** — the reusable
  pipeline: Nextflow workflows, modules, `bin/` scripts, profile configs.
- **`DeltaGain_Fungi`** (this repo) — `docs/` (results website), `reports/`
  (Quarto report source), and project-specific params/results.

**Known gap from the split, not yet resolved:** DeltaGain's `conf/` profiles
still hardcode Fungi_BFD-specific paths (`Fungi_BFD_runs/samples.csv`,
`nr_cluster_seq.dmnd`, etc.) directly in the pipeline repo, rather than this
repo supplying them via a `-params-file`. Until that refactor happens, this
repo doesn't yet drive its own runs — DeltaGain's launch scripts
(`run_genome_classify.sh`, `run_protein_novelty.sh`) still do, from within
the `DeltaGain` checkout, and `results/`/`work/` are still there too (a real
run was in flight at split time). Once that run completes, expect: those
results move here, and a `params-fungi.yaml` in this repo becomes the source
of truth for Fungi_BFD-specific paths, invoked as e.g.:

```
nextflow run 1KFG/delta-gain -r <tag> -profile genome_classify \
    -params-file params-fungi.yaml
```

## Website

Published at [dgf.fungalgenomes.org](https://dgf.fungalgenomes.org) via
GitHub Pages (`main`/`/docs`). `docs/index.md` is a placeholder landing page
deployed to confirm Pages/DNS/CNAME connectivity before any real report
content exists — see the design doc for the reporting plan.

`reports/` is the Quarto project source (`_quarto.yml` → `output-dir:
../docs/reports`); render with:

```
module load quarto
cd reports && quarto render
```

## Design doc

`/bigdata/stajichlab/shared/projects/BFD/Ideas/protein_novelty_esm2_plan.md`
— the living plan (shared with DeltaGain; not duplicated here).
