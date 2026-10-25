# DeltaGain_Fungi — Handoff / Status (2026-10-24)

## Update 2026-10-24: 2000-genome protein screen (read this first)

The sections below are the 2026-09-08 pilot record. They are still correct
for the pilot. This section replaces them where they disagree.

**State of the input data**
- `Fungi_BFD_runs/input/pep/` now holds 23,448 files named
  `*.proteins.fa.gz` (gzipped). The pilot read `*.proteins.fa`. The old
  "46% annotated" figure is out of date: 23,447 of 23,680 `samples.csv`
  rows now match a pep file (229 not annotated, 4 ambiguous tags).
- The BFD annotation run is still in progress (about 10 more days of
  compute, per the user). The pep files may change again. If they change
  after a run starts, that run's manifest no longer matches the files.
- Final genome tables will land in `Fungi_BFD_runs/tables/` as parquet when
  the full run ends. The current contents are a snapshot from
  2026-07 to 2026-09-10. Row counts and genome coverage of that snapshot are
  not verified (no `duckdb` in the default python). Do not use it for the
  accumulation curves.

**Code changes (DeltaGain repo, all committed)**
- `bin/build_pep_manifest.py`: matches `<tag>.proteins.fa.gz`, falls back to
  `.fa`. Uses one `os.listdir` instead of about 23k `os.path.exists` calls
  (the per-file version ran over 7 minutes on GPFS; the listing takes 0.4 s).
  `--n-test N` now keeps a **seeded random sample** (`--seed`, param
  `n_test_seed`, default 42) of all matched rows. It used to keep the first
  N rows in `samples.csv` order. Same seed and input give an identical
  manifest (tested); seed 7 gives a different one.
- `bin/combine_bfd_proteins.py`: reads `.gz` input.
- New `modules/compress/COMPRESS_ZSTD`: publishes zstd-compressed copies
  (level 6) of the combined FASTA, provenance TSV, both cluster TSVs and
  both representative FASTAs. Producers no longer publish the plain BFD
  files. Downstream steps still read the plain files from `work/`.
  - No cached container has a `zstd` binary (checked python_3.12, diamond,
    seqkit, mmseqs2). The process runs on the host with the shared install at
    `params.zstd_bin`
    (`/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin/zstd`).
    If that binary is missing it falls back to `pigz`, then `gzip`
    (output suffix `.gz`).
  - Level choice: on 200 MB of pilot FASTA, levels 3/6/10 gave
    98.4/100.4/96.1 MB. Protein text gains little from higher levels. Level 19
    was killed on the login node (signal 9), so it is not measured.
  - **Bug found on the first real run:** Nextflow stages inputs as symlinks,
    and `zstd` silently skips symlinked input files. The first version made
    no output and failed at `zstd -t`. Fixed by reading from stdin
    (`zstd ... -o out < in`). A `-stub-run` cannot catch this. The stub test
    only touches the output file.
- `diamond_chunk_size` raised 20,000 -> 60,000. Basis: 24 successful pilot
  DIAMOND tasks at 20,000 queries took 1,100-3,000 s each (from work-dir
  `.command.begin` to `.exitcode`). Three times that targets 1-1.5 h per
  task. This assumes runtime is linear in query count. DB load time is a fixed
  cost, so the real runtime may be lower. Check the first completed tasks.

**Run in progress**
- Config: `n_test = 2000` (about 21M raw proteins, extrapolated from the
  pilot's 10.7k proteins per genome), seed 42, output to
  `results/protein_novelty_n2000/` (the pilot's `results/protein_novelty/`
  is untouched).
- Job 30452834 (second submission). The first, 30451574, was cancelled after
  26 min because of the zstd symlink bug above. `-resume` reuses the cached
  manifest and combine steps.
- **Measured so far (2026-10-24, job 30452834):**
  - Raw proteins: 20,373,920 from 2000 genomes (10.2k per genome; the
    estimate was 21M). After linclust (95% id / 90% cov): 10,615,362
    representatives, which is 52.1% of raw. The pilot kept 62.4%
    (666,654 / 1,067,931), so the ratio changed with scale or sample.
  - `MMSEQS_LINCLUST_BFD`: 5 min 34 s, MaxRSS 25.6 GB. The request was
    32 cpu / 256 GB / 24 h, which is about 10x too much memory. For the full
    set (about 10x the proteins), 32 cpu / 64-128 GB is a reasonable next
    request. That is an extrapolation, not a measurement.
  - `COMBINE_BFD_PROTEINS`: 8 min 23 s, MaxRSS 107 MB (single-threaded
    python, gzip input).
  - `COMPRESS_ZSTD`: 8 cpu; the 20M-protein combined FASTA took 68 s. All
    published `.zst` files pass `zstd -t`. Sizes: combined FASTA 5.46 GB,
    representatives 2.81 GB, clusters 158 MB, provenance 28 MB.
  - `SPLIT_FASTA`: 177 chunks of 60,000. First 3 `DIAMOND_SEARCH_UNIPROT`
    tasks finished in 860 s on average (max 1,029 s), well under the 1-1.5 h
    target. `DIAMOND_NOVELTY_SCREEN` (the larger `nr_cluster` database) had
    not finished any task when this was written. If it also finishes in
    under about 45 min, raise the chunk size again for the full run.
- **`DIAMOND_CLUSTER` lost-task problem (2026-10-24):** Nextflow marked
  `DIAMOND_CLUSTER` as failed ("terminated for an unknown reason") after the
  20-minute `exitReadTimeout`, while SLURM job 30452869 was still running on
  hz02. Nextflow started a retry (30454465) that waits in the queue
  (`AssocGrpCpuLimit`: the 170+ pending DIAMOND chunk jobs fill the account
  CPU limit). This is the same exit-status problem as in the pilot, now at
  about 180 concurrent jobs plus 500+ other jobs from the BFD run. The branch
  is comparison-only, so the novelty bins do not depend on it. Two copies may
  run and write to different work dirs. Not fixed. Options: raise
  `exitReadTimeout`, or lower concurrent chunk jobs with `queueSize`.
- **Cause of the lost tasks and duplicate jobs (found 2026-10-24):**
  the chunk jobs were submitted to `epyc` but were moved to other
  partitions (`stajichlab`, `preempt`, `highclock`). Nextflow polled only
  the submit partition, lost sight of the moved jobs, and resubmitted them.
  Result: 106 retries ("terminated for an unknown reason") and 104 extra
  chunk jobs (458 submitted for 354 needed). Some retries duplicate work
  that already finished. Fix: `executor.queueGlobalStatus = true` in
  `DeltaGain/nextflow.config` (same setting as `Fungi_BFD/nextflow`). It
  applies from the next launch, not to the run in progress.
- 257 pending DIAMOND jobs were put on hold with `scontrol hold` at the
  user's request (running jobs left alone, head job 30452834 still alive).
  Release with `scontrol release <jobid>`. I did not verify that Nextflow
  tolerates a long hold.
- `DIAMOND_CLUSTER` resources (32 cpu / 128 GB / 48 h) are still a guess.
- The accumulation-curve branch is skipped. The four `bfd_duckdb_*` /
  `genome_classification_tsv` / `bfd_version_tag` params are unset on
  purpose, until the final tables exist.
- Plain copies of the large files stay in `work/protein_novelty/`. Run
  `nextflow clean` after the run to free the space.

**Functional data needed before any "novel" call**
A no-hit protein is only a candidate. To call it novel, these are needed:
1. `asm_stats` and `busco_genome` parquet (final), for the accumulation
   curves and genome QC.
2. Per-protein Pfam/InterPro (Stage 5 input). `tables/pfam.parquet` exists
   as a snapshot.
3. Signal peptide, TM, TargetP and secondary-metabolite calls, to separate
   real secreted novel proteins from artefacts. Snapshots exist in `tables/`.
4. A repeat/TE overlap flag per gene. Not found in `tables/`. EarlGrey runs
   are in `Fungi_BFD_runs/do_earlgrey` (not inspected).
5. A contamination/taxonomy check per genome. Not found in `tables/`.
Items 4 and 5 are open. Without them, no-hit proteins can include
transposon-derived genes and contaminant sequence.

**Next steps**
1. When the run ends: check the `.zst` outputs exist and pass `zstd -t`,
   read task resources from `sacct`, and compute the strong/weak/no-hit
   fractions at 2000 genomes. Compare them to the pilot (71.0 / 22.6 / 6.5%).
   Do not assume the pilot fractions hold.
2. After the full BFD run: wire the final parquet paths and
   `bfd_version_tag`, run the accumulation curves, and tune the E-value and
   coverage cutoffs on the final data.
3. Decide the full-run chunk size and resources from step 1's measurements.
   At 20,000 per chunk the full set would make about 12,500 chunks, and at
   60,000 about 2,600 per search (both extrapolations).
4. The older items below (ARI 0.72 question, Stage 3b controls, Stage 5, the
   reports page, the DIAMOND `--reps`+`--out` bug report) are still open. The
   bug-report repro script is in session scratch and may be gone.

## Where things actually stand (2026-09-08 pilot record)

**Stage 1 (genome classification, `DeltaGain/genome_classify.nf`):** code
complete. **Real `--n_test 100` pilot completed successfully 2026-09-08**
(job 28200405, 121/121 tasks). Results in
`results/genome_classify/genome_classification/genome_classification.tsv`:

| cluster_class | count |
|---|---|
| redundant_strain | 53 |
| singleton_isolated | 21 |
| unresolved_below_ani_floor | 14 |
| label_mismatch | 11 |
| novel_strain_same_species | 1 |

ANI spread across the 68 within-component pairwise comparisons: min 87.49%,
median 99.13%, max 99.99%. First empirical read: most redundant-strain calls
sit far above the 95% threshold (median 99.13%), and the one genuinely
borderline pair (87.49%) correctly fell into `unresolved_below_ani_floor`
rather than being forced into the wrong bucket — a good first signal the
95%/90%/80% thresholds are reasonable, though n=100 is still small for a
confident call at full 22k scale. The 11% `label_mismatch` rate is itself a
real finding worth a look (potential mislabeled genomes/metadata QC, per the
design doc's expectation for this step).

**Two real bugs found and fixed getting this pilot to run** (neither
surfaced at the earlier `--n_test 2` validation, which only ever produced
singleton components and so never exercised SKANI_TRIANGLE or the
multi-genome path):
1. `run_genome_classify.sh`/`run_protein_novelty.sh` had
   `module load singularity 2>/dev/null || module load apptainer ...` —
   loading `singularity` succeeds without providing the `apptainer` binary
   the config actually needs (`apptainer.enabled=true`), so the fallback
   never fired and every task failed with "env: 'apptainer': No such file
   or directory." Fixed: load `apptainer` first in both scripts.
2. `bin/find_components.py` had two bugs from the same root cause (`mash
   sketch`'s internal sketch ID is the input FASTA filename, not the `-o`
   prefix, so `mash dist`'s query/ref columns are `<ASMID>.fa.gz`, not
   `<ASMID>.msh`): (a) stripping only a literal `.msh` suffix left `.fa.gz`
   in place, doubling it in downstream genome paths
   (`.../<ASMID>.fa.gz.fa.gz`) and making skani skip every genome in every
   multi-member component; (b) `all_asmids` was seeded from the full
   samples.csv (23,683 rows) rather than only genomes that actually
   appeared in the mash-dist input, so every genome NOT part of an
   `--n_test`-limited run showed up as a false `singleton_isolated` (100
   real genomes analyzed came out as 23,604 "singletons" before the fix).
   Both fixed in `find_components.py`; see its inline comments for detail.

**Stage 3 (protein novelty, `DeltaGain/protein_novelty.nf`):** BFD-side
pilot **completed for real 2026-09-08** (job 28206550, COMPLETED/exit 0,
3h54m). Results in `results/protein_novelty/`:

- **UniProt-fungi reference side**: real dereplicated set already built —
  `DeltaGain/results/protein_novelty/uniprot_fungi_nr.fasta` (15,089,757
  sequences, OLD pre-split location). `params.uniprot_fungi_nr_fasta` now
  makes the workflow reuse this file directly by default instead of
  re-deriving it (13GB FTP download + ~17.5M-sequence linclust would
  otherwise be redone from scratch every run for no reason) — set to `''`
  explicitly to force a real re-derivation. `DIAMOND_MAKEDB` builds a
  DIAMOND db from it for the BFD-side search.
- **BFD-side** (wired 2026-09-08): `BUILD_PEP_MANIFEST` maps ASMID -> pep
  FASTA via the exact `makeSampleTag` convention
  `Fungi_BFD/nextflow/modules/common/utils.nf` uses (ported to Python, not
  guessed) -> `COMBINE_BFD_PROTEINS` (prefixes IDs with ASMID for
  provenance) -> `MMSEQS_LINCLUST` (reused, generalized to take a `prefix`)
  -> `SPLIT_FASTA` -> `DIAMOND_NOVELTY_SCREEN` (vs. `nr_cluster_seq.dmnd`)
  + `DIAMOND_SEARCH_UNIPROT` (vs. the UniProt-fungi db above) ->
  `BIN_NOVELTY_HITS` (strong/weak/no-hit per E-value + coverage).
- **Annotation coverage confirmed 2026-09-08:** 10,922/23,683 samples.csv
  rows (46%) have a matching pep file in
  `Fungi_BFD_runs/input/pep/` today; the rest aren't annotated yet. 6
  species+strain tags collide across 2 ASMIDs each (flagged
  `ambiguous_tag`, not guessed at) — see `bfd_pep_manifest.tsv`'s output
  each run for the exact list.

**Real Stage 3 pilot results (100 genomes, 1,067,931 raw BFD proteins):**

`bfd_novelty_bins.tsv` (666,655 dereplicated proteins screened against both
`nr_cluster_seq.dmnd` and UniProt-fungi):

| status | count | % |
|---|---|---|
| strong_hit | 473,000 | 71.0% |
| weak_hit | 150,648 | 22.6% |
| no_hit | 43,006 | 6.5% |

`bfd_clustering_comparison.tsv` (MMseqs2 linclust vs. DIAMOND cluster,
same 1,067,931-protein input, both at nominal 95%/90% thresholds — see
"DIAMOND DeepClust methodology question" below):

| metric | MMseqs2 linclust | DIAMOND cluster |
|---|---|---|
| n_clusters | 666,654 | 776,574 |
| n_singletons | 551,469 | 698,701 |
| mean cluster size | 1.60 | 1.38 |
| max cluster size | 136 | 90 |

**Adjusted Rand Index = 0.7155** — moderate-to-good agreement but clearly
not identical; DIAMOND cluster is noticeably more conservative (more/smaller
clusters) at the same nominal thresholds. Open question for whoever picks
this up next: is this a real algorithmic difference, or do DIAMOND's
`--approx-id`/`--mutual-cover` semantics just not map onto MMseqs2's
`--min-seq-id`/`-c --cov-mode 0` as directly as assumed? Worth a closer look
before trusting either number as "the" dereplication result at full scale.

**DIAMOND DeepClust methodology question (raised by user 2026-09-08):**
diamond/2.2.6 (module now installed; container built at
`/bigdata/stajichlab/shared/singularity_cache/diamond_2.2.6--he361c42_0.sif`
from `quay.io/biocontainers/diamond:2.2.6--he361c42_0`) already ships the
cascaded-clustering primitives behind DIAMOND DeepClust (Buchfink et al.
2026 Nat Methods) as `diamond cluster --cluster-steps/--round-coverage`
plus `realign`/`recluster`/`reassign` -- no separate "deepclust" binary
needed. Wired as `DIAMOND_CLUSTER` + `EXTRACT_CLUSTER_REPS` +
`COMPARE_CLUSTERINGS`, run alongside (not replacing) `MMSEQS_LINCLUST_BFD`.
Not yet applied to the UniProt-fungi reference set itself -- deliberately
deferred pending this same comparison's outcome (see design-doc decision
above); the ~0.72 ARI is a first signal that it might be worth extending,
but not yet conclusive at n=100.

## What's next (priority order)

1. **Both Stage 1 and Stage 3 pilots are done for real (see results above)
   — decide next scale-up step.** n=100 gave clean, sane-looking results on
   both; consider a bigger intermediate run (500-1000 genomes?) before the
   full ~22k / ~10.9k currently-annotated set, both to get a more solid read
   on the ANI-threshold and clustering-method questions and to catch any
   scale-dependent issues 100 genomes/1M proteins wouldn't surface (e.g. the
   ~34-concurrent-DIAMOND-task GPFS exit-status-timeout retries hit at this
   pilot's scale — see Bugs below — will only get more frequent at 10x-100x
   the concurrency).

2. **Decide on the MMseqs2-linclust-vs-DIAMOND-cluster question** (ARI=0.72
   at pilot scale) before scaling either dereplication step further — pin
   down whether the disagreement is a real algorithmic difference or a
   parameter-equivalence mismatch (`--approx-id`/`--mutual-cover` vs.
   `--min-seq-id`/`-c --cov-mode 0`) before trusting either number at full
   scale, and decide whether to extend the comparison to the UniProt-fungi
   reference set too.

3. **Revisit the E-value/coverage thresholds in `bin/bin_novelty_hits.py`**
   now that a real hit-score distribution exists (71% strong / 22.6% weak /
   6.5% no-hit at pilot scale) — the current `evalue-strong=1e-5`/
   `cov-strong=70%` defaults are still placeholders, not derived from this
   dataset; the real distribution is now available to inform that call.

4. **Stage 3b controls** (artefact/contamination/fast-evolving-family
   filters) and **Stage 5 domain/fold sorting** (Pfam/InterPro + Foldseek) —
   both still design-doc-only, no code yet. The 6.5% real no-hit bucket
   (43,006 proteins) is the concrete candidate-novel set these controls
   would filter, once written.

5. **Reports website.** `docs/`/`reports/` are still a placeholder landing
   page — real results now exist for both Stage 1 (`cluster_class` counts)
   and Stage 3 (novelty bins, clustering comparison); a first real report
   page is low-effort to add.

6. **File the DIAMOND `--reps`+`--out` bug report** — repro script at
   `/scratch/jstajich/28198688/claude-1181/-bigdata-stajichlab-shared-projects-BFD-DeltaGain-Fungi/4ef395c8-b544-4707-a44e-bc7492252cde/scratchpad/diamond_cluster_reps_out_bug_repro.sh`
   (session-scratch, copy it somewhere durable first) — confirmed
   reproducible on diamond 2.2.6.186, both containerized and native.

## Open decisions still needing a human judgment call

- **ANI thresholds** (95% cluster / 90% outlier / 80% mash prefilter floor)
  — use as-is for the pilot (per 2026-09-08 direction), revisit empirically
  once the pilot's real ANI distribution is in hand (item 1 above).
- **E-value/coverage cutoffs** for the novelty screen — initial defaults are
  now written into `bin/bin_novelty_hits.py` (1e-5 / 70% coverage), flagged
  as placeholders pending real hit-score data (item 3 above).

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
- The Fungi_BFD_runs `input/pep/*.proteins.fa` filenames are NOT ASMID-keyed
  — they're built from `makeSampleTag(SPECIES, STRAIN)` (note: `SPECIES`
  column, not `SPECIES_IN` — the two differ on ~76/2000 sampled rows, mostly
  `sp.`/`aff.` qualifiers). Reusing the exact Groovy sanitization rule from
  `Fungi_BFD/nextflow/modules/common/utils.nf:58-114` (ported to Python in
  `bin/build_pep_manifest.py`) was necessary — a naive re-derivation would
  have missed the asterisk-handling and character-collapsing edge cases.
- **DIAMOND 2.2.6 bug**: `diamond cluster --out X --reps Y` (both flags in
  the same invocation) crashes ("Error opening file .../input.tsv") right
  after computing the clustering, even on trivial 4-sequence inputs;
  `--out` alone works fine. Confirmed 3x (container, native module, tiny and
  pilot-scale real data). Worked around by extracting representative
  sequences ourselves from `--out`'s table (`bin/extract_cluster_reps.py`,
  now its own process, `EXTRACT_CLUSTER_REPS`). Standalone repro script at
  `/scratch/jstajich/28198688/claude-1181/-bigdata-stajichlab-shared-projects-BFD-DeltaGain-Fungi/4ef395c8-b544-4707-a44e-bc7492252cde/scratchpad/diamond_cluster_reps_out_bug_repro.sh`
  (session-scratch — copy elsewhere before it's cleaned up) if this ever
  gets reported upstream.
- The `quay.io/biocontainers/diamond` image has **no python3 at all** — any
  Python postprocessing step run inside a DIAMOND-labeled process fails with
  `env: can't execute 'python3'`. This is why `EXTRACT_CLUSTER_REPS` is a
  separate process (python container) rather than tacked onto the end of
  `DIAMOND_CLUSTER`'s script.
- DIAMOND's sensitivity modes (`--sensitive`/`--more-sensitive`/
  `--very-sensitive`/`--ultra-sensitive`) are standalone boolean flags, NOT
  a `--sensitivity <value>` option — `DIAMOND_NOVELTY_SCREEN` and
  `DIAMOND_SEARCH_UNIPROT` both had this wrong since they were first
  written (never actually run for real until 2026-09-08's Stage 3 pilot;
  "Error: Invalid option: sensitivity" on every chunk). Fixed:
  `--${params.diamond_sensitivity}` builds the flag name directly.
- At ~30+ concurrent DIAMOND search tasks (SPLIT_FASTA's per-pilot chunk
  count), several genuinely completed (real `exit 0`, correct output) but
  Nextflow still logged them `[ERROR]`/"terminated for an unknown reason":
  `nextflow.config`'s `exitReadTimeout = '20 min'` expired before Nextflow
  could see the SLURM exit-status file — almost certainly GPFS
  metadata-visibility lag under heavy concurrent load, not a real task
  failure. Nextflow's automatic retry handled every instance correctly
  (retry succeeded), but the run's own interim `[FAILED]` summary line looks
  alarming until you check `sacct`/the actual output files. Not fixed (no
  code change needed — retry already covers it), but worth knowing before
  scaling past ~30 concurrent DIAMOND chunks: consider raising
  `exitReadTimeout` and/or reducing concurrent chunk count at 10x-100x this
  pilot's scale if retries become frequent enough to matter.
