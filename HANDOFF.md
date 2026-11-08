# DeltaGain_Fungi — Handoff / Status (2026-11-08)

## Update 2026-11-08: enzymes_UPO brought under version control, reviewed, and a coordination warning

Session: https://claude.ai/code/session_01DDHWqZxJoeYpu2sDQS3BGF (laptop, `novelty-enzymes-50`).
Concurrent session doing the compute: https://claude.ai/code/session_01CoSBcRoaY8ycLQsGMS4hbp (on
`skylark`, not addressable by message from the laptop).

### Why this update exists

On 2026-11-07 the UPO project was rewritten on the HPCC — new framing, two measured negative
results, a new tool — and **none of it was in git**. About eleven hours of documents and code sat
untracked on shared storage. It is now committed.

### What was committed

| Commit | Contents |
|---|---|
| `c090229` | The 2026-11-07 HPCC work: `DESIGN.md`, rewritten `PLAN.md` (+ `PLAN.superseded.md`), `CALIBRATION_MATRIX.md`, `POCKET_TRANSFER.md`, `RESULT_T2_negative.md`, `SUBSTRATE_CALIBRATION.md`, `G1_BFD_SUMMARY.md`, `upo_sites/`, `upo_reg/`, small evidence data. 37 files |
| `4dcd560` | `REVIEW_PLAN_superseded.md` — an independent review of the original plan that had existed only on the laptop and had never been seen by the work that followed |
| `7d1e89a` | `SYNTHESIS_REVIEW.md` — all plans and reviews assessed against the full BFD run |
| this one | `COORDINATION.md` and this handoff section |

A scoped `analysis/enzymes_UPO/.gitignore` keeps out `data/bfd_full/` (~85 MB, regenerable),
`logs/`, the 9.3 MB publisher supplement and its figures, and the `.lineage.tsv` cache. The two
tables transcribed **by hand** from the supplement figures are committed, since they cannot be
recovered without redoing that work.

### Where the science stands

The project ran two activity gates and **both came back negative**:

- **T2** (aromatic ring hydroxylation, the original target): the "aliphatic pocket" does not
  predict it. Pearson r = **+0.176** across the 19-enzyme calibration panel.
- **T1** (peroxygenase-vs-peroxidase propensity, the thing it pivoted *to*): raw r = +0.748, but
  **partial r = +0.129 given clade**, and +0.001 within the ten non-Clade-V enzymes. Clade carries
  the signal; the pocket is nearly collinear with it and adds almost nothing.

The response was to narrow the claim rather than rescue the model. The deliverable is now a
**clade annotator with site evidence**, and the refusal to emit a per-protein activity score is
written into the tool's own output header with a measured reason. That is the right call and it is
the most valuable thing in the directory.

Full BFD run completed: **97,520 PF01328 proteins across 15,137 genomes** (8,314 genomes carry
none). Max copy number 144, in *Sphaerobolus stellatus*.

### The numbers to be careful with

- **Only 8.05 % of clade calls fall inside the caller's validated identity range.** The backbone is
  19 sequences; BFD proteins sit at 25–42 % identity to anything characterized. Everything else is
  labelled extrapolation, correctly.
- **The real working denominator is ~11,200** (catalytic-core AND validated-range), or **~6,700**
  excluding provisional Clade IV — not 97,520 and not 34,561.
- **Clade IV supplies 38.8 % of all validated-range calls** (4,968 of 12,818) at a mean identity of
  62.1 %, ~20 points above every other clade — and it is the one clade that *cannot* be validated,
  having a single characterized member. Not yet explained. Check before any candidate list is drawn
  from it.
- **Clade I (33,424 proteins, 34 % of all calls) carries the universal Glu180 at only 16.8 %**,
  while all five characterized Clade I enzymes have it. `G1_BFD_SUMMARY.md` §6 reads this as the
  low-identity Clade I bin absorbing proteins that are not UPOs. That caveat needs to travel with
  the lineage table in §4, where it currently does not.

### Open items, highest value first

1. **Zamocky et al. 2022 (PMC9138132) and RedoxiBase appear in no current document.** Their ~174
   curated sequences would raise the clade backbone from 19, which is the only cheap route out of
   the extrapolation problem, and would supply the missing outgroup and the prior-art baseline in
   the same action. Highest-value unexploited item in the project.
2. **`DESIGN.md` is stale.** Its §10 still says "pivot to T1" and calls the S8 digitization "the
   gate on whether T1 is viable". That gate was run and T1 failed. `PLAN.md` is correct;
   `DESIGN.md`, the longest document and the one framed as the intent document, is not.
3. **Clade IV check** (above), before dereplication picks representatives.
4. **Gene-model QC was dropped in the rescope** while exposure grew from ~800 genomes to 15,137 of
   self-generated models. BFD's `gene_transcripts` already carries `is_partial`,
   `has_start_codon`, `has_stop_codon`; joining them is one query and would let the
   catalytic-machinery numbers be reported QC-clean versus flagged. Until then "65 % lack the full
   catalytic set" is not safely a biological statement.
5. **Thresholds are uncalibrated.** The 30 % identity floor and `--min-margin` are explicitly
   labelled heuristics. Low-identity Clade I proteins lacking Glu180 are a ready-made decoy set for
   calibrating them.

### Coordination warning — two sessions are in this directory

The `skylark` session is running compute right now and cannot be messaged from the laptop;
`COORDINATION.md` in `analysis/enzymes_UPO/` is the channel, and the user is relaying.

As of 2026-11-08 10:10 PT it had jobs **31865316** (`upo_nr`, DIAMOND vs nr, running 34 min) and
**31865695** (`upo_prom`, promoter extraction) in flight, and had added `G1_FINDINGS.md`,
`derep.sh`, `diamond_nr.sh` and five `upo_reg/bin/` scripts since commit `c090229`.

**Flag for the user:** `PLAN.md` §2 lists **D6 as "the one blocking decision"** — G3 regulation
scope, all-BFD versus RNA-seq-followable clades — and §4 says the two options "are different
projects, not different settings". Job 31865695 is processing **7,055 target genes in 3,274
genomes**, i.e. the broad option. D6 appears to have been settled by starting rather than by
decision. If that is intended, `PLAN.md` should stop listing it as open.

Note also that `G1_FINDINGS.md` (skylark) and `SYNTHESIS_REVIEW.md` (laptop) were written
independently and overlap. They do not contradict each other, but they should probably be merged.

## Update 2026-11-07: enzyme-discovery side projects (read this first for `analysis/`)

This section covers only `analysis/`. It does not change the delta-gain funnel status below.

Two side projects started from `analysis/novelty_enzymes/BACKGROUND.md`, a report written by an agent.
Treat BACKGROUND.md claims as unverified; at least one is wrong (see "Corrections" below).
Committed in `1d32ec5` and `d133436`, pushed to `main`.

### Project 1 — FunBGCeX repeat (`analysis/novelty_enzymes/`)

Plan: `FUNBGCEX_PLAN.md`. Nothing has run.

Claim to test: BFD annotation of genomes NCBI left unannotated finds biosynthetic gene clusters that the
published 1,990-genome FunBGCeX run could not see, because that tool needs GenBank input with translations.

Measured facts that make this feasible:
- Funannotate already wrote `genome_annotation/<out>/predict_results/<out>.gbk.gz` for **23,450 genomes**,
  with `/translation` present (checked in one genome only). No format conversion is needed.
- **antiSMASH has not run for BFD.** The modules `BFD/ANTISMASH` and `BFD/MERGE_ANTISMASH` exist, but there
  is no output under `results/function/` and no table. The paper's headline contrast is FunBGCeX vs antiSMASH
  (1,050 vs 325 Pyr4 clusters), so this comparator must be decided.
- NCBI baseline genomes in `1KFG/2026/NCBI_fungi/source/NCBI_ASM/<acc>/` have `.fna`, `.gff` and `.faa` but
  **no `.gbff`**, so the NCBI-vs-BFD comparison needs a conversion step or a separate download.

Next step: Phase 0 of the plan — install FunBGCeX, fetch Zenodo 10.5281/zenodo.8126803, and run ~20 genomes
that appear both in the paper's Supplementary Data 4 and in BFD, once per annotation source. That measures
our gene-model error rate in cluster context before any large run.

### Project 2 — UPOs (`analysis/enzymes_UPO/`)

Plan: `PLAN.md` (user decisions are at the top), with `CENSUS.md`, `REFERENCES.md` and `sql/`.
Tool: `upo_id/` (built and tested; see `METHOD.md` and `README.md`).

Census (measured, `sql/census.sql`, `sql/novel_lineages.sql`):
- Join key is `samples.ASMID = prior_annotation_status.folder`. 24 of 23,680 genomes do not match; all 24 are
  Dikarya, so no early-diverging genome is lost.
- 827 early-diverging genomes (Microsporidia excluded); 813 have proteins; 10.17 M proteins.
- BFD Pfam at the gathering threshold finds **807 PF01328 proteins** there; 309 are in BFD-only genomes.
- Of the 23,451 genomes with BFD gene models: 17,864 had no prior annotation, 5,553 were NCBI-annotated.
- Newly sampled Dikarya: 6 classes, 22 orders and 107 families have **no** NCBI-annotated genome.
  Lecanoromycetes (lichens) is the largest new sample at **340 BFD-only genomes** vs 137 NCBI-annotated.

**The first real result, and it is negative for the original framing.** Running `upo_id` on the 807 proteins
gives 665 CPO-like Family I, 133 that fail the motif checks, 6 low-confidence UPO Family II and 3 UPO
Family I. **All 9 UPO calls sit in NCBI-annotated genomes. None is in a BFD-only genome.** So for the three
phyla in the original request, BFD adds no new UPO. The 6 Family II calls are Mucoromycota
(*Jimgerdemannia*, *Endogone*), about 300 aa, and lack the full long-type signature.

Next step: run `upo_id` on the newly sampled Dikarya lineages, starting with the 340 BFD-only
Lecanoromycetes genomes. That is a short job and it is the cheapest test of whether the project still has a
novelty claim. Not checked: whether lichen assemblies carry photobiont or bacterial contigs, which would
affect gene calls.

### Corrections to BACKGROUND.md

- "113 putative UPOs from more than 800 genomes" is a mis-citation. 113 is the number of seed sequences in
  Faiza et al. 2019 (BMC Evol Biol 19:76).
- The EHD/EGD motif conflict is settled. Five crystal structures show EHD is the Mg site in short UPOs
  (MroUPO, HspUPO, CviUPO), so EHD does not exclude a UPO.
- PF01328 alone is not a UPO call. It covers 97,520 proteins across BFD.
- UPObase (`upobase.bioinformaticsreview.com`) does not resolve. The user decided to build our own reference
  set rather than ask the authors.

### Open work

1. **Tree placement** replaces HMM score competition in `upo_id` (user decision 7). Not implemented.
2. **Residue-probing literature** (Allozymes case study; PMC12117455 and its references) to define which
   positions to score. Not started.
3. **CglUPO and CciUPO accessions** are unverified (noted for a later lookup).
4. **No plan or tool has had an independent review.** Three attempts were stopped by model safeguards
   (Fable 5.1 twice, Opus 5.5 once), all with `[bio]`. Request IDs: `req_011CgpSgRiZJzvitHNnm88ev`,
   `req_011CgpTutXr2zigZNSVDRXhs`, `req_011CgpUDDkhzHdk3rJ1ypPUL`. A later attempt to start the
   residue-literature agent was also stopped. Run reviews from a fresh session, or report the false positives.

### Conventions used here

- `upo_id` needs `module load hmmer/3.4` (plus `mafft/7.505` for `build`) and `/usr/bin/python3.12`.
  The default `python3` on HPCC is 3.9 and will not work.
- Files that no tool reads at runtime are zstd -19. HMM files stay plain because `hmmsearch` and `hmmalign`
  cannot read compressed input. `.gitignore` carries an exception so the generic `results/` rule does not
  drop `upo_id/results/`.

## Update 2026-11-05: strategy change — baseline funnel replaces the universal search (read this first)

Full plan, definitions, workflow diagram and measures of success:
`DeltaGain/docs/DELTA_GAIN_PLAN.md`. This section replaces the 11-04 section where they disagree.

- **Full DIAMOND run stopped** (about 15 % of 2,614 searches done). It searched every representative against all of
  NCBI nr clustered (`nr_cluster_seq.dmnd`, `--taxonlist 4751`, so the full 195 GB index is still scanned), which
  measures "not yet public", not delta gain. Both `work/` directories were deleted; raw hits were kept zstd-compressed
  in `results/protein_novelty_n2000/raw_hits/` (177 nr + 177 UniProt) and
  `results/protein_novelty_full/raw_hits_partial/nr_cluster/` (196 finished chunks). Full linclust output is kept in
  `results/protein_novelty_full/`.
- **New design:** baseline = NCBI proteins of genomes with prior annotation; query = BFD proteins of genomes with none.
  Nested MMseqs2 linclust tiers (95/70/50 %) on baseline + query together; novelty counted per protein. Scripts:
  `DeltaGain/bin/{classify_genomes_by_prior_annotation,build_funnel_input,delta_funnel_summary}.py`,
  `DeltaGain/launch/delta_funnel.sbatch`. Smoke test passed (40 + 40 genomes, 2m13s). Full run not launched yet.
- **Prior-annotation table** `prior_annotation_status.tsv` (27,401 NCBI assemblies): 5,511 baseline_ok, 158 annotated
  without protein file, 1,228 annotated never downloaded (1,184 are same-strain duplicates dropped by NCBI_fungi's
  dedup), 18,145 no prior annotation (the new set), 2,359 not downloaded. A gap download (`scripts/fetch_gap.sh`, 202
  accessions) was running; 8 accessions have annotation info but no protein file. Refresh the table after it finishes.
- **Re-bin of the 2000-genome hits** (`scripts/dg_rebin.py`): the old rule is reproduced (37 of 10.6M differ);
  best-class moves 39 % of "weak" to strong (strong 78.0 -> 84.9 %, no-hit unchanged at 4.35 %). Genomes with prior NCBI
  annotation have 47.5 % of queries with a 99-100 % top hit vs 12.6 % for unannotated genomes: self-hit inflation.
- **Hardware:** never run compute on `abu_dhabi` nodes (all of `batch`; c02/c21 in `short`/`preempt`/`stajichlab`). No
  AVX2, so MMseqs2/DuckDB die with "Illegal instruction". Use `-p epyc,intel` with
  `--constraint="ryzen|broadwell|cascade|sapphire"`.
- **Next:** finish the gap download and refresh the status table -> launch the funnel -> fill the measured counts in the
  plan -> write the survivor-search stage -> Pfam family walk. Open questions are listed in section 13 of the plan.

## Update 2026-11-04: full-set run, pangenome scripts, run-state (read this first)

This section replaces the 10-24 section below where they disagree.

**Deadline:** the cluster shutdown reservation starts **2026-11-10 08:00** (all
nodes, `20261110-Shutdown`). Jobs whose time limit overlaps it stay pending
("ReqNodeNotAvail, Reserved for maintenance"). Head jobs are submitted with
`-t` of 5 days or less so they end before it. Nextflow `-resume` keeps finished
work, so anything unfinished resumes after the shutdown.

**State of the input data:** BFD annotation is finished per the user. The
`Fungi_BFD_runs/input/pep/` count was still 23,448 files with none modified
after 2026-10-24 14:08 when checked on 11-04. The user is building the final
merge/parquet tables in `Fungi_BFD_runs/tables/` (still the 2026-09-10
snapshot when last looked at). The accumulation/pangenome analysis needs the
new `asm_stats`, `busco_genome`, the Stage 1 classification and
`bfd_version_tag`; those params are still unset.

**Runs (as of 2026-11-04 ~17:00; check `squeue`/`sacct`, these were just submitted)**
- **2000-genome run, relaunch:** job 31677214, `-resume`, outdir
  `results/protein_novelty_n2000`, launched from `DeltaGain_Fungi/`. The earlier
  head (31662005) ended FAILED after 7h09m because `COMPARE_CLUSTERINGS` was
  OOM-killed (exit 137, 8 GB request; the standalone run peaked at 12.5 GB).
  Memory is now 32 GB in `conf/profile_protein_novelty.config`. With
  `errorStrategy 'finish'` that one failure stopped new submissions, so
  `bfd_novelty_bins.tsv` did not exist yet. Not verified: how many
  DIAMOND chunks had finished (268 exit files were seen in work/, which includes
  both search types and older runs).
- **Full-set run (all annotated genomes, `--n_test 0`):** launched from
  `runs/full/` (its own launchDir, because two Nextflow heads cannot share one
  `.nextflow` cache) with `-c full.config`, outdir `results/protein_novelty_full`,
  DIAMOND cluster branch off. Job 31676864 failed twice over:
  1. `COMBINE_BFD_PROTEINS` was sent to the `short` partition (2 h cap) with an
     8 h limit, so sbatch refused it. Fixed in `runs/full/full.config`
     (queue `batch`); same for `SPLIT_FASTA` and `BIN_NOVELTY_HITS`.
  2. `FETCH_UNIPROT_FUNGI (reviewed)`, the REST smoke test, got HTTP 429 from
     UniProt and failed. New switch `--run_uniprot_smoke_test false`.
  Job 31676864's head was left alive to finish `COMBINE_BFD_PROTEINS`
  (job 31676885). Relaunch is job 31677215 with `--dependency=afterany:31676864`
  and `-t 4-12:00:00`, `-resume`.
- **Correction (2026-11-04 evening):** `--run_uniprot_smoke_test false` and
  `--run_diamond_cluster false` did NOT work in the relaunches above. Nextflow
  26.04.6 passes `--flag false` as the string 'false' (truthy), so the REST
  smoke test failed again (429) and `DIAMOND_CLUSTER` started on the full set.
  I cancelled that DIAMOND_CLUSTER job (31690829). Fixed in DeltaGain commit
  "Parse boolean switches explicitly" (`flagOn()`). Not yet confirmed that a
  relaunch honours it.
- **2000-genome run, third launch:** the second head (31677214) got through all
  DIAMOND chunks, then `BIN_NOVELTY_HITS` was OOM-killed (137, 4 GB, 10.6M
  reps). Memory raised to 32 GB (a guess above the failure point, not a
  measurement). Relaunched as job 31703284 (`-resume`).
- **Full run, relaunch:** job 31703285, dependent on head 31677215 ending.
  `MMSEQS_LINCLUST_BFD` (31690831) was still running under 31677215. Full-run
  `BIN_NOVELTY_HITS` raised to 160 GB on `epyc`. `bin_novelty_hits.py` holds
  every query's best hit in per-target dicts plus a list of all query ids;
  measured 495 MB per 1M queries (2 targets, shared id strings; real figure
  probably ~700 MB), so ~60-90 GB at ~126M representatives (the representative
  count is itself extrapolated). If this ever needs to be smaller, each query
  appears in exactly one chunk per target, so the script could stream chunk by
  chunk. Also note it keeps the lowest-E-value hit and classifies only that,
  so a high-coverage hit with a slightly worse E-value can be masked.
- **2000-genome result (done 2026-11-04 21:00, job 31703284):**
  `results/protein_novelty_n2000/bfd_novelty_bins.tsv`, 10,615,362 queries (equals
  the linclust representative count). Overall: strong 80.0%, weak 15.7%,
  no_hit 4.31% (457,750). Pilot (100 genomes) was 71.0 / 22.6 / 6.5%. Per
  genome, no-hit representatives fell from ~430 (pilot) to ~229. The samples
  differ (pilot = first rows, this = seeded random), thresholds are still
  placeholders, and no-hit proteins are candidates only (no repeat or
  contamination filter yet). UniProt-fungi alone rescues only 3,896 proteins
  (0.04%) that `nr_cluster` fungi missed. The masking issue (lowest-E-value row
  only) can move strong<->weak but not no_hit.
- **Full run measured:** `MMSEQS_LINCLUST_BFD` took 1h25m with MaxRSS 224 GB
  (the 400 GB request was right; the 64-128 GB guess in the 10-24 note was not).
  1,307 chunks of 60,000, so ~78M representatives, which makes 2,614 DIAMOND
  searches. First 7 finished: UniProt median 21 min, nr_cluster 35 min (n=1).
  At ~30 min and ~23 concurrent, roughly 2.4 days; Nextflow `queueSize` keeps
  ~500 submitted at a time.
- **Streaming binning (new, off by default):** `bin/bin_novelty_hits_stream.py`
  + `BIN_NOVELTY_HITS_STREAM`, selected with `--bin_streaming true`. Memory is
  bounded by chunk size. With `--classify best-evalue` (default) it reproduces
  the original script exactly (tested); `--classify best-class` classifies every
  row and keeps the best class, which changes results and is not the default.
  The running full-run head was launched before this existed, so it will use
  the original `BIN_NOVELTY_HITS` (160 GB on epyc, enough for ~78M queries by the
  0.5-0.7 GB per 1M estimate). To use streaming, resume with `--bin_streaming true`.
- **Sizing in `runs/full/full.config` is extrapolated, not measured:**
  linclust 400 GB (25.6 GB MaxRSS at 20M proteins, scaled about 11.7x linearly;
  the 10-24 note of 64-128 GB does not follow from that measurement),
  `COMBINE_BFD_PROTEINS` 8 h (8m23s at 20M), `SPLIT_FASTA` 6 h / 8 GB and
  `BIN_NOVELTY_HITS` 12 h / 64 GB (both guesses).
- **Capacity:** the account CPU cap (`AssocGrpCpuLimit`) was already full of the
  2000-genome chunk jobs on 11-04. The two runs compete for it. The user chose
  to let both run (no holds).

**Code changes (DeltaGain, committed locally on HPCC, not pushed)**
- `26cd0e1`: `bin/species_pangenome.py`, `bin/compare_species_pangenomes.py`,
  their tests, and `params.run_diamond_cluster` (default true).
- `55f9b9d`: `params.run_uniprot_smoke_test` (default true).
- DeltaGain_Fungi `75b6e41`: `compare_clusterings` memory, `runs/full/full.config`,
  `.venv-pangenome/` in `.gitignore`.

**Species pangenome analysis (tested on synthetic data only; reviewed by Fable 2026-11-04)**
Question from the user: build a consolidated pangenome per species, then
compare pangenomes across species. Both scripts were rewritten after a review
(commit message in DeltaGain has the full list).
- `bin/species_pangenome.py`: per species (grouped by `samples.csv` SPECIES) with
  at least `--min-genomes` (default 10): pangenome size; core / soft-core
  (>= ceil(0.95n), blank when n < 20) / shell / cloud (<= max(1, floor(0.15n)),
  singletons always cloud); `rarefied_pan_size` (+sd) and `rarefied_core` at
  `--rarefy-n` (default 10); alpha of the new-cluster rate fitted on the
  permutation-mean curve for N = 2..n, with `alpha_n_points` and `alpha_r2`.
  alpha < 1 = open. A species that adds nothing after genome 1 is "closed"
  (`no_new_after_first`), fewer than 3 usable points is "unresolved". Placeholder
  names ("sp.", "cf.", "aff.", "uncultured", "unclassified") are skipped. Alpha is
  not comparable across very different n; use the rarefied columns.
- `bin/compare_species_pangenomes.py`: exact pairwise shared clusters and Jaccard
  (sparse products), plus `*_narrow` versions excluding clusters in more than
  `--max-species-per-cluster` (100) species; same-genus / same-family flags;
  per-taxon open/closed/unresolved counts and medians. It fails if its genome
  counts disagree with the species table (pass the same `--exclude-asmids` to
  both steps). Lineage conflicts go to `lineage_conflicts.tsv`.
- Scale test (synthetic): 3,000 genomes, 13.2M pairs, 50 species: 12 s / 1 GB for
  the species step with 4 processes, 4 s / 1 GB for the comparison. Extrapolated
  (not measured) to ~230M pairs: minutes, ~15-20 GB.
- Input is `incidence_matrix.parquet` from `BUILD_CLUSTER_INCIDENCE_MATRIX`, which
  only runs once the four accumulation params above are set. Neither script has
  run on real data.
- **Still open:** clustering identity (95/90 linclust; a looser per-species cutoff
  needs another clustering and was not launched); whether to drop Stage 1
  `redundant_strain` genomes first (needs a full-scale classification; only the
  100-genome pilot exists); the `--min-genomes` / `--rarefy-n` defaults are
  placeholders; taxon labels are keyed by name only (homonyms merge); an
  alpha uncertainty estimate (bootstrap over genomes) is not implemented; a
  random 2000-genome sample leaves most species with too few genomes, so use the
  full run. Cross-species sharing at 95% identity mostly reflects HGT,
  contamination or recent splits, so deep pairs will be near 0.
- Python env: `uv` env at `.venv-pangenome/` (Python 3.12, numpy, pyarrow, pandas,
  duckdb, scipy, pytest). Run with
  `PYTHONPATH=DeltaGain/bin .venv-pangenome/bin/python ...`.

**Next steps**
1. Check both relaunched heads start and the full run's `COMBINE_BFD_PROTEINS`
   is reused from cache. Read real MaxRSS/elapsed from `sacct` for linclust
   and the combine step before trusting the extrapolated sizing.
2. When the 2000-genome run ends: compute strong/weak/no-hit fractions and
   compare with the pilot (71.0 / 22.6 / 6.5%).
3. Before 11-10 08:00: note which full-run tasks are done. After the shutdown,
   resume the full run from `runs/full/` with the same command.
4. After the user's tables are built: set the four accumulation params, run
   `species_pangenome.py` then `compare_species_pangenomes.py` on the full run.

## Earlier update 2026-10-24: 2000-genome protein screen

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

### DIAMOND cluster vs MMseqs2 linclust (2000 genomes, measured 2026-10-24)

The orphan `DIAMOND_CLUSTER` job 30452869 finished (exit 0, 47 min 58 s, MaxRSS 21 GB).
Job 30476464 compressed its table and ran `compare_clusterings.py`.

| metric | mmseqs_linclust | diamond_cluster |
|---|---|---|
| members | 20,373,920 | 20,373,920 |
| clusters | 10,615,362 | 12,947,923 |
| singletons | 8,473,455 | 11,478,069 |
| max size | 1477 | 1464 |

- Adjusted Rand Index: 0.7291 (pilot: 0.7155). Both tables hold the same members.
- DIAMOND cluster gives more, smaller clusters. I did not test why. The two tools use different identity and coverage definitions.
- Files: `results/protein_novelty_n2000/bfd_proteins_diamond_clusters.tsv.zst` (87 MB, passes `zstd -t`, 20,373,920 lines) and `bfd_clustering_comparison.tsv`.
- The compare job used 12.5 GB peak memory (64 GB requested) and took 4 min 26 s.

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
