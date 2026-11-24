# pathway_census: pathway gene census for BFD genomes

Shared package that counts the genes of one biosynthetic pathway in every BFD genome. Carotenoid is
the pilot pathway. No real-data run exists yet; all tests use synthetic fixtures. The real-data
checklist is in `docs/superpowers/plans/2026-11-20-pilot-lessons.md`.

Pipeline (one CLI per stage). Every CLI that reads the pathway YAML, except
`list_unverified_refs.py`, refuses a `verified: false` pathway unless `--provisional`, and then tags
its outputs `.PROVISIONAL`. `merge_absence.py` reads no pathway; it carries the tag of its inputs.

1. `pathways/<id>.yaml` + `refs/<id>/`: definition and references (`pathway_census.defs`, `bin/list_unverified_refs.py`).
2. `bin/build_step_hmms.py`, `bin/search_step.py`: step HMMs and candidate search (`hits`).
3. `bin/assign_steps.py`: step labels from candidates plus a placements TSV (`assign`).
4. `pathway_census.quality.quality_flags`: quality, contamination and haplotype flags (library only).
5. `bin/build_matrix.py`: genome x step copy matrix and phylum summary (`matrix`).
6. `bin/absence_check.py` + `bin/merge_absence.py`: tblastn check of absences (`absence`).
7. `bin/paired_compare.py`: deposited versus BFD annotation on the same assemblies (`paired`).
8. `bin/triage.py`: triage bins (`bins`). `chapters/A_carotenoid.qmd`: report scaffold.
The genome table (`pathway_census.genomes`) feeds stages 5-7.

Run tests with the system Python 3.12 (default `python3` is 3.9):

    /usr/bin/python3.12 -m pytest analysis/pathway_census/tests -v

`pytest` and `pyyaml` were already present under `/usr/bin/python3.12`, so no venv was needed.
If they go missing, create `analysis/pathway_census/.venv` (git-ignored) and install the `pyproject.toml` dependencies (`pyyaml pandas numpy duckdb`, plus `pytest` from the `test` extra) there.

## Genome table (Task 2)

`pathway_census.genomes.load_genome_table(tables_dir, prior_status_tsv)` returns one row per `species.parquet` row (left joins).
Query: `sql/genomes.sql`. Schemas and join key are recorded at the top of `pathway_census/genomes.py`.
Join key: BFD tables join on `ASMID`; the prior-annotation TSV joins on the accession prefix of `ASMID` (`GCA_/GCF_` + number + version) = TSV `accession`.
Columns: `locustag, asmid, folder, species, phylum, class_, genus, genus_source, prior_status, prior_folder_state, annotation_provider, busco_complete, n50, downloaded, has_proteins, has_deposited_proteins`. `asmid` comes from `species.parquet`. `folder` comes from the prior-annotation TSV (`<NA>` without a TSV row); it names the deposited NCBI_ASM folder.

Regenerate (writes `results/genomes.tsv.zst`, git-ignored):

    cd analysis/pathway_census && /usr/bin/python3.12 -m pathway_census.genomes

MEASURED 2026-11-20 on the real tables:
- Rows: 24,363 (= `species.parquet`; 24,363 unique LOCUSTAG). `busco_genome` has 23,582 rows and `asm_stats` 23,584; 814 genomes have null `busco_complete` and `n50`.
- Rows per `prior_status` (4 of the 5 TSV statuses occur among BFD genomes): no_prior_annotation_present 18,055; baseline_ok 5,625; no_prior_annotation_not_downloaded 673; annotated_folder_missing_protein 10; annotated_NOT_downloaded 0; unknown 0.
  The TSV has 27,401 rows; 3,038 of them are not in `species.parquet`.
- Empty GENUS: 349 rows. genus_source counts: genus_field 24,014; parsed_species 49 (includes 8 Antarctolithica); unresolved 300. Unresolved rows get a non-null label "<first token> sp." (e.g. "Ascomycota sp."). Genus is never null.
- `has_proteins` = LOCUSTAG in `seqdb/genomes.tsv` and the BFD `proteins.fa.gz` exists on disk: True 23,379; False 984. Checked: 23,379 of the 24,363 species LOCUSTAGs are in seqdb/genomes.tsv (23,412 rows) and in `gene_info`; all 23,379 protein files exist.
- `downloaded` = BFD genome FASTA (`<out>.scaffolds.fa.gz` beside the proteins file) exists on disk: True 23,379; False 984 (the same genomes as has_proteins=False). Genome FASTA availability is derived from this file path, not from the prior-annotation TSV. Genomes outside seqdb have no BFD FASTA path in the tables; a copy may exist elsewhere (e.g. NCBI_ASM). That was not checked.
- `has_deposited_proteins` (TSV `protein_file_ok == 1`): True 5,634; False 18,729.
- `prior_folder_state` (TSV): present 23,690; not_downloaded 673.
- Crosstab prior_status vs has_proteins=True: baseline_ok 5,535 of 5,625; no_prior_annotation_present 17,834 of 18,055; annotated_folder_missing_protein 10 of 10; no_prior_annotation_not_downloaded 0 of 673.

## Carotenoid definition (Task 3, DRAFT for user review)

- Schema v2 (`pathway_census.defs`): `domains` is renamed `domains_any` (any-of; the loader rejects `domains` with "use domains_any"). New optional step keys: `route`, `reactions`, `decoys` (same fields as a ref), `lineage_note`. New optional ref keys: `gene_name`, `region` (`[start, end]`, 1-based inclusive; checked by `validate` against `refs/<id>/<step>.faa`; skipped only when the refs directory does not exist; with a refs directory, a region without an accession, a missing FASTA or an accession absent from it is an error; `load_pathway` assumes the YAML is in `<root>/pathways/` and takes `refs_dir` otherwise). New top-level `routes` list (`id`, `name`, `step_ids`, `note`). `validate` reports unknown route steps, a step whose `route` does not list it in `step_ids` (`Route.step_ids` is the membership; a step may be in several routes and may have no `route`), duplicate decoy accessions, undeclared step routes, duplicate reaction ids in a step, and a decoy accession that is also a ref of the same step.
- `pathways/carotenoid.yaml` (step set v2, 2026-11-21): 14 steps in 7 routes (`core`, `neurosporaxanthin`, `retinal`, `astaxanthin`, `trisporoid`, `trisporoid_tsp3`, `trisporoid_late`); each route is a step order that `kb.path_possible` can follow, except `trisporoid_late` (its compounds are not in ChEBI), pathway `verified: false`. User review of the step list is required before any census run (spec 5.8).
- `refs/carotenoid/<step>.faa`: reference proteins from UniProt (retrieved 2026-11-21), one file per step. The five fused phytoene synthase / lycopene cyclase proteins are cut per domain at the UniProt "Region" features: `phytoene_synthase.faa` holds the synthase regions and `lycopene_cyclase.faa` the cyclase regions. Full sequences: `refs/carotenoid/full_length/fused.faa`.
- `region` in the YAML is in coordinates of the FASTA record (`defs.validate` checks it against `<step>.faa`). For the cut records it covers the whole record; the coordinates on the full protein are in the `region` column of `provenance.tsv`.
- `refs/carotenoid/provenance.tsv`: one row per (step, ref) with gene name, region on the full protein, source and verified flag.
- `refs/carotenoid/decoys/<step>.faa` and `decoys/provenance.tsv`: proteins of the same family with another function, with the reason for each. They are outgroups for tree placement, not references.
- `refs/carotenoid/outgroups/<step>.faa` and `outgroups/provenance.tsv`: true homologs that do the reaction of the step, kept only as outgroups for tree placement (column `role` = `tree_outgroup`). They are not decoys. The absence check and the decoy-calibrated `min_score` do not read this directory, and the YAML has no key for it. Today it holds one entry: `lycopene_cyclase` Q9HNE5.
- `kb/carotenoid/*.tsv`: compounds, reactions, participants, steps and the (protein, reaction) links with evidence tier, PMID and one quoted sentence (see `kb/README.md`).
- `curation_scripts/carotenoid/`: one-off build helpers that fetched the records and wrote these files. Not run by the tests (see its README).
- `refs/carotenoid/UNVERIFIED.md`: output of `bin/list_unverified_refs.py pathways/carotenoid.yaml` plus a hand-written list of gaps, held tiers, excluded (T4) steps and open decisions.

## Candidate retrieval (Task 4)

Code: `pathway_census/hits.py`, `bin/build_step_hmms.py`, `bin/search_step.py`, `slurm/search_step.sbatch`.

Tool modules found on UCR HPCC (`module avail`, 2026-11-20).
- hmmer: 3.3.2, 3.4 (default), plus -mpi builds and old 1.8.5 / 2.3.2. Used: `hmmer/3.4` (HMMER 3.4, Aug 2023).
- mafft: 7.490, 7.505. Used: `mafft/7.505`.
- ncbi-blast: 2.2.22+ through 2.17.0+ (not used in Task 4; recorded as requested).
- mmseqs2: `mmseqs2/17-b804f` (also 13-45111, 15-6f452). Used for dereplication.
- cd-hit: 4.8.1 (available, not used).
`hits.run_tool` runs a tool from PATH, or under `bash -lc 'module load hmmer/3.4 mafft/7.505 mmseqs2/17-b804f && ...'`
(override with `PATHWAY_CENSUS_MODULES`). Tests that need these tools skip with a reason if they cannot be loaded. They never mock a search.

Domtbl parsing: `score` is the per-domain score (domtbl column 14), `evalue` is the domain i-Evalue (column 13).
`candidates_for_step` returns one row per protein (its best-scoring domain) with score >= `min_score`. `min_score` has no default.
A zero-hit search returns an empty DataFrame with all eight columns.
`bin/search_step.py` adds a `step_id` column (the searched step) to the candidates TSV, so the file feeds `assign_steps` directly (final-review fix I4, part). `candidates_for_step` itself is unchanged. Output names still carry the step only, not the genome: see the lessons-note checklist (step 3) before a real run.

### Provisional mode
`build_step_hmms.py` and `search_step.py` refuse a pathway with `verified: false` (exit non-zero) unless `--provisional` is passed.
With `--provisional`, output names contain `.PROVISIONAL` (e.g. `ggpp_synthase.PROVISIONAL.hmm`, `.domtbl`, `.candidates.tsv`),
the HMM directory gets `PROVISIONAL.txt`, and the candidates TSV begins with a `# PROVISIONAL` line.
The carotenoid pathway is `verified: false` until the user reviews it (spec 5.8), so no real HMMs or searches were run in Task 4.

### min_score rule (PROPOSED, not decided)
For each step, the pilot will score the step's own reference proteins against the step HMM and set `min_score`
to the lowest reference score (so every reference is retrievable). The value is then recorded per step. Task 7 sets it.
Open: refs are the HMM's training set, so their scores are optimistic. A leave-one-out check may be needed. No data yet.

### Dereplication (PROPOSED values)
`hits.dereplicate(fasta_in, fasta_out, identity, coverage=0.8)` runs `mmseqs easy-cluster` (`--min-seq-id identity -c coverage --cov-mode 0`) and writes one representative per cluster.
Rule (PROPOSED): a family whose protein set has more than 50,000 proteins is dereplicated at 90% identity before the HMM search.
Smaller families are searched against the full BFD protein set without dereplication or a Pfam prefilter (spec 5.4 step 2).
Both numbers are proposals with no measured data behind them. Revisit after the first real family size is counted.
The rule is not enforced in code; the caller applies it.

### SLURM
`slurm/search_step.sbatch` takes all paths from environment variables (`PC_ROOT`, `PC_PATHWAY`, `PC_STEP`, `PC_HMM_DIR`, `PC_PROTEINS`, `PC_OUT_DIR`, `PC_MIN_SCORE`, optional `PC_PROVISIONAL=1`). It stages proteins in `$SCRATCH` and copies results to `PC_OUT_DIR`. Not yet submitted.

## Step assignment (Task 5)

Code: `pathway_census/assign.py` (`assign_steps`), `bin/assign_steps.py`. Tests: `tests/test_assign.py`.

### upo_place decision (read 2026-11-20)
`analysis/enzymes_UPO/upo_place` does NOT fit as is. Evidence (from its README and `bin/upo_place.py`):
- Input: PF01328 proteins, one fixed 187-leaf reference tree (`data/refset_zamocky2022/tree/`), one HMM from a 224-column alignment, and `clade_members.tsv` (Clades I-V, Zamocky Figure 1).
- Four subcommands (argparse in `bin/upo_place.py`): `match-cols`, `prune`, `assign`, `loo`.
- `assign` takes `--clades` and `--min-lwr` and reads any jplace file. It is reusable per step with a new `clade_members.tsv`. Output: clade call (LWR sum inside a named clade >= `--min-lwr`, default 0.5, else `unassigned`) and a `long_pendant` flag.
- UPO-specific parts: `match-cols` (224-column trimming) and the fixed reference paths in the slurm scripts. There is no step table and no per-step reference selection, so the pipeline as a whole does not run unchanged on carotenoid steps.
The reusable method: `hmmalign` to match columns, `epa-ng` placement, `assign`. Each carotenoid step needs its own reference alignment, tree and `clade_members.tsv`.

### Replacement tool (follow-up, not implemented)
`module avail epa-ng pplacer iqtree fasttree raxml` on 2026-11-20 found: iqtree 1.6.12, 2.1.1, 2.1.3, 2.2.0-2.2.2.6, 3.0.1; fasttree 2.1.11; raxml 8.2.12; raxml-ng 1.1.0, 1.2.0, 2.0.2 (also RAxML_NG 1.1.0-2.0.2). No `epa-ng` or `pplacer` module. `upo_place` gets epa-ng through `pixi exec` (pinned in its `slurm/env.sh`); that route was not tested here.
Plan (PROPOSED): per step, build tree with iqtree/raxml-ng on the step references plus family homologs, place with epa-ng via pixi, then apply the `upo_place.py assign` rule. No wrapper was written because no module provides a placement tool. `assign_steps` reads a documented placements TSV instead.

### Interfaces
- `placements` TSV columns: `protein_id, step_id, clade, support` (support = LWR sum in the clade). Empty clade = outside every clade.
- `candidates`: Task 4 candidate tables concatenated, with an added `step_id` column (the step whose HMM gave the row). Only this table supplies domain hits.
- Build `candidates` by concatenating per-step outputs: `pd.concat([candidates_for_step(d, min_score).assign(step_id=sid) ...])`. A missing `step_id` column raises ValueError.
- Output columns `protein_id, step_id, support, label, fusion_candidate`. One row per (protein, step).
- CLI `bin/assign_steps.py` refuses `verified: false` without `--provisional`; with it, `--out X.tsv` becomes `X.PROVISIONAL.tsv` and the file starts with a `# PROVISIONAL` line.
- Label rule: `assigned` only if the (protein, step) pair is in `candidates` (HMM evidence) AND a placement row for it has a clade and support >= `min_support` (inclusive; NaN fails). Otherwise `unassigned`.
- Fusion: HMM hits to step A and step B with B in `A.fused_with` (or reverse) and alignments overlapping < 50% of the shorter (PROPOSED threshold, no data) set `fusion_candidate = True` on both rows. Fusion only decides which steps a protein may carry. A fused protein is `assigned` for a step only with its own qualifying placement; with none it stays `unassigned` for those steps and keeps its rows. A protein with qualifying placements for both steps counts once toward each.
- A placement for a pair not in `candidates` is kept as `unassigned` (documented choice; no HMM evidence).
- A step not in the pathway (placements or candidates) raises `ValueError`. Support below `min_support` gives `unassigned` (kept). `min_support` has no default.
- Every candidate (protein, step) pair without a qualifying placement gets an `unassigned` row (support NaN if no placement). Rows are never dropped.

### Family size
Reference counts per step (MEASURED, `refs/carotenoid/*.faa`): largest is 5 (phytoene_synthase, lycopene_cyclase). BFD family sizes are UNVERIFIED: no search has been run (pathway `verified: false`).

### Deferred
Running on real candidates (brief step 6) is deferred until the user reviews `carotenoid.yaml`. No counts of assigned, unassigned or fused proteins exist.

## Quality, contamination and haplotype flags (Task 6)

`pathway_census.quality.quality_flags(...)` returns `protein_id, step_id, domain_ok, motif_ok, truncated, contam_suspect, hap_dup_suspect`. Rules and PROPOSED thresholds are in the module docstring. Flags are nullable booleans. `<NA>` means the needed input was not supplied. The unit tests use small synthetic fixtures. No real-data run has been done.

What is NOT implemented, and what a later run needs:
- No per-gene taxonomy. BFD `tables/` (listed 2026-11-20) has no per-gene taxonomy or lineage table. `swissprot_subject.parquet` has `ncbi_taxid` per SwissProt entry, and `swissprot_best.parquet` has a best hit per protein. A fungal vs non-fungal call from them needs a taxid-to-lineage map that is not in BFD. Also, the memory notes call the SwissProt table "old format". UNVERIFIED whether it is usable. It was not used.
- Contamination uses a caller-supplied `best_hits` table (`protein_id, taxon_class` = fungal / non_fungal). A later run must make it: search proteins of the scaffold neighbourhoods (not only the pathway genes) against a reference set labelled fungal vs bacterial/algal. That reference set does not exist yet. The search is not run by this module.
- The haplotype flag needs a caller-supplied `pair_identity` table (`protein_a, protein_b, pident, qcov`, fractions) from an all-vs-all search inside each genome, covering neighbour genes. It is not produced here. It judges gene neighbourhood duplication only. It does not align scaffolds.
- `motif_ok` is `<NA>` unless `motif_hits` is supplied. The carotenoid definition has no catalytic motifs yet.
- `gene_info` must be joined with protein `length` from `gene_proteins.parquet`. Without it `truncated` is `<NA>`.
- `domain_ok` needs `expected_domains` (step_id -> hmm names, from `domains_any` in the YAML). The check is any-of: True if at least one listed domain is present, False if none, NA if the list is empty or the step has no entry. Names are normalised before the match (uppercase, Pfam version suffix after a dot dropped, so `PF00494.22` equals `PF00494`); a non-list entry raises ValueError. It matches `hmm_name` in the `domains` frame, so the domain search must be run with the matching HMM library (for example Pfam ids from `pathways/carotenoid.yaml`).

Task 7 must filter with `quality.clean_mask(flags)` and never `fillna(False)`: `<NA>` means unknown, not clean. Task 7 must also report `quality.unknown_counts(flags)` in its output. Contamination is `<NA>` when fewer than `min_neighbours_with_hits` neighbours have a best hit.

## Copy matrix (Task 7)

Code: `pathway_census/matrix.py`, `bin/build_matrix.py`. Columns and rules are in the `matrix.py` docstring.
- Full genome x step cross product, zeros included. `has_proteins=False` genomes get zero rows and are excluded from `n_genomes`.
- `copy_raw` = assigned rows. Mutually exclusive buckets, priority flagged > unknown > filtered: `copy_flagged` = any required flag known True; `copy_unknown` = no flag True and some required flag NA; `copy_filtered` = all required flags known False (default require: truncated, contam_suspect). `copy_filtered + copy_unknown + copy_flagged = copy_raw` per row. `hap_dup_suspect` is counted in `n_hap_suspect` / `n_hap_unknown` only.
- `phylum_summary` uses `copy_filtered` by default (`count_col="copy_raw"` for the unfiltered table). `n_unknown` = genomes with unevaluated copies.
- `phylum_summary` columns: `phylum, step_id, n_genomes, n_with_step, median_copy, max_copy, frac_with_step, n_unknown, copy_raw_total, copy_filtered_total, copy_unknown_total, copy_flagged_total`. The four `*_total` columns are bucket sums over the genomes with `has_proteins` True; they do not depend on `count_col`.
- Without a contamination table, `contam_suspect` is NA for every copy. It is in the default `require`, so every copy is `copy_unknown` and `copy_filtered` is 0 everywhere (reviewer MEASURED on a synthetic case: `copy_raw` 1, `copy_filtered` 0, `frac_with_step` 0.0). Chapter A shows raw and filtered totals side by side.
- The real-data run (brief step 5) is deferred. Do not read PROVISIONAL outputs as results.

## Absence check (Task 8)

Code: `pathway_census/absence.py`, `bin/absence_check.py`, `slurm/absence_check.sbatch`. Tests: `tests/test_absence.py`; fixtures from the seeded `tests/data/make_absence_fixtures.py`.
Tool: `tblastn` from `module load ncbi-blast/2.17.0+` (BLAST 2.17.0+, build Jul 1 2025; MEASURED with `tblastn -version`, 2026-11-20). Real tblastn is used in tests; nothing is mocked.
- `select_absences(matrix, relatives)`: copy_raw == 0, has_proteins True, and at least one listed relative with copy_raw > 0. `has_proteins` False genomes are not selected (their zero is not evidence). `unevaluated_genomes(matrix)` lists them.
- `same_genus_relatives(genomes)`: groups by genus only where `genus_source != "unresolved"`; never self; single-genome genus gives `[]`.
- `check_absence(genome_fna, step_refs, min_bitscore, min_query_cover=0.7)` returns `present_unannotated | fragmented | absent | unconfirmed`. Rules and out-of-scope cases (stop codon in frame) are in the module docstring.
- PROPOSED, no measured data: `min_query_cover = 0.7`, `min_bitscore = 50`.
- `downloaded == False` gives `unconfirmed`, reason `no_bfd_fasta`. The CLI does not look for a FASTA elsewhere (for example NCBI_ASM). The FASTA path comes from `genomes.scaffolds_path` applied to the proteins path in `seqdb/genomes.tsv`.
- CLI refuses `verified: false` without `--provisional`; outputs then carry `.PROVISIONAL`.
- SLURM sizing is UNVERIFIED (target 1-1.5 h per job; no job has run). Chunk with `PC_LOCUSTAGS`.
- Chunk outputs are named `<pathway>_absence.chunk<ID>[.PROVISIONAL].tsv` and `<pathway>_unevaluated_genomes.chunk<ID>...tsv`; ID = `SLURM_ARRAY_TASK_ID`, else the `PC_LOCUSTAGS` basename minus its last suffix. Merge them with `bin/merge_absence.py --out results/carotenoid_absence_summary.tsv <dir>/carotenoid_absence.chunk*.tsv` (output name PROPOSED; Chapter A reads it). The script refuses mixed provisional and non-provisional chunks (a `.PROVISIONAL` name or a `PROVISIONAL` header line marks a chunk provisional), chunks with different `#` header lines or columns, and a (locustag, step_id) pair in two chunks. Provisional chunks give `carotenoid_absence_summary.PROVISIONAL.tsv` with a `# PROVISIONAL` first line. Chapter A refuses a file with a `PROVISIONAL` header line even under the canonical name. `unevaluated_genomes` is the same list in every chunk (it is not chunked): keep one copy.
- `evalue` (tblastn, default 1e-3) is PROPOSED. In-frame stop: one synthetic case MEASURED as present_unannotated (see module docstring).
- Real-data run (brief step 5) is DEFERRED by controller ruling. No counts per outcome exist.

## Paired deposited-versus-BFD comparison (Task 9)

Code: `pathway_census/paired.py`, `bin/paired_compare.py`. Tests: `tests/test_paired.py`. Details are in the `paired.py` docstring.

What the paired test shows: for assemblies that have both a deposited annotation and a BFD annotation, it shows how many more or fewer copies of each pathway step the BFD annotation gives. The assembly is the same in both arms. So the difference is the effect of the annotation.

What it does NOT show: it says nothing about the genomes with no prior annotation (BFD-scoped MEASURED, Task 2: 18,055 `no_prior_annotation_present` and 673 `no_prior_annotation_not_downloaded`; spec 5.6 gives the TSV-wide 18,145). It does not show which annotation is correct. It does not show absence of a gene. A real-data difference counts as a finding only if it has the same sign across annotation providers and survives the BUSCO control. Otherwise report it as not established.

- Pair rule: `prior_status == baseline_ok`, BFD proteins exist, deposited proteins flag true, the deposited protein file exists and was searched. Other genomes are counted in `n_unpaired` under the first failing reason: `not_baseline`, `no_bfd_proteins`, `no_deposited_proteins`, `deposited_file_missing`.
- Keys: both matrices use the caller-supplied `locustag`. Deposited protein IDs are not parsed. Counts are kept per source and never merged by protein ID.
- Run the same HMMs and the same assignment on both sources. Build the deposited matrix with a genome table copy where `has_proteins` = the deposited-file flag. The deposited file path is `paired.deposited_protein_path(root, folder)`: first `<root>/<folder>/<folder>_protein.faa.gz`; `<root>/<folder>/<folder>_protein.faa` only if the `.gz` file does not exist and the plain file does. `folder` is the prior-annotation TSV `folder` column, carried by `load_genome_table`. `count_pairs(genomes, deposited_root)` needs the `folder` column; a null folder counts as `deposited_file_missing`.
- Spec vs measured layout: spec 5.6 said `<accession>/protein.faa`. The MEASURED layout differs (final review, reviewer): 0 of 24,807 folders in NCBI_ASM are named by bare accession, and in the first 200 GCF folders 0 have `protein.faa` and all have `<folder>/<folder>_protein.faa.gz`. Only the first 200 GCF folders were checked for the file name.
- `paired_table`: `delta = copy_bfd - copy_deposited` (copy_raw). `delta_filtered` is set only when both unknown counts are 0, else NA. Unknown counts are kept in `unknown_deposited` / `unknown_bfd`.
- `summarize_paired`: bootstrap over genomes (resample genomes with replacement), 95% interval from the 2.5 and 97.5 percentiles. Rows per (step, provider) plus one `ALL` row per step. Providers below `min_pairs` are returned in `dropped_providers`, not removed silently. PROPOSED defaults, no measured basis: `min_pairs=20`, `n_boot=2000`, `seed=20261120`.
- `summarize_paired_by_busco`: mean delta per step per BUSCO bin, default edges 0, 80, 95, 100 (PROPOSED), plus an explicit `unknown` bin for null BUSCO. A step with zero copies in both sources appears with mean 0.
- CLI refuses `verified: false` without `--provisional`; outputs then carry `.PROVISIONAL`.
- DEFERRED by controller ruling: the real pair count (brief step 1; Fable's earlier figure of 5,625 `baseline_ok` genomes is the genome-table count, not a pair count, and was not re-run) and the real-data run (brief step 6). No real counts or deltas exist.

### DEFERRED follow-up (G6): the deposited arm of the real-data pass

Not implemented. `paired.py` compares two finished matrices only. Building the deposited matrix is open work.
- Input: `paired.deposited_protein_path(deposited_root, folder)` for each pair genome (deposited root `/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM`; `folder` from the genome table). A table `protein_id -> locustag` (the caller supplies it; deposited protein IDs do NOT start with the LOCUSTAG prefix before the first underscore (8 hex characters in all but 3 BFD genomes)).
- Output: a matrix in `build_matrix` format (same pathway steps, same HMMs and assignment rules as the BFD arm), keyed by the BFD `locustag`, with `has_proteins` set from the deposited-file flag (`has_deposited_proteins` and the file exists).
- Design choice, to decide in the real-data pass (two options): (a) rewrite deposited protein IDs to `<locustag>_<deposited id>` before assign / quality / matrix; or (b) add a `genome_of` mapping parameter (protein_id -> locustag) to `matrix.build_matrix` and `quality.quality_flags`. Option (a) needs no code change but alters protein IDs in every output; option (b) keeps IDs and changes two functions. Neither is chosen.
- Change surface (places that take the genome from the `LOCUSTAG_` prefix of a protein ID, i.e. the text before the first underscore; it is 8 hex characters in all but 3 BFD genomes (MEASURED: 3 of 24,363 LOCUSTAGs have 9 characters: FF38F669A, F3DC5624A, FCB937CBA); found by grep on 2026-11-20):
  - `pathway_census/matrix.py:77` `a["locustag"] = a.protein_id.str.split("_").str[0]`, and the error text at line 28.
  - `pathway_census/quality.py:104` `{p.split("_")[0] for p in assigned.protein_id}`; also the gene-id rule in the quality.py docstring (line 14: `F77FFF1D_000001-T1` -> `F77FFF1D_000001`, `gene_info` keyed by gene id), which a deposited ID format may not follow.
  - `pathway_census/assign.py`: no prefix parsing found. `absence.py` and the CLIs `bin/*.py`: no prefix parsing found (grep for `split("_")`, `8-hex`, `[:8]`). UNVERIFIED beyond that grep; the real pass must re-check.
- `load_genome_table` now carries `asmid` and `folder` (final-review fix I7). `results/genomes.tsv.zst` must be regenerated before the real pass to get these columns.

### Reserved values, flags and outputs
- Providers named `ALL` and `(missing)` are reserved. `summarize_paired` raises ValueError if a real provider has either name.
- A genome flag that is NA counts as not True. NA and False share one reason code. `n_unpaired["unknown_flag"]` counts the NA cases; it is a subset of the reason counts, not an extra reason.
- Bootstrap resampling is per step, not joint across steps. `dropped_providers` has column `provider`; the summary has `annotation_provider`.
- Every output file carries `# params: min_pairs=... n_boot=... seed=... busco_bins=...` (the JSON has a `params` key). `--busco-bins` sets the edges.

## Triage bins (Task 10)

The bins are triage. They are not a ranking. They are not an activity prediction (spec 5.5). No weights are fitted and there is no additive score.
The rules are written in the docstring of `pathway_census/bins.py`, and they were fixed before any real data was binned. Summary:
- Bins: `tier1`, `tier2`, `tier3`, `excluded`, `unknown`.
- `excluded`: a required flag is known bad (truncated True, contam_suspect True, domain_ok False, motif_ok False). `unknown`: no flag known bad, at least one required flag NA.
- Four boolean terms: `divergent`, `lineage_expansion`, `in_cluster`, `expressible`. tier1 = at least `tier1_min` true, tier2 = at least `tier2_min`, else tier3.
- `hap_dup_suspect` True demotes tier1 to tier2 with a reason. NA does nothing and is shown.
- All parameters are PROPOSED (no data basis). `window_bp=20000` is UNVERIFIED against the known cluster span in the reference species.
- CLI `bin/triage.py --pathway P.yaml --copies copies.tsv[.gz|.zst] --out triage.tsv`. Refuses `verified: false` without `--provisional` (output then gets `.PROVISIONAL`). Booleans are strict (`True`/`False`/empty). Output has `# params:` as its first line (then a PROVISIONAL line if tagged): read it with `pd.read_csv(..., sep="\t", comment="#")`. A `.zst` out name is compressed.
- Input is an already-joined frame (columns in the `bins.py` docstring). `signalp` is positives only: the caller passes `signal_peptide=False` for proteins absent from the table (ASSUMPTION, unverified). `soluble_steps` is a parameter because the pathway YAML has no soluble flag.
- DEFERRED by controller ruling: the real-data run (brief step 6). No real bin counts exist.

## Batched BFD search (census steps 2 to 4)

Scripts: `bin/score_refs_decoys.py`, `bin/make_min_scores.py`, `bin/make_chunks.py`, `bin/search_chunk.py`.
SLURM: `slurm/score_refs_decoys.sbatch`, `slurm/search_chunk.sbatch`.
Tests: `tests/test_make_min_scores.py`, `tests/test_make_chunks.py`, `tests/test_search_chunk.py`.
The census plan is `docs/superpowers/plans/2026-11-22-census-step-list.md`.
All scripts run under `/usr/bin/python3.12`. A SLURM job reads its paths from environment variables. No script uses `BASH_SOURCE`.

### Run order for one pathway

1. `bin/build_step_hmms.py` builds the step HMMs (`--provisional` while `verified: false`).
2. `bin/score_refs_decoys.py` (or `sbatch slurm/score_refs_decoys.sbatch`) scores every reference and decoy against each step HMM, with a leave-one-out score per reference. Output in `results/step2/<pathway>/`: `step_scores*.tsv`, `leave_one_out*.tsv`, `targets_at_or_above_min*.tsv`.
3. `bin/make_min_scores.py` turns `step_scores` into `min_scores*.tsv`.
4. `python -m pathway_census.genomes` writes `results/genomes.tsv.zst`.
5. `bin/make_chunks.py` writes the chunk manifest.
6. One `sbatch slurm/search_chunk.sbatch` per chunk runs `bin/search_chunk.py`.

### `make_min_scores.py`: the `min_score` rule (owner, 2026-11-22)

- A step with two or more references uses its lowest leave-one-out score.
- A step with one reference uses its highest decoy score as a floor.
- A step without a value for its basis stops the script with an error. There is no default.
- `--scale STEP=FACTOR` (repeatable, 0 < FACTOR <= 1) lowers one step for an exploratory run. The column `rule_min_score` keeps the unscaled value. A search keeps every hit at or above `min_score`, so a lowered value gives a superset. Filter on `score` to get the rule result back.
- Output columns: `step, min_score, basis, n_refs, rule_min_score`.

### `make_chunks.py`: the chunk manifest

- Input: `results/genomes.tsv.zst` and `seqdb/genomes.tsv` (locustag to proteins path).
- A genome is listed only if `has_proteins` is true (case-insensitive) and the seqdb has a path. Any other value is excluded. Nothing is guessed.
- Excluded genomes go to `--excluded` with a reason. They are never absences.
- The order is shuffled with `--seed` (default 1), so every chunk has a mixed set of genomes. The same seed gives the same chunks, whatever the input order.
- Output: `chunk_id, locustag, proteins_path` (zstd). A size of 4000 gives 6 chunks for the 23,376 genomes listed on 2026-11-22.

### `search_chunk.py`: one chunk, all step HMMs

- It checks first, and fails before any staging work, for: `SCRATCH` unset, a step in `min_scores` without an HMM, and an unverified pathway without `--provisional`.
- It joins the protein files of the chunk into one FASTA in `$SCRATCH/chunkNNNNN/`, then runs one `hmmsearch` per step.
- It counts proteins per genome and counts ids that do not start with `<locustag>_` (`n_bad_id_prefix`). BFD ids carry the locustag, so ids are unique across genomes. A non-zero count needs a look before the candidate table is used.
- A missing, empty or truncated protein file gives `status` `missing` for that genome. The chunk goes on. A `missing` genome is never an absence.
- Outputs in `--out-dir`, names carry the chunk id:
  - `candidates.chunkNNNNN[.PROVISIONAL].tsv.zst`: one row per (protein, step) at or above `min_score`, columns of `hits.COLUMNS` plus `step_id`. The first line is `# PROVISIONAL` for a provisional run.
  - `timing.chunkNNNNN.tsv`: seconds per stage, genome, protein and candidate counts.
  - `searched.chunkNNNNN.tsv`: the searched-genome manifest. It is written last. A chunk without this file is not complete and must be re-run.
- Re-run rule: a chunk is complete if and only if its `searched.chunkNNNNN.tsv` exists. To re-run a chunk, delete its three output files and submit it again. Chunks are independent.
- Resource numbers are MEASURED on one chunk of 200 genomes (16 CPUs, 1 min 31 s, 1.48 GB). They are UNVERIFIED for a 4,000-genome chunk until the first full chunk reports.

### Is Nextflow needed yet?

Not yet (PROPOSED, no data against it). The run is 6 independent jobs per pathway. A chunk is complete when its manifest file exists, so a loop over missing manifests is a safe re-run. The later steps (trees, placement, matrices) have dependencies that Nextflow would handle well. Revisit when the pipeline has more than about 20 job types, or when it runs for several pathways at once and a failed job needs automatic retry. The scripts keep paths in arguments and write complete-marker files, so a later Nextflow wrapper does not need changes to them.

## Step 5 pilot (phytoene_synthase; step-generic)

Scripts: `bin/step5_place.py` (subcommands `extract`, `select-reps`, `clades`, `placements`, `summarize`; the step id is an argument). SLURM: `slurm/step5_prepare.sbatch`, `slurm/step5_place.sbatch`. Tests: `tests/test_step5_place.py`. The tools hold no step name. `lycopene_cyclase` runs with `PC_STEP=lycopene_cyclase` and `PC_OUTGROUP_FA=<outgroup FASTA>` (an outgroup file exists only for that step). Reused without edits (imported from `analysis/enzymes_UPO`): `upo_place.py` (`assign`, `match-cols`, `edge_bipartitions`, `to_newick`) and `refset/bin/tree_outline.py` (`parse_newick`, `reroot`).

### Run order
1. `PC_STEP=<step> sbatch slurm/step5_prepare.sbatch` (run from the worktree root). It cuts domains, clusters (mmseqs2 `easy-cluster`, `--min-seq-id 0.5 -c 0.8 --cov-mode 0`), draws representatives, aligns (mafft `--auto`), trims (clipkit `smart-gap`), builds the tree (raxml-ng, `LG+G4`, 5 parsimony + 5 random starts, seed 1) and writes `clade_members.tsv`.
2. Read `clade_report.json`. If `refs_form_clade_free_of_decoys` is false, stop: that is a result (`step5_place.py clades` exits with code 3).
3. `PC_STEP=<step> PC_MIN_SUPPORT=<x> sbatch slurm/step5_place.sbatch`. It builds the reference HMM (`hmmbuild --symfrac 0`), aligns all domains (`hmmalign`, insert columns dropped), places with epa-ng (pixi, pinned 0.3.8), applies the `upo_place assign` rule, writes `placements.tsv`, runs `bin/assign_steps.py --provisional` and writes `SUMMARY.tsv`.
Environment variables: `PC_STEP` (required), `PC_PATHWAY` (carotenoid), `PC_ROOT`, `PC_SEQDB`, `PC_MAX_REPS`, `PC_OUTGROUP_FA`, `PC_MIN_LWR`, `PC_MAX_PENDANT`, `PC_MIN_SUPPORT` (required for part 2, no default). Output folder: `results/step5/<pathway>/<step>/` (git-ignored). `RUN.md` in that folder lists the exact commands and tool versions.

### Inputs
- Candidates: `results/search/carotenoid_full/candidates.chunk0000{0..5}.PROVISIONAL.tsv.zst`, rows with `step_id == <step>`.
- References, decoys, outgroup: `refs/carotenoid/<step>.faa`, `refs/carotenoid/decoys/<step>.faa`, and (lycopene_cyclase only) `refs/carotenoid/outgroups/<step>.faa`.
- Proteins: `seqdb/proteins/<sha1(LOCUSTAG)[:2]>/<LOCUSTAG>.proteins.fa.gz` (bgzip with `.fai`). MEASURED: the layout holds for the 10,047 genomes needed; 0 of 10,863 proteins failed to extract.
- Phylum per genome: `results/genomes.tsv.zst`.

### Rules and parameters
- Domain cut: best row per protein (highest score), `ali_from..ali_to` padded by 20 residues each side, clamped to the protein ends. The pad of 20 is PROPOSED. Coordinates are in `domain_coords.tsv`.
- Representatives: one per mmseqs cluster. Cluster threshold 0.5 identity and 0.8 coverage are the values from the brief (PROPOSED). Cap 300 (PROPOSED): if there are more clusters, draw 300 with `random.Random(1)`, proportional by phylum, at least 1 per phylum (largest-remainder allocation, never more than the stratum size). The cap did not apply in the phytoene_synthase run (191 clusters).
- Rooting: the outgroup file if given, else the decoys (the smallest bipartition side that holds them all). epa-ng takes an unrooted tree; the clade rule uses bipartitions and does not depend on the root.
- Clades (`clade_members.tsv`): `target` = the smallest bipartition side that holds all references and no decoy or outgroup, with the representatives inside it. `decoy` = every leaf that is not in `target` and not an outgroup (decoys plus the representatives outside `target`). `outgroup` = the outgroup leaves. A placement counts for a clade only when its edge lies inside it (the `upo_place assign` rule). Representatives are members because a query placed next to a representative must count for that clade.
- `--min-lwr 0.5` (PROPOSED, taken from upo_place, not validated for this step), `--max-pendant 1.0` (flag only; PROPOSED, from upo_place). `PC_MIN_SUPPORT=0.5` for `assign_steps` (PROPOSED; with `--min-lwr 0.5` it removes nothing more).
- `placements.tsv` writes `clade = target` only when the call is `target`. A `decoy`, `outgroup` or `unassigned` call gives an empty clade, because `assign_steps` treats any non-empty clade as in-step. `placements_detail.tsv` keeps the real call, the support, the stem LWR, the pendant length and `long_pendant`.
- `assign_steps` cannot be limited to one step by an option. It was run on a candidates table that holds only this step (`candidates_step.tsv`) and the placements of this step.

### Re-run rule
Part 1 is complete when `clade_report.json` exists. Part 2 is complete when `SUMMARY.tsv` exists. To re-run a part, delete its outputs (or set `PC_OUT` to a new folder) and submit again. Part 2 needs part 1 from the same folder. Changing `PC_MAX_REPS` or the seed changes the tree: re-run both parts. Both jobs ran on the short queue (MEASURED for phytoene_synthase: part 1 31 min 43 s, 16 CPUs, of which raxml-ng used 30 min with 4 threads; part 2 5 min 0 s).

### Reading the phytoene_synthase result (MEASURED, 2026-11-22)
See `results/step5/carotenoid/phytoene_synthase/SUMMARY.tsv`. The clade `target` holds 186 of the 200 leaves. It is the smallest clade with all 5 references, and it contains 181 of 191 representatives from 5 phyla. So the placement separates the candidates from the 4 squalene synthase decoys only. It does not separate real phytoene synthase from other family members inside `target`. No candidate domain has 50% identity or more to any squalene synthase decoy, and the step 2 maximum decoy score was 49.3, below the 98.8 cutoff. So squalene synthase is not among the 10,863 candidates (MEASURED); the 155 candidates called `decoy` sit on a branch of 10 representatives whose best hit is a phytoene synthase reference, not an ERG9 decoy. What the 10,863 candidates are is UNVERIFIED.

### v2: named tree members, non-fungal outgroups, lineage subclades (phytoene_synthase, 2026-11-22)

New options. Every default is unchanged. MEASURED on the v1 tree: `make-refset` and `clades` with the v1 inputs write byte-identical `ref_set.faa`, `clade_members.tsv` and `clade_report.json`. The v1 folder `results/step5/carotenoid/phytoene_synthase/` was not touched. The v2 folder is `results/step5/carotenoid/phytoene_synthase_v2/` (`RUN.md` has commands, tools, jobs and failed submissions).
- `step5_place.py make-refset`: writes the reference-alignment FASTA from references, decoys, named members, outgroups and representatives, and stops on a duplicate id. The cap of 300 applies to candidate representatives only. Members, outgroups, references and decoys are never in the draw.
- `step5_place.py clades --members-fa F --outgroup-fa F --outgroup-role tree_outgroup --labels labelled_leaves.tsv`. `--outgroup-fa` (alias of `--outgroup`) may hold several sequences; the tree is rooted on their MRCA (the smallest side that holds them all). The report gives `outgroups_monophyletic`, the roles inside the outgroup MRCA, and where the decoys sit. A named member that is also a reference, decoy or outgroup stops the run. The `role` column has `named_member` and (with `--outgroup-role tree_outgroup`) `tree_outgroup`.
- Named members are never decoys: a member outside `target` gets no clade label, and it is left out of the `decoy` set. They do not enter the cutoff or HMM code. The `target` definition is unchanged (smallest clade with all references and no decoy or outgroup); the report adds `target_named_members` and `named_members_outside_target`.
- `step5_place.py make-labels`: `labelled_leaves.tsv` (leaf, role, class, evidence) from the references (class from a TSV `accession<TAB>class`, a required input: nothing is guessed) and the members (class from `tree_members/provenance.tsv`).
- Lineage subclades (PROPOSED rule): for each class with at least 2 labelled leaves (references plus named members), take the smallest clade that holds them all. It is a named subclade `<Class>` only if it holds no labelled leaf of another class and no decoy or outgroup leaf. Otherwise no subclade is made and `lineage_subclades.tsv` gives the reason. The subclade labels go in the `fig1_clade` column of `clade_members.tsv` (the second clade system of the `upo_place assign` rule), so one `assign` run gives both the target call and the subclade call. `placements --system fig1` writes the subclade as the clade call (`placements_subclade.tsv`); `--system our` (default) keeps `target`.
- `step5_place.py subclade-report`: `subclade_counts.tsv`, `nearest_labelled.tsv` (nearest reference or named member by tree distance for each representative, with the number of candidates in its cluster) and `nearest_labelled_by_class.tsv`.
- Env vars for `slurm/step5_prepare.sbatch`: `PC_MEMBERS_FA`, `PC_MEMBERS_PROV`, `PC_REF_CLASSES`, `PC_OUTGROUP_ROLE`, `PC_RAXML_THREADS`, `PC_RAXML_WORKERS`, `PC_RAXML_STARTS`. Both sbatch files now carry `--constraint="ryzen|cascade|broadwell"`: the abu_dhabi nodes (`c01`-`c29`) gave `Illegal instruction` in `mmseqs easy-cluster` (MEASURED, job 32377314). `sbatch --export` splits at commas: pass a value with a comma through the shell environment with `--export=ALL`.
- raxml-ng 2.0.4 aborted with "core oversubscription" at `--threads 8` (jobs 32377295 and 32377308; 16 CPUs = 8 physical cores). `--threads 4` works.

v2 result (MEASURED; `SUMMARY.tsv`, `clade_report.json`). Inputs: 220 leaves = 5 references, 4 decoys, 13 named members, 7 outgroups, 191 representatives. Trimmed alignment: 14.65% of columns removed by clipkit. Tree: LG+G4, best log-likelihood -113638.6 from 3 searches (only 3 starting trees ran; see `RUN.md`).
- The 7 outgroups are not monophyletic. Their MRCA side has 11 leaves: the 7 outgroups and the 4 squalene synthase decoys. The 4 decoys form a clade of their own (4 leaves).
- `target`: 185 leaves = 5 references, 13 of 13 named members, 167 representatives. `decoy` side: 28 leaves (4 decoys, 24 representatives). v1 had 186 / 181 / 14 / 10.
- Subclades (PROPOSED rule): Dothideomycetes 3 labelled, 0 representatives, 813 candidates; Eurotiomycetes 4 labelled (34 leaves), 30 representatives, 1,376 candidates; Mucoromycetes 3 labelled (10 leaves), 7 representatives, 232 candidates; Ustilaginomycetes 2 labelled (4 leaves), 2 representatives, 61 candidates. No subclade: Sordariomycetes (5 labelled; the smallest clade that holds them also holds the 3 Dothideomycetes members B2WAQ3, E4ZUB5, Q0V6M5) and Tremellomycetes (1 labelled leaf).
- Candidates (10,863): `target` 10,815 and `decoy` 48 (v1: 10,708 and 155). Of the 10,815 target calls, 2,482 fall in a named subclade and 8,333 in none. No subclade call lies outside `target`. Support (LWR sum in target): minimum 0.670, 10,814 of 10,815 at least 0.7, 10,812 at least 0.9. long_pendant 0. assign_steps (`--min-support 0.5`): 10,815 assigned, 48 unassigned, 10,003 genomes with an assigned protein.
- So the named subclades hold 23% of the candidates. 77% of the candidates sit in `target` with no named subclade. The 12 of 13 named members are annotated only by similarity to P37295, so a subclade name here is a tree grouping of named leaves, not evidence of function for a candidate.

## Gain and loss on the genus tree

Question: are the carotenoid steps gained or lost on the same branches of the fungal tree? Code: `bin/make_calls_table.py` (calls table; a TSV config switches a step between `hits` and a placed `assigned.tsv`), `bin/gainloss_lib.R`, `bin/run_gainloss.R`, `bin/run_pagel.R`, `bin/gainloss_compare_trees.R`, `slurm/gainloss_run.sbatch`, `slurm/gainloss_pagel.sbatch`; figures `reports/carotenoid/figures/bin/fig9_tree_presence.R`, `fig10_gain_loss_counts.R`, `fig11_co_events.R`. Tests: `tests/test_make_calls_table.py`, `tests/test_gainloss_r.py` (runs `tests/gainloss_toy_check.R`; skipped when Rscript cannot load ape and castor; use `module unload R; module load R/4.2.0`). Full record, commands, decisions and failures: `results/gainloss/carotenoid/RUN.md`.

Summary of what was measured (calls are PROVISIONAL; 1,639 of 1,654 genus representatives are tree tips with data, 4 are NA, 11 are not in the tree):

- Tree: FastTree support tree (1,643 tips, substitutions per site, unrooted) was chosen over the ASTRAL consensus (rooted, coalescent-unit internal lengths, lower Mk likelihood for all four steps, 70.9% of nodes with support at or above 0.8 against 93.9%). The root (Opisthosporidia split) is PROPOSED. Changing the root to Blastocladiomycota or Chytridiomycota changed event counts by at most 3 branches.
- Steps: phytoene_synthase and lycopene_cyclase (placed, with named members), phytoene_desaturase and cart (placed on small reference trees: 15 and 27 leaves, so a weaker filter). ARD fits better than ER for every step (delta AIC -76 to -167). With a 0.8 probability rule, gains per step are 5 to 13 and losses 39 to 117 (ARD); most events are on tip edges.
- Co-events: phytoene_synthase, phytoene_desaturase and lycopene_cyclase share 7 to 9 gain branches and 65 to 71 loss branches (expected under the permutation null: under 0.1 and 5 to 7; q at the floor of 0.0004). Most shared events are tip edges, the null ignores tree structure, and synthase and cyclase are one fused protein in most fungi. cart shows no co-gain and at most weak co-loss. Pagel tests give p below 1e-8 for all 6 pairs; this does not separate shared events from shared clade structure.
- Not done: null simulations, a completeness filter, ggpp_synthase and floor steps. See RUN.md.

### Follow-up: BUSCO filter, simulation null, ASTRAL comparison

Record and tables: `results/gainloss/RUN_followup.md`. Code: `bin/filter_calls_busco.py` (tips below a BUSCO cutoff set to NA for all steps; also a lineage-relative option), `bin/gainloss_null_lib.R`, `bin/run_null_sim.R`, `bin/summarise_null.R` (parametric bootstrap null for co-events, with an optional prevalence-matched variant), `bin/run_astral_compare.R`, `bin/compare_busco_versions.py`, `slurm/gainloss_null.sbatch`, figure `reports/carotenoid/figures/bin/fig13_null_co_events.R`. Tests: `tests/test_filter_calls_busco.py`, `tests/test_gainloss_null.py`.

- BUSCO cutoff 70 (primary) removes 69 of 1,654 representatives; 65 removes 59. Core loss counts change by 3% to 13% and core co-loss counts by 6% to 14% against the unfiltered run. 31 of the 69 removed tips are Microsporidia (30) and Rozellomycota (1). A 65 to 70 floor cannot catch an incomplete yeast assembly (Saccharomycetes sit at 72 to 78).
- The simulation null as specified (fitted ARD rates, stationary root, N = 2,000) simulates far fewer present tips than the data have (mean 121 to 150 against 683 to 737 for the core steps) and is not a fair null. A prevalence-matched null (present-tip count within 10% of observed) reproduces the observed gains and losses. Against it, core co-loss is 59 to 71 observed against about 7 to 10 expected (internal edges: 15 to 17 against about 1), p = 0.0005 (floor, N = 2,000); cart pairs are within chance (p 0.26 to 0.54). Synthase x cyclase is one fused gene and not a meaningful test.
- The ASTRAL tree gives similar core co-event counts (co-loss 57 to 61 against 65 to 71) with its own branch lengths; with all lengths set to 1 the counts change more.

### New primary run: Microsporidia excluded, Blastocladiomycota root

Record, tables and the side-by-side comparison: `results/gainloss/carotenoid_noMicro_blasto/RUN.md` and `results/gainloss/carotenoid_noMicro_blasto_compare/`. Code: `bin/drop_phylum_calls.py`, `bin/prune_tree_tips.R`, `bin/run_gainloss_roots.R`, `bin/compare_runs.py`, `bin/events_on_taxon_edges.R`; figure scripts `fig9_tree_presence_nomicro.R`, `fig10_gain_loss_counts_nomicro.R`, `fig11_co_events_nomicro.R`, `fig13_null_co_events_nomicro.R`. The 30 Microsporidia tips are dropped and the tree is rooted on the 6 Blastocladiomycota tips (they form a clade of exactly 6 tips). None of the old run's events sat on Microsporidia edges; the counts change by a few branches through the model refit. Core co-loss 62 to 69 and co-gain 8 to 9 still exceed the prevalence-matched simulation null (internal-edge co-loss 13 to 16 against about 1.2); cart pairs are within chance.

## Step 5 for the single-reference steps

Record with commands, floors, jobs and failures: `results/step5/carotenoid/single_ref_steps_RUN.md` (2026-11-23, PROVISIONAL). Run folders: `results/step5/carotenoid/{cars,cars_v2,cars_v2best,acaa,acaa_best,tsp3,tsp3_v2,tsp2,tdh,tdh_id80}/`. Inputs, screens and submit scripts: `results/step5/carotenoid/single_ref_inputs/`. carx and crtr were not placed (see the record).

New tools (existing step 5 tools are unchanged):
- `bin/filter_candidate_folder.py --in-dir D --out-dir D2 --step S --floor F`: candidate folder of one step at a floor (same file format). Folders `results/search/<step>_floor<N>/`.
- `bin/screen_outgroups.py` and `bin/select_by_screen.py`: alignability screen (best domain >= 10 bits and covering >= 0.5 of the HMM; `best` and `union` coverage) and selection of the proteins that pass.
- `bin/step5_single_ref_clades.py` with `slurm/step5_single_ref_clades.sbatch`: `target` rule for ONE reference. The default rule (smallest side with the references) gives the reference leaf alone, so no placement can be `target`. New rule (PROPOSED): the largest side with the reference and no decoy, outgroup or named member. Run it between `step5_prepare.sbatch` and `step5_place.sbatch`. It keeps the default files as `clade_members.smallest.tsv` and `clade_report.smallest.json`.
- `slurm/step5_prepare_ident.sbatch`: `step5_prepare.sbatch` with `PC_MIN_SEQ_ID` (clustering identity, default 0.5).
- `bin/single_ref_nearest.py` (nearest labelled leaf, distance to the reference, genus and class, calls by class), `bin/target_overlap.py` (overlap of target sets and of hits), `bin/hits_by_class_floors.py` (share of genomes with a hit per class).
- Both `step5_prepare.sbatch` and `step5_place.sbatch` read `PC_CAND_DIR`. Export it for both parts. After each run, `candidates_step.tsv` must have as many rows as `placements.tsv`.
- Reading rule: with this target rule, `target` is bounded by the labelled leaves that were put in the tree. It is not a set of orthologs. A larger target can follow from fewer labelled neighbours.
Tests: `tests/test_single_ref_tools.py`.

### Fifth step: ggpp_synthase (placed) on the new primary

Record: `results/gainloss/carotenoid_noMicro_blasto_5steps/RUN.md`. ggpp_synthase is present in 1,540 of 1,609 tips with data (96%); the model finds 0 gains, 60 losses (57 on tip edges), 0 uncertain branches, and ER fits as well as ARD. Its co-loss with the other steps (3 to 4 branches, all tip edges) and its co-gain (none) are within chance in the tree-aware null. Absent calls are mostly single tips; grouped absences are in Pichiomycetes and Agaricomycetes. Code: `bin/gainloss_null_models.R`, `bin/run_null_sim_models.R`, `bin/summarise_null_models.R`, `bin/ggpp_absence_analysis.R`.

## Steps 1 and 2 for melanin and mycosporine

Note with the sign-off table and owner decisions: `docs/superpowers/notes/2026-11-23-melanin-mycosporine-step2.md`. Branch `pathway-census-melanin-myco`. All outputs are PROVISIONAL and in git-ignored `results/`.

Commands (run from `analysis/pathway_census`; `<p>` is `melanin_dhn`, `melanin_ldopa` or `mycosporine`):

```
module load hmmer/3.4
module load mafft/7.505
/usr/bin/python3.12 bin/build_step_hmms.py --pathway pathways/<p>.yaml --refs-dir refs/<p> --out-dir results/hmms/<p> --provisional
# mycosporine: add  --steps ddgs omt atp_grasp_ligase  (two steps have no reference)
PC_ROOT=$PWD PC_PATHWAY_ID=<p> PC_OUT_DIR=$PWD/results/step2/<p> PC_TAG=PROVISIONAL \
  sbatch --constraint="ryzen|cascade|broadwell" slurm/score_refs_decoys.sbatch
/usr/bin/python3.12 bin/make_min_scores.py --step-scores results/step2/<p>/step_scores.PROVISIONAL.tsv --out results/step2/<p>/min_scores.PROVISIONAL.tsv
sbatch -p short --constraint="ryzen|cascade|broadwell" --cpus-per-task=4 --mem=8G --time=00:20:00 --wrap="module load hmmer/3.4; /usr/bin/python3.12 bin/target_scores_at_min.py --refs-dir refs/<p> --hmm-dir results/hmms/<p> --min-scores results/step2/<p>/min_scores.PROVISIONAL.tsv --out results/step2/<p>/target_scores_at_min.PROVISIONAL.tsv --tag PROVISIONAL"
```

Results (MEASURED): 18 HMMs (9 + 6 + 3). Step 2 jobs 32407157 (dhn), 32407158 (ldopa), 32407159 (mycosporine); add-on jobs 32407363, 32407364, 32407365.

- `bin/target_scores_at_min.py` (new, `tests/test_target_scores_at_min.py`) lists every decoy and other-step reference at or above the rule `min_score`. `score_refs_decoys.py` counts at the lowest full-HMM reference score, so its `targets_at_or_above_min` table is empty for these pathways. Use `target_scores_at_min.PROVISIONAL.tsv`.
- `melanin_ldopa`: `hppd` and `aspulvinone_e_synthase` are single-reference steps with no decoy, so `make_min_scores.py` stops with an error. `min_scores.PROVISIONAL.tsv` was made from `step_scores.floor_defined.PROVISIONAL.tsv` (the same table without those two rows).
- `mycosporine`: `mg_glucosyltransferase` and `glutamicol_formation` have no reference, so no HMM.
- No search, no tree, no lowered cutoff (`--scale` not used).

Follow-up 2026-11-23 (note section 9): decoys were added for `hppd`, `aspulvinone_e_synthase`, `ayg1`, `scytalone_dehydratase`, `hmga` with `curation_scripts/melanin_myco_decoys/add_decoys.py`. `hppd` and `aspulvinone_e_synthase` now have a floor, so `step_scores.floor_defined` is no longer needed. New `make_min_scores.py` options: `--floor-override STEP=VALUE` with `--override-basis TEXT` (mycosporine `omt` uses `--floor-override omt=100.6 --override-basis "owner override 2026-11-23: highest decoy score; revisit in second round"`), and `--allow-undefined` (empty `min_score` and the basis "no floor defined: ..." for a step with no rule value). Run `score_refs_decoys.sbatch` from a directory you can write to (it writes its log to the current directory).

## Step 3 for melanin and mycosporine

Record: `results/step3/melanin_myco/RUN.md` (commands, jobs, timings, candidates per step, projection). Note: `docs/superpowers/notes/2026-11-23-melanin-mycosporine-step2.md`, section 11. All outputs are PROVISIONAL. Step 3 is one timing job on 200 genomes. The full search was not run.

- No change to `bin/search_chunk.py`. It takes one `--hmm-dir` and one `--min-scores` table, so the 18 HMMs of the three pathways go in one directory (symlinks) with one combined `min_scores` table (the three step 2 tables joined; step names are unique). `--pathway` only gates on `verified`.
- Timing job 32415026 (MEASURED): 200 genomes, 2,099,680 proteins, 16 CPUs, 18 HMMs. Job 4 min 21 s, script 209.3 s, MaxRSS 3.1 GB, 1,367 candidate rows (1,299 proteins), 0 bad id prefixes. The first job (32412179, `short`) waited about 3 h on the 384-CPU user limit and was cancelled; the timing job ran on `stajichlab`.
- Projection (PROPOSED, UNVERIFIED): 6 chunks of 4,000 genomes, about 64 to 75 min per job, about 47 GB (up to 62 GB), so submit with `--mem=96G`; about 159,773 candidate rows in total. The default `--mem=32G` of `slurm/search_chunk.sbatch` is too small for 4,000 genomes.
- Unusual counts: `tyrosinase` 439 rows, `hppd` 154, `omt` 125 (all within 60 bits of the floor), against 0 for `aspulvinone_e_synthase`.

Second commit of step 3 (2026-11-23): UniProt Q8NKC1 (*S. pombe* COMT2, reviewed, 281 residues) is a new `omt` decoy (`curation_scripts/melanin_myco_decoys/add_omt_decoy.py`). Step 2 for mycosporine was re-run (job 32415136). The highest `omt` decoy score is 120.7 (was 100.6). The `omt` floor follows the owner rule "floor = highest decoy score" through `--floor-override omt=120.7`. The timing job stays on the 100.6 table (kept in `results/step3/melanin_myco/`). In the timing chunk 84 of 125 `omt` rows lie between 100.6 and 120.7. Old step 2 tables: `results/step2/mycosporine.pre_Q8NKC1/`.
