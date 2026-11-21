# pathway_census

Shared pathway census package (carotenoid pilot). Task 1: definition loader (`pathway_census.defs`).

Run tests with the system Python 3.12 (default `python3` is 3.9):

    /usr/bin/python3.12 -m pytest analysis/pathway_census/tests -v

`pytest` and `pyyaml` were already present under `/usr/bin/python3.12`, so no venv was needed.
If they go missing, create `analysis/pathway_census/.venv` (git-ignored) and install `pyyaml pytest` there.

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

- `pathways/carotenoid.yaml`: 8 steps, pathway `verified: false`. User review of the step list is required before any census run (spec 5.8).
- `refs/carotenoid/<step>.faa`: reference proteins from UniProt (retrieved 2026-11-20), one file per step. Fused PSY/cyclase proteins appear in both step files.
- `refs/carotenoid/provenance.tsv`: one row per (step, ref) with source and verified flag.
- `refs/carotenoid/UNVERIFIED.md`: output of `bin/list_unverified_refs.py pathways/carotenoid.yaml` plus a hand-written list of open decisions.

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
- `domain_ok` needs `expected_domains` (step_id -> hmm names). It matches `hmm_name` in the `domains` frame, so the domain search must be run with the matching HMM library (for example Pfam ids from `pathways/carotenoid.yaml`).

Task 7 must filter with `quality.clean_mask(flags)` and never `fillna(False)`: `<NA>` means unknown, not clean. Task 7 must also report `quality.unknown_counts(flags)` in its output. Contamination is `<NA>` when fewer than `min_neighbours_with_hits` neighbours have a best hit.

## Copy matrix (Task 7)

Code: `pathway_census/matrix.py`, `bin/build_matrix.py`. Columns and rules are in the `matrix.py` docstring.
- Full genome x step cross product, zeros included. `has_proteins=False` genomes get zero rows and are excluded from `n_genomes`.
- `copy_raw` = assigned rows. Mutually exclusive buckets, priority flagged > unknown > filtered: `copy_flagged` = any required flag known True; `copy_unknown` = no flag True and some required flag NA; `copy_filtered` = all required flags known False (default require: truncated, contam_suspect). `copy_filtered + copy_unknown + copy_flagged = copy_raw` per row. `hap_dup_suspect` is counted in `n_hap_suspect` / `n_hap_unknown` only.
- `phylum_summary` uses `copy_filtered` by default (`count_col="copy_raw"` for the unfiltered table). `n_unknown` = genomes with unevaluated copies.
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
- Chunk outputs are named `<pathway>_absence.chunk<ID>[.PROVISIONAL].tsv` and `<pathway>_unevaluated_genomes.chunk<ID>...tsv`; ID = `SLURM_ARRAY_TASK_ID`, else the `PC_LOCUSTAGS` basename minus its last suffix. Merge with: header from one file, then `tail -q -n +2` of the rest (skip the leading `#` lines first), (output name PROPOSED; Chapter A reads it) e.g. `for f in *_absence.chunk*.tsv; do grep -v '^#' "$f"; done | awk 'NR==1||$1!="locustag"' > results/carotenoid_absence_summary.tsv`. `unevaluated_genomes` is the same list in every chunk (it is not chunked): keep one copy.
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
- Input: `paired.deposited_protein_path(deposited_root, folder)` for each pair genome (deposited root `/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM`; `folder` from the genome table). A table `protein_id -> locustag` (the caller supplies it; deposited protein IDs do NOT use the 8-hex BFD prefix).
- Output: a matrix in `build_matrix` format (same pathway steps, same HMMs and assignment rules as the BFD arm), keyed by the BFD `locustag`, with `has_proteins` set from the deposited-file flag (`has_deposited_proteins` and the file exists).
- Design choice, to decide in the real-data pass (two options): (a) rewrite deposited protein IDs to `<locustag>_<deposited id>` before assign / quality / matrix; or (b) add a `genome_of` mapping parameter (protein_id -> locustag) to `matrix.build_matrix` and `quality.quality_flags`. Option (a) needs no code change but alters protein IDs in every output; option (b) keeps IDs and changes two functions. Neither is chosen.
- Change surface (places that assume the 8-hex `LOCUSTAG_` prefix; found by grep on 2026-11-20):
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
