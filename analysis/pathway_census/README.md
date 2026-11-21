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
- `copy_raw` = assigned rows. `copy_filtered` = assigned rows known False for every `require` flag (default truncated, contam_suspect). `copy_unknown` = assigned rows with a required flag NA. `hap_dup_suspect` is counted in `n_hap_suspect` / `n_hap_unknown` only.
- `phylum_summary` uses `copy_filtered` by default (`count_col="copy_raw"` for the unfiltered table). `n_unknown` = genomes with unevaluated copies.
- The real-data run (brief step 5) is deferred. Do not read PROVISIONAL outputs as results.
