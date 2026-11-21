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
