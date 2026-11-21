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
- Empty GENUS: 349 rows. 49 parsed from SPECIES (`parsed_species`; includes 8 Antarctolithica). 300 left with null genus (SPECIES is a higher rank, e.g. "Ascomycota sp.").
- `has_proteins` = TSV `protein_file_ok == 1` (prior-annotation protein file). It does not mean BFD lacks proteins.
