-- Proposed BFD DuckDB tables for antiSMASH lite results. Run: duckdb -c ".read build_bgc_parquet.sql"
-- from analysis/bgc_subset/ with IN_DIR and OUT_DIR below. Reads the per-genome gz tables written by
-- extract_antismash_subset.py (also the BFD module outputs <LOCUSTAG>.antismash.{regions,cds,known,known_pairs}.tsv.gz).
-- No sequences are stored. Key = LOCUSTAG (same as results/function/*), joins to gene_info.gene_id.
CREATE OR REPLACE TEMP TABLE r AS SELECT * FROM read_csv('results/*.regions.tsv.gz', delim='\t', header=true, auto_detect=true);
CREATE OR REPLACE TEMP TABLE c AS SELECT * FROM read_csv('results/*.cds.tsv.gz', delim='\t', header=true, auto_detect=true, types={'locus_tag':'VARCHAR','protein_id':'VARCHAR'});
CREATE OR REPLACE TEMP TABLE k AS SELECT * FROM read_csv('results/*.known.tsv.gz', delim='\t', header=true, auto_detect=true);
-- genome label -> LOCUSTAG (prefix of the gene id). Test only: the BFD module already names files by LOCUSTAG.
CREATE OR REPLACE TEMP TABLE lt AS SELECT genome, any_value(split_part(locus_tag,'_',1)) AS locustag FROM c GROUP BY genome;
CREATE OR REPLACE TEMP TABLE best AS
  SELECT genome, record, region, arg_max(mibig_accession, (similarity, -rank)) AS best_mibig, max(similarity) AS best_similarity,
         arg_max(rank, (similarity, -rank)) AS best_rank
  FROM k GROUP BY genome, record, region;
COPY (
  SELECT lt.locustag, r.record, r.region::USMALLINT AS region, r.start::UINTEGER AS start, r."end"::UINTEGER AS "end", r.length::UINTEGER AS length,
         r.products, r.categories, (r.contig_edge)::BOOLEAN AS contig_edge, r.core_start::UINTEGER AS core_start, r.core_end::UINTEGER AS core_end,
         r.candidate_kinds, r.n_cds::USMALLINT AS n_cds, r.n_biosynthetic::USMALLINT AS n_biosynthetic, r.n_additional::USMALLINT AS n_additional,
         b.best_mibig, b.best_similarity::UTINYINT AS best_similarity, b.best_rank::UTINYINT AS best_rank
  FROM r JOIN lt USING (genome) LEFT JOIN best b ON b.genome=r.genome AND b.record=r.record AND b.region=r.region
  ORDER BY lt.locustag, r.record, r.region) TO 'db_test/bgc_cluster.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);
COPY (
  SELECT lt.locustag, c.record, c.region::USMALLINT AS region, c.locus_tag AS gene_id, c.start::UINTEGER AS start, c."end"::UINTEGER AS "end",
         (c.strand='-')::BOOLEAN AS is_minus, c.gene_kind, c.functions
  FROM c JOIN lt USING (genome) ORDER BY lt.locustag, c.record, c.region, c.start) TO 'db_test/bgc_cluster_gene.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);
COPY (
  SELECT lt.locustag, k.record, k.region::USMALLINT AS region, k.rank::UTINYINT AS rank, k.mibig_accession, k.hit_description AS compounds, k.hit_cluster_type,
         k.similarity::UTINYINT AS similarity, k.core_gene_hits::USMALLINT AS core_gene_hits, k.gene_hits::USMALLINT AS gene_hits, k.n_entry_proteins::USMALLINT AS n_entry_proteins,
         k.synteny_score::SMALLINT AS synteny_score, k.blast_score::FLOAT AS blast_score
  FROM k JOIN lt USING (genome) ORDER BY lt.locustag, k.record, k.region, k.rank) TO 'db_test/bgc_mibig_hit.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);
