-- BFD tables for antiSMASH lite results (schema v2, 2026-11-22). Run from this folder:
--   duckdb -c ".read build_bgc_parquet.sql"
-- Inputs: per-genome gz tables written by extract_antismash_subset.py (--minimal --topn 5), read from IN_DIR below
-- (in the BFD layout: results/function/antismash/*/*.antismash.<table>.tsv.gz; the `genome` column of every table
-- is the LOCUSTAG because the module passes meta.locustag as the genome id). Outputs: four zstd parquet tables in db_test/.
-- Join key between the tables: (locustag, record, region). gene_id joins to gene_info.gene_id (and LOCUSTAG).
-- No sequences are stored. Zero-cluster genomes give header-only files and no rows; an EMPTY (0-byte) file makes read_csv fail, so
-- skip files with size 0 when building the glob.
-- TEST MODE (these short ids are not LOCUSTAGs): the lt table derives the LOCUSTAG from the gene id prefix. In BFD replace the
-- lt statement with:  CREATE OR REPLACE TEMP TABLE lt AS SELECT DISTINCT genome, genome AS locustag FROM r;
CREATE OR REPLACE TEMP TABLE r  AS SELECT * FROM read_csv('IN_DIR/*.regions.tsv.gz',     delim='\t', header=true, auto_detect=true, types={'contig_edge':'BOOLEAN'});
CREATE OR REPLACE TEMP TABLE c  AS SELECT * FROM read_csv('IN_DIR/*.cds.tsv.gz',         delim='\t', header=true, auto_detect=true, types={'locus_tag':'VARCHAR','protein_id':'VARCHAR'});
CREATE OR REPLACE TEMP TABLE k  AS SELECT * FROM read_csv('IN_DIR/*.known.tsv.gz',       delim='\t', header=true, auto_detect=true);
CREATE OR REPLACE TEMP TABLE p  AS SELECT * FROM read_csv('IN_DIR/*.known_pairs.tsv.gz', delim='\t', header=true, auto_detect=true, types={'query_locus_tag':'VARCHAR','mibig_protein':'VARCHAR'});
CREATE OR REPLACE TEMP TABLE lt AS SELECT genome, any_value(split_part(locus_tag,'_',1)) AS locustag FROM c GROUP BY genome;
CREATE OR REPLACE TEMP TABLE best AS
  SELECT genome, record, region, arg_max(mibig_accession, (similarity, -rank)) AS best_mibig, max(similarity) AS best_similarity,
         arg_max(rank, (similarity, -rank)) AS best_rank
  FROM k WHERE rank <= 5 GROUP BY genome, record, region;

-- 1. bgc_cluster: one row per antiSMASH region. Sorted by locustag, record, region.
COPY (
  SELECT lt.locustag, r.record, r.region::USMALLINT AS region, r.start::UINTEGER AS start, r."end"::UINTEGER AS "end", r.length::UINTEGER AS length,
         r.products, r.categories, r.contig_edge, r.core_start::UINTEGER AS core_start, r.core_end::UINTEGER AS core_end,
         r.candidate_kinds, r.n_cds::USMALLINT AS n_cds, r.n_biosynthetic::USMALLINT AS n_biosynthetic, r.n_additional::USMALLINT AS n_additional,
         b.best_mibig, b.best_similarity::UTINYINT AS best_similarity, b.best_rank::UTINYINT AS best_rank
  FROM r JOIN lt USING (genome) LEFT JOIN best b ON b.genome=r.genome AND b.record=r.record AND b.region=r.region
  ORDER BY lt.locustag, r.record, r.region) TO 'db_test/bgc_cluster.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);

-- 2. bgc_cluster_gene: one row per gene inside a region. Sorted by locustag, record, region, start.
COPY (
  SELECT lt.locustag, c.record, c.region::USMALLINT AS region, c.locus_tag AS gene_id, c.start::UINTEGER AS start, c."end"::UINTEGER AS "end",
         (c.strand='-') AS is_minus, c.gene_kind, c.functions
  FROM c JOIN lt USING (genome) ORDER BY lt.locustag, c.record, c.region, c.start) TO 'db_test/bgc_cluster_gene.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);

-- 3. bgc_mibig_hit: one row per kept MIBiG hit (rank 1 to 5). Sorted by locustag, record, region, rank.
COPY (
  SELECT lt.locustag, k.record, k.region::USMALLINT AS region, k.rank::UTINYINT AS rank, k.mibig_accession, k.hit_description AS compounds, k.hit_cluster_type,
         k.similarity::UTINYINT AS similarity, k.core_gene_hits::USMALLINT AS core_gene_hits, k.gene_hits::USMALLINT AS gene_hits, k.n_entry_proteins::USMALLINT AS n_entry_proteins,
         k.synteny_score::SMALLINT AS synteny_score, k.blast_score::FLOAT AS blast_score
  FROM k JOIN lt USING (genome) WHERE k.rank <= 5 ORDER BY lt.locustag, k.record, k.region, k.rank) TO 'db_test/bgc_mibig_hit.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);

-- 4. bgc_mibig_pair: our gene -> MIBiG gene, one row per pair. Sorted by locustag, record, region, rank, gene_id.
-- Query coordinates and strand are not stored (join bgc_cluster_gene on gene_id); the MIBiG annotation text is not stored.
COPY (
  SELECT lt.locustag, p.record, p.region::USMALLINT AS region, p.rank::UTINYINT AS rank, p.mibig_accession, p.query_locus_tag AS gene_id,
         p.mibig_gene, p.mibig_protein, p.perc_ident::UTINYINT AS perc_ident, p.perc_coverage::FLOAT AS perc_coverage, p.blast_score::FLOAT AS blast_score, p.evalue::DOUBLE AS evalue
  FROM p JOIN lt USING (genome) WHERE p.rank <= 5 ORDER BY lt.locustag, p.record, p.region, p.rank, p.query_locus_tag) TO 'db_test/bgc_mibig_pair.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);
