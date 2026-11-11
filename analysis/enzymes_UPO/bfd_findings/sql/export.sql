-- Gene-model and taxonomy tables for bfd_findings/bin/report_tables.py.
-- Run from Fungi_BFD_runs/tables with duckdb 1.1.3; ${IDS} and ${OUT} are substituted by the driver.
-- LOCUSTAG is the protein/transcript ID prefix before the first '_'.
.mode csv
.separator "\t"
create temp table ids as select protein_id from read_csv('${IDS}', delim='\t', header=true);
create temp table tax as
  select s.LOCUSTAG, s.ASMID, s.PHYLUM, s.SUBPHYLUM, s.CLASS, s."ORDER" as ORDER_, s.FAMILY,
         s.GENUS, s.SPECIES, a.contig_count, a.N50_bp, a.total_length_bp
  from read_parquet('species.parquet') s left join read_parquet('asm_stats.parquet') a using (ASMID);

-- 1. Every PF01328 gene model: coordinates, BFD completeness flags, exon count, taxonomy.
copy (
  select t.transcript_id as protein_id, t.gene_id, t.chrom, t.start, t."end", t.strand,
         t.is_partial, t.has_start_codon, t.has_stop_codon, t.CDS_length,
         e.n_exons, x.*
  from read_parquet('gene_transcripts.parquet') t
  join ids on t.transcript_id = ids.protein_id
  left join (select transcript_id, count(*) as n_exons from read_parquet('gene_exons.parquet')
             where transcript_id in (select protein_id from ids) group by 1) e
         on e.transcript_id = t.transcript_id
  left join tax x on x.LOCUSTAG = split_part(t.transcript_id, '_', 1)
) to '${OUT}/upo_models.tsv' (header, delimiter '\t');

-- 2. Genome-wide gene-model completeness, every BFD genome (the baseline for the UPO rates).
copy (
  select split_part(transcript_id, '_', 1) as LOCUSTAG, count(*) as n_models,
         count(*) filter (where has_start_codon = 'TRUE' and has_stop_codon = 'TRUE') as n_complete,
         count(*) filter (where has_start_codon = 'FALSE' and has_stop_codon = 'TRUE') as n_no_start,
         count(*) filter (where has_start_codon = 'TRUE' and has_stop_codon = 'FALSE') as n_no_stop,
         count(*) filter (where has_start_codon = 'FALSE' and has_stop_codon = 'FALSE') as n_neither,
         count(*) filter (where has_start_codon is null or has_stop_codon is null) as n_flags_null,
         count(*) filter (where CDS_length % 3 <> 0) as n_cds_not_mult3
  from read_parquet('gene_transcripts.parquet') group by 1
) to '${OUT}/genome_qc.tsv' (header, delimiter '\t');

-- 3. Taxonomy of every genome (denominators: genomes per class).
copy (select * from tax) to '${OUT}/genome_tax.tsv' (header, delimiter '\t');
