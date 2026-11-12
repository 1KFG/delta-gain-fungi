-- Model sets for bin/check_codons.py: BFD flags vs genome sequence.
-- Run from Fungi_BFD_runs/tables with duckdb 1.1.3; ${IDS} and ${OUT} are substituted by the driver.
--   muc_upo         every PF01328 model from a Mucoromycota genome
--   asco_upo        500 random Ascomycota PF01328 models (control)
--   muc_genome      3000 random gene models from Mucoromycota genomes (is the genome-wide rate real?)
--   glomero_genome  1000 random gene models from Glomeromycota genomes
--   asco_genome     1000 random gene models from Ascomycota genomes (baseline)
.mode csv
.separator "\t"
select setseed(0.42);
create temp table ids as select protein_id from read_csv('${IDS}', delim='\t', header=true);
create temp table tx as
  select t.transcript_id, t.has_start_codon, t.has_stop_codon, s.ASMID, s.PHYLUM,
         t.transcript_id in (select protein_id from ids) as is_upo
  from read_parquet('gene_transcripts.parquet') t
  join read_parquet('species.parquet') s on s.LOCUSTAG = split_part(t.transcript_id, '_', 1)
  where s.PHYLUM in ('Mucoromycota', 'Glomeromycota', 'Ascomycota');
create temp table sel as
  select 'muc_upo' as "set", * from tx where is_upo and PHYLUM = 'Mucoromycota'
  union all select 'asco_upo', * from (select * from tx where is_upo and PHYLUM = 'Ascomycota' order by random() limit 500)
  union all select 'muc_genome', * from (select * from tx where PHYLUM = 'Mucoromycota' order by random() limit 3000)
  union all select 'glomero_genome', * from (select * from tx where PHYLUM = 'Glomeromycota' order by random() limit 1000)
  union all select 'asco_genome', * from (select * from tx where PHYLUM = 'Ascomycota' order by random() limit 1000);
copy (select "set", transcript_id, ASMID, has_start_codon, has_stop_codon from sel)
  to '${OUT}/codon_models.tsv' (header, delimiter '\t');
copy (select c.transcript_id, c.chrom, c.start, c."end", c.strand
      from read_parquet('gene_CDS.parquet') c where c.transcript_id in (select transcript_id from sel))
  to '${OUT}/codon_cds.tsv' (header, delimiter '\t');
