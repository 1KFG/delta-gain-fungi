-- Gene-model QC flags and assembly IDs for the PF01328 proteins, for candidates/bin/shortlist.py.
-- Run from Fungi_BFD_runs/tables with duckdb 1.1.3; ${IDS} and ${OUT} are substituted by the driver.
.mode csv
.separator "\t"
create temp table ids as select protein_id from read_csv('${IDS}', delim='\t', header=true);
copy (select t.transcript_id, t.is_partial, t.has_start_codon, t.has_stop_codon, t.CDS_length
      from read_parquet('gene_transcripts.parquet') t join ids on t.transcript_id = ids.protein_id)
  to '${OUT}/geneqc.tsv' (header, delimiter '\t');
copy (select LOCUSTAG, ASMID from read_parquet('samples.parquet'))
  to '${OUT}/asm.tsv' (header, delimiter '\t');
