create view st as select * from read_csv('/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/prior_annotation_status.tsv', delim='\t', header=true, all_varchar=true);
create table h as select split_part(protein_id,'_',1) LOCUSTAG, protein_id, count(*) ndom, min(domain_i_evalue) ev, max(full_seq_score) sc, max(pfam_len) plen, min(hmm_from) hf, max(hmm_to) ht from '${SCRATCHPAD}/pf01328.parquet' group by 1,2;
create table hx as select h.*, s.PHYLUM, case when st.status='baseline_ok' then 'NCBI' when st.status like 'no_prior%' then 'BFD_only' else coalesce(st.status,'NO_MATCH') end ann
 from h join 'samples.parquet' s using(LOCUSTAG) left join st on s.ASMID=st.folder
 where s.PHYLUM in ('Mucoromycota','Glomeromycota','Chytridiomycota','Zoopagomycota','Basidiobolomycota','Monoblepharomycota','Aphelidiomycota','Calcarisporiellomycota');
create table hy as select hx.*, gp.length, (sp.protein_id is not null) has_sp from hx left join 'gene_proteins.parquet' gp using(protein_id) left join (select distinct protein_id from 'signalp.signal_peptide.parquet') sp using(protein_id);
copy hy to '${SCRATCHPAD}/early_pf01328_hits.parquet';
.mode markdown
select PHYLUM, ann, count(*) n, median(length) med_len, min(length) mn, max(length) mx, sum((length<250)::int) lt250, sum((length between 250 and 350)::int) l250_350, sum((length>350)::int) gt350, round(avg(has_sp::int)*100,1) pct_sp, round(avg(((ht-hf+1)>=0.8*plen)::int)*100,1) pct_hmm80, any_value(plen) plen from hy group by all order by 1,2;
