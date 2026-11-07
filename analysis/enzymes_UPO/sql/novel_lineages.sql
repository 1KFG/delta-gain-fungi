-- Lineages (class/order/family) in BFD with no NCBI-annotated genome (status baseline_ok).
-- Restricted to genomes with BFD gene models (gene_info). Run from Fungi_BFD_runs/tables.
create view st as select * from read_csv('/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/prior_annotation_status.tsv', delim='\t', header=true, all_varchar=true);
create table lt as select distinct LOCUSTAG from 'gene_info.parquet';
create table t as select s.*, coalesce(st.status,'NO_MATCH') status from 'samples.parquet' s left join st on s.ASMID=st.folder where s.LOCUSTAG in (select LOCUSTAG from lt);
.mode markdown
select status, count(*) n from t group by 1 order by 2 desc;
select PHYLUM, CLASS, count(*) n_bfd_genomes, count(distinct "ORDER") n_orders, count(distinct GENUS) n_genera from t
 where CLASS in (select CLASS from t group by 1 having sum((status='baseline_ok')::int)=0) group by all order by 3 desc;
select PHYLUM, count(*) orders_no_ncbi, sum(n) genomes from (select PHYLUM, "ORDER", count(*) n from t group by all having sum((status='baseline_ok')::int)=0) group by 1 order by 3 desc;
select PHYLUM, count(*) families_no_ncbi, sum(n) genomes from (select PHYLUM, FAMILY, count(*) n from t group by all having sum((status='baseline_ok')::int)=0) group by 1 order by 3 desc;
select CLASS, status, count(*) n from t where CLASS in ('Lecanoromycetes','Lichinomycetes','Arthoniomycetes','Coniocybomycetes','Candelariomycetes') group by all order by 1,2;
select GENUS, SPECIES, status from t where GENUS ilike 'Antarct%' or GENUS ilike 'Cryomyces%' or GENUS ilike 'Friedmanniomyces%';
