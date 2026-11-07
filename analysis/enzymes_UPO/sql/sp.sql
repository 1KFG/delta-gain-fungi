create view st as select * from read_csv('/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/prior_annotation_status.tsv', delim='\t', header=true, all_varchar=true);
create table j as select s.*, case when st.status='baseline_ok' then 'NCBI' when st.status like 'no_prior%' then 'BFD_only' else coalesce(st.status,'NO_MATCH') end ann from 'samples.parquet' s left join st on s.ASMID=st.folder;
create table b as select j.*, bg.complete_pct, bg.fragmented_pct, a.total_length_bp, a.N50_bp from j left join 'busco_genome.parquet' bg using(ASMID) left join 'asm_stats.parquet' a using(ASMID) where PHYLUM in ('Mucoromycota','Glomeromycota','Chytridiomycota','Zoopagomycota','Basidiobolomycota','Mortierellomycota','Kickxellomycota','Blastocladiomycota','Neocallimastigomycota','Entomophthoromycota','Monoblepharomycota');
.mode markdown
select PHYLUM, ann, count(*) genomes, count(distinct SPECIES) species,
 count(distinct SPECIES) filter (where ann='BFD_only' and SPECIES not in (select SPECIES from j where ann='NCBI')) new_species_vs_ncbi,
 round(median(complete_pct),1) busco_c_med, sum((complete_pct<80)::int) busco_lt80, round(median(total_length_bp)/1e6,1) med_Mb, round(median(N50_bp)/1e3,1) med_N50_kb
from b where ann in ('NCBI','BFD_only') group by all order by 1,2;
