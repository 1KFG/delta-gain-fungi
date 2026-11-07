create view st as select * from read_csv('/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/prior_annotation_status.tsv', delim='\t', header=true, all_varchar=true);
create view sm as select s.*, coalesce(st.status,'NO_MATCH') status, (g.LOCUSTAG is not null) has_pfam from 'samples.parquet' s left join st on s.ASMID=st.folder left join '${SCRATCHPAD}/pfam_genomes.parquet' g on s.LOCUSTAG=g.LOCUSTAG;
create view h as select split_part(protein_id,'_',1) LOCUSTAG, protein_id, min(domain_i_evalue) ev, max(full_seq_score) sc from '${SCRATCHPAD}/pf01328.parquet' group by 1,2;
create view per as select sm.LOCUSTAG, sm.PHYLUM, sm.status, sm.has_pfam, count(h.protein_id) n_all, count(case when h.ev<1e-5 then 1 end) n_e5 from sm left join h using(LOCUSTAG) group by all;
.mode markdown
select PHYLUM, case when status='baseline_ok' then 'NCBI' when status like 'no_prior%' then 'BFD_only' else status end ann,
 count(*) genomes, sum(has_pfam::int) w_pfam, sum(n_all) pf01328_prot, sum(n_e5) prot_e5, sum((n_e5>0)::int) genomes_w_hit_e5,
 median(n_e5) filter (where has_pfam) med_per_genome, max(n_e5) mx
from per group by all order by PHYLUM, ann;
select min(ev), max(ev), max(sc), min(sc) from h;
