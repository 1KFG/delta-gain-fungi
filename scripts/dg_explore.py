import duckdb
b='/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/results/protein_novelty_n2000/'
c=duckdb.connect(); c.execute("PRAGMA threads=4")
c.execute(f"""create table g as select regexp_extract(query_id,'^(.*)__',1) asmid,
 overall_status s, nr_cluster_status ns, try_cast(nr_cluster_evalue as double) ev, try_cast(nr_cluster_scov as double) scov
 from read_csv('{b}bfd_novelty_bins.tsv',delim='\t',header=true,quote='',all_varchar=true)""")
c.execute(f"create table cl as select * from read_csv('/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/results/genome_classify/genome_classification/genome_classification.tsv',delim='\t',header=true,all_varchar=true)")
c.execute("""create table pg as select asmid, count(*) n, avg((s='no_hit')::int) nohit, avg((s='weak_hit')::int) weak from g group by 1""")
def show(t,q):
    print('==',t); print(c.execute(q).fetchdf().to_string(index=False))
show('genomes in run / per-genome nohit quantiles',"select count(*) genomes, sum(n) reps, median(n) med_reps, quantile_cont(nohit,[.1,.25,.5,.75,.9,.99]) q from pg")
show('share of all no-hit reps held by top 5/10/25% genomes',"""with t as (select asmid, nohit*n k, sum(nohit*n) over () tot, row_number() over (order by nohit*n desc) r, count(*) over () N from pg)
select 'top5%' g, sum(k)/max(tot) frac from t where r<=0.05*N union all select 'top10%',sum(k)/max(tot) from t where r<=0.10*N union all select 'top25%',sum(k)/max(tot) from t where r<=0.25*N""")
show('by cluster_class',"""select cluster_class, count(*) genomes, round(avg(nohit),4) mean_nohit, round(median(nohit),4) med_nohit, round(avg(weak),3) mean_weak from pg join cl using(asmid) group by 1 order by 2 desc""")
show('by accession type',"select substr(asmid,1,3) acc, count(*) genomes, round(avg(nohit),4) mean_nohit, round(median(nohit),4) med from pg group by 1")
show('by genus (n>=15), top mean nohit',"""select split_part(existing_species,' ',1) genus, count(*) genomes, round(avg(nohit),4) mean_nohit from pg join cl using(asmid) group by 1 having count(*)>=15 order by 3 desc limit 10""")
show('by genus (n>=15), lowest mean nohit',"""select split_part(existing_species,' ',1) genus, count(*) genomes, round(avg(nohit),4) mean_nohit from pg join cl using(asmid) group by 1 having count(*)>=15 order by 3 limit 5""")
show('nohit vs reps/genome correlation',"select corr(nohit,n) r_nreps, corr(nohit,ln(n)) r_lnreps from pg")
show('weak-hit evalue/scov among weak (nr)',"select quantile_cont(ev,[.1,.5,.9]) ev_q, quantile_cont(scov,[.1,.5,.9]) scov_q from g where ns='weak_hit' and ev is not null")
