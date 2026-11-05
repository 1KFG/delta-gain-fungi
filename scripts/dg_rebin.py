"""Re-examine the 2000-genome DIAMOND-vs-nr hits.
 (1) reproduce the old binning (lowest-E-value row) and check it against bfd_novelty_bins.tsv
 (2) best-class re-bin: how much does 'weak' shrink?
 (3) self-hit evidence: top-hit identity for genomes WITH prior NCBI annotation (their proteins are in nr)
     vs genomes WITHOUT it.
"""
import duckdb, sys
F='/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/'
R=F+'results/protein_novelty_n2000/'
c=duckdb.connect(); c.execute("PRAGMA threads=16"); c.execute("PRAGMA memory_limit='100GB'"); c.execute(f"PRAGMA temp_directory='/scratch/jstajich/duck_rebin'")
cols="{'q':'VARCHAR','s':'VARCHAR','pident':'DOUBLE','len':'BIGINT','mm':'BIGINT','go':'BIGINT','qs':'BIGINT','qe':'BIGINT','ss':'BIGINT','se':'BIGINT','ev':'DOUBLE','bits':'DOUBLE','qc':'DOUBLE','sc':'DOUBLE','title':'VARCHAR'}"
c.execute(f"""create table agg as select q, count(*) nrows,
  arg_min(case when ev<=1e-5 and qc>=70 and sc>=70 then 2 when ev<=1e-3 then 1 else 0 end, struct_pack(a:=ev,b:=-bits)) old_cls,
  max(case when ev<=1e-5 and qc>=70 and sc>=70 then 2 when ev<=1e-3 then 1 else 0 end) best_cls,
  arg_max(pident, bits) top_pident, arg_max(qc, bits) top_qc, arg_max(sc, bits) top_sc, min(ev) min_ev
  from read_csv('{R}raw_hits/nr_cluster/*.zst', delim='\t', header=false, quote='', escape='', columns={cols}, compression='zstd')
  group by q""")
print('queries with >=1 hit:', c.execute("select count(*) from agg").fetchone()[0])
c.execute(f"""create table u as select query_id q, regexp_extract(query_id,'^(.*)__',1) asmid, nr_cluster_status old_status
  from read_csv('{R}bfd_novelty_bins.tsv', delim='\t', header=true, quote='', all_varchar=true)""")
c.execute(f"""create table st as select * from read_csv('{F}prior_annotation_status.tsv', delim='\t', header=true, all_varchar=true)""")
c.execute("""create table gc as select u.asmid, coalesce(case when s.status='baseline_ok' then 'baseline_annotated'
   when s.status like 'annotated%' then 'annotated_pending' when s.status like 'no_prior%' then 'new_unannotated' end,'not_in_ncbi_table') cls
   from (select distinct asmid from u) u
   left join st s on s.folder=u.asmid or s.accession=regexp_extract(u.asmid,'^(GC[AF]_[0-9]+\\.[0-9]+)',1)""")
c.execute("""create table m as select u.q, u.asmid, g.cls gclass, u.old_status, coalesce(a.old_cls,-1) old_cls, coalesce(a.best_cls,-1) best_cls, a.top_pident, a.top_qc, a.top_sc
   from u left join agg a using(q) left join gc g using(asmid)""")
def show(t,q):
    print('\n==',t); print(c.execute(q).fetchdf().to_string(index=False))
show('(1) reproduction: bins-file status vs recomputed old-rule class',
 "select old_status, case old_cls when 2 then 'strong_hit' when 1 then 'weak_hit' else 'no_hit' end recomputed, count(*) n from m group by 1,2 order by 1,2")
show('(2) best-class re-bin vs old rule (queries)',
 "select case old_cls when 2 then 'strong' when 1 then 'weak' else 'none' end old_rule, case best_cls when 2 then 'strong' when 1 then 'weak' else 'none' end best_class, count(*) n from m group by 1,2 order by 1,2")
show('(2b) totals old vs best-class',
 "select 'old' rule, sum((old_cls=2)::int) strong, sum((old_cls=1)::int) weak, sum((old_cls<1)::int) none from m union all select 'best-class', sum((best_cls=2)::int), sum((best_cls=1)::int), sum((best_cls<1)::int) from m")
show('(3) genome class sizes', "select gclass, count(distinct asmid) genomes, count(*) queries from m group by 1 order by 2 desc")
show('(3b) by genome class: old-rule fractions and top-hit identity',
 """select gclass, count(*) queries,
   round(avg((old_cls=2)::int),4) strong, round(avg((old_cls=1)::int),4) weak, round(avg((old_cls<1)::int),4) no_hit,
   round(avg((top_pident>=99)::int),4) top_pident_ge99, round(avg((top_pident>=95)::int),4) top_pident_ge95,
   round(median(top_pident),1) median_pident from m group by 1 order by 2 desc""")
show('(3c) per-genome: share of queries with top hit >=99% identity, quantiles by class',
 """select gclass, count(*) genomes, quantile_cont(f,[.1,.5,.9]) q10_50_90 from
   (select gclass, asmid, avg((top_pident>=99)::int) f from m group by 1,2) group by 1""")
show('(3d) identity bins of top hit, new vs baseline genomes (fraction of all queries in class)',
 """select gclass, case when top_pident is null then 'no hit' when top_pident>=99 then '99-100' when top_pident>=95 then '95-99' when top_pident>=70 then '70-95' when top_pident>=50 then '50-70' else '<50' end bin,
   count(*) n, round(count(*)*1.0/sum(count(*)) over (partition by gclass),4) frac from m where gclass in ('baseline_annotated','new_unannotated') group by 1,2 order by 1,2""")
c.execute(f"copy (select gclass, asmid, count(*) queries, avg((old_cls<1)::int) nohit_old, avg((old_cls=1)::int) weak_old, avg((best_cls=1)::int) weak_best, avg((top_pident>=99)::int) frac_ge99 from m group by 1,2) to '{F}reports/rebin_selfhit_n2000_per_genome.tsv' (delimiter '\t', header)")
