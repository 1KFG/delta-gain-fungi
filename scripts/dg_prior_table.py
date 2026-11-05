import json,os,csv,collections
L='/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/lib/ncbi_accessions.json'
A='/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/'
out='/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/'
txt=open(L).read()
try:
    d=json.loads(txt); recs=d['reports'] if isinstance(d,dict) else d
except Exception:
    recs=[json.loads(l) for l in txt.splitlines() if l.strip()]
folders={}
for h in os.listdir(A):
    p=h.split('_'); folders[p[0]+'_'+p[1]]=h
rows=[];c=collections.Counter()
for r in recs:
    a=r['accession']; ai=r.get('annotation_info') or {}
    pc=(ai.get('stats') or {}).get('gene_counts',{}).get('protein_coding')
    f=folders.get(a); hasp=False
    if f:
        pf=A+f+'/'+f+'_protein.faa.gz'
        hasp=os.path.exists(pf) and os.path.getsize(pf)>1024
    ann=bool(ai) and (pc is None or pc>0)
    status=('baseline_ok' if ann and hasp else 'annotated_NOT_downloaded' if ann and not f else 'annotated_folder_missing_protein' if ann else 'no_prior_annotation_present' if f else 'no_prior_annotation_not_downloaded')
    c[status]+=1
    rows.append([a,f or '',ai.get('provider',''),ai.get('release_date',''),pc if pc is not None else '',int(hasp),status,(r.get('assembly_info') or {}).get('release_date','') ,(r.get('source_database') or '')])
with open(out+'prior_annotation_status.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(['accession','folder','annotation_provider','annotation_release','protein_coding','protein_file_ok','status','assembly_release','source_database']); w.writerows(rows)
with open(out+'prior_annotation_gap_accessions.txt','w') as fh:
    for r in rows:
        if r[6] in('annotated_NOT_downloaded','annotated_folder_missing_protein'): fh.write(r[0]+'\n')
print(dict(c)); print('total',len(rows))
# present folders not in the JSON
js={r['accession'] for r in recs}
print('folders not in json:',sum(1 for k in folders if k not in js))
