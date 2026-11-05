import csv,os
N='/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/'
D='/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/'
rows=list(csv.DictReader(open(N+'lib/ncbi_accessions.csv')))
print('csv rows',len(rows), 'mtime',os.path.getmtime(N+'lib/ncbi_accessions.csv'))
byacc={r['ACCESSION']:r for r in rows}
st=list(csv.DictReader(open(D+'prior_annotation_status.tsv'),delimiter='\t'))
import collections
c=collections.Counter()
nodl=[];nop=[]
for r in st:
    s=r['status']
    if s=='annotated_NOT_downloaded': c[('notdl', r['accession'] in byacc)]+=1; 
    if s=='annotated_NOT_downloaded' and r['accession'] in byacc: nodl.append(byacc[r['accession']])
    if s=='annotated_folder_missing_protein': c[('noprot', r['accession'] in byacc)]+=1; nop.append(r)
print(dict(c))
# write sub-CSV for the not-downloaded ones that are in the csv
w=csv.DictWriter(open(D+'gap_notdownloaded.csv','w',newline=''),fieldnames=rows[0].keys()); w.writeheader(); w.writerows(nodl)
# folders for noprot: list accession,folder
with open(D+'gap_noprotein_folders.tsv','w') as f:
    for r in nop: f.write(r['accession']+'\t'+r['folder']+'\n')
print(len(nodl),len(nop)); print(nodl[:2]); print(nop[:2])
