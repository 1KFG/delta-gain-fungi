import json,csv,collections,re
N='/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/lib/'
D='/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/'
d=json.load(open(N+'ncbi_accessions.json'))['reports']
kept={r['ACCESSION'] for r in csv.DictReader(open(N+'ncbi_accessions.csv'))}
st={r['accession']:r['status'] for r in csv.DictReader(open(D+'prior_annotation_status.tsv'),delimiter='\t')}
def key(a):
    o=a.get('organism',{}); inf=o.get('infraspecific_names',{}) or {}
    s=inf.get('strain') or inf.get('isolate') or ''
    return (re.sub(r'\W+','',o.get('organism_name','').lower()), re.sub(r'\W+','',s.lower()))
grp=collections.defaultdict(list)
for a in d: grp[key(a)].append(a['accession'])
c=collections.Counter(); ex=[]
for a in d:
    acc=a['accession']
    if acc in kept: continue
    ai=a.get('annotation_info')
    if not ai: continue
    others=[x for x in grp[key(a)] if x in kept]
    if not others: c['annotated_dropped_no_kept_sibling']+=1; continue
    o=others[0]
    c['sibling_'+st.get(o,'?')]+=1
    if st.get(o,'').startswith('no_prior'): ex.append((acc,o))
print(dict(c)); print(ex[:5])
