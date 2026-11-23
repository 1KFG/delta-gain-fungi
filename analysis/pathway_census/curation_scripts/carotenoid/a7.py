from h import *
import urllib.parse, time
def rh(key, q, limit=15):
    url = "https://www.rhea-db.org/rhea/?query=" + urllib.parse.quote(q) + f"&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit={limit}"
    try:
        t = raw("rhea", key, url, as_json=False)
    except KbFetchError as x:
        t = f"ERROR {x}"
    print(f"--- {key}: {q}"); print(t)
for i in "30975 15585 31519 31515 55660 55664 34475 32219 32239 32887 17653 56772".split():
    rh("id_"+i, f"rhea:{i}", 3)
for q in ["torulene", "3,4-didehydrolycopene", "neurosporaxanthin", "echinenone", "canthaxanthin", "phoenicoxanthin", "adonirubin", "3-hydroxyechinenone", "torularhodin", "trisporate", "trisporin", "12'-apo-beta-carotenal", "13-apo-beta-carotenone", "4-ketotorulene", "all-trans-retinal beta-carotene", "gamma-carotene"]:
    rh("q_"+q.replace(" ","_").replace("'","p"), q)
# NCBI protein CAC16261
t = raw("ncbi_protein", "CAC16261_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=CAC16261&rettype=gp&retmode=text", as_json=False)
print(t[:2500])
