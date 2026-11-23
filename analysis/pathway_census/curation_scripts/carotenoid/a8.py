from h import *
import urllib.parse, re
t = raw("ncbi_protein", "CAC16261_gp", "", as_json=False) if False else None
gp = raw("ncbi_protein", "CAC16261_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=CAC16261&rettype=gp&retmode=text", as_json=False)
seq = re.sub(r"[^a-z]", "", gp.split("ORIGIN")[1].split("//")[0]).upper()
q = uniprot_entry("Q3HR17")["sequence"]; b = uniprot_entry("B5LEC8")["sequence"]
print("CAC16261 len", len(seq), "== Q3HR17:", seq == q, "; diffs vs B5LEC8:", [(i+1, x, y) for i,(x,y) in enumerate(zip(seq,b)) if x!=y], len(b))
def cs(term):
    url = "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/?term=" + urllib.parse.quote(term) + "&page=1&size=8"
    try:
        d = raw("chebi_search", term, url)
    except KbFetchError as x:
        print("ERR", term, x); return
    print(f"--- {term}")
    for r in d.get("results", [])[:8]:
        s = r.get("_source", r)
        print("  ", s.get("chebi_accession"), "|", s.get("name"), "|", s.get("formula"), "|", s.get("charge"), "| stars", s.get("stars"))
for term in ["neurosporaxanthin", "phoenicoxanthin", "adonirubin", "3-hydroxyechinenone", "torularhodin", "4-ketotorulene", "3-hydroxy-4-ketotorulene", "12'-apo-beta-carotenal", "13-apo-beta-carotenone", "trisporin", "dihydrotrisporin", "methyl trisporate", "trisporic acid", "torulene", "geranylgeranyl diphosphate", "farnesyl diphosphate", "isopentenyl diphosphate", "4'-apo-beta-carotenal", "3,4-didehydrolycopene"]:
    cs(term)
