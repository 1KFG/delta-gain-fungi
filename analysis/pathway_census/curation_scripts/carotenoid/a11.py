from h import *
import urllib.parse
def cs(term):
    url = "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/?term=" + urllib.parse.quote(term) + "&page=1&size=8"
    try:
        d = raw("chebi_search", term, url)
    except KbFetchError as x:
        print("ERR", term, str(x)[-40:]); return
    print(f"--- {term}", [(r.get('_source', r).get("chebi_accession"), r.get('_source', r).get("name"), r.get('_source', r).get("formula")) for r in d.get("results", [])[:8]])
for t in ["zzzqqqxx", "ketotorulene", "trisporin B", "trisporic acid B", "methyl trisporate B", "trisporol", "torulene", "4-keto-torulene", "trisporate"]:
    cs(t)
