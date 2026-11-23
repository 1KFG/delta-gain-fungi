from h import *
import urllib.parse, time
F = "accession,reviewed,gene_names,organism_name,length,ec,protein_name"
def s(key, q, size=25):
    url = "https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) + f"&fields={F}&format=tsv&size={size}"
    t = raw("uniprot_search", key, url, as_json=False)
    print(f"--- {key}: {q}"); print(t)
for pm in ["19584000", "21854466", "22445959", "16416328", "17822929"]:
    s("lit_"+pm, f"(lit_pubmed:{pm})")
e = uniprot_entry("A0A0D1E9J4"); print({k:v for k,v in e.items() if k!='sequence'})
# domain coordinates
for acc in "P37295 Q8X0Z1 Q9P854 Q9UUQ6 Q7Z859".split():
    d = raw("uniprot_full", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json")
    print("==", acc, d["sequence"]["length"])
    for f in d.get("features", []):
        if f["type"] in ("Region", "Domain", "Transmembrane", "Chain"):
            print("  UP", f["type"], f["location"]["start"]["value"], f["location"]["end"]["value"], f.get("description"))
    ip = raw("interpro", acc, f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}?page_size=100")
    for r in ip["results"]:
        m = r["metadata"]
        for p in r["proteins"]:
            for loc in p["entry_protein_locations"]:
                fr = [(x["start"], x["end"]) for x in loc["fragments"]]
                print("  IP", m["source_database"], m["accession"], m["name"], fr)
    time.sleep(0.5)
