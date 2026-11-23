from h import *
import re
gp = raw("ncbi_protein", "ADU04395_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ADU04395&rettype=gp&retmode=text", as_json=False)
print(gp.split("FEATURES")[0])
seq = re.sub(r"[^a-z]", "", gp.split("ORIGIN")[1].split("//")[0]).upper()
u = uniprot_entry("G8Z9M3")
print("len", len(seq), u["length"], "identical:", seq == u["sequence"])
d = raw("uniprot_full", "G8Z9M3", "https://rest.uniprot.org/uniprotkb/G8Z9M3.json")
print([ (x["database"], x["id"]) for x in d["uniProtKBCrossReferences"] if x["database"] in ("EMBL","RefSeq")])
print([r.get("citation",{}).get("title") for r in d["references"]])
d = raw("uniprot_full", "A0ABR3AZC2", "https://rest.uniprot.org/uniprotkb/A0ABR3AZC2.json")
print([ (x["database"], x["id"], [p["value"] for p in x.get("properties",[])]) for x in d["uniProtKBCrossReferences"] if x["database"] in ("EMBL","RefSeq")])
print([(r.get("citation",{}).get("title"), r.get("citation",{}).get("publicationDate")) for r in d["references"]])
for p in ["22445959","25504221","8828234","18931040","16928467"]:
    print(p, pubmed_abstract(p)["title"][:70])
