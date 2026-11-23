from h import *
import time
accs = """P24322 Q9P885 P37295 Q8X0Z1 Q9P854 Q9UUQ6 Q7Z859 P21334 Q8X0Z0 P54982 O13506 S0ENQ2 A1KQY3 A7UXI1 Q1K615 F6IBC7 B5LEC8 Q3HR17 B6D970 G8Z9M3 A0ABR3AZC2 A5HXP5 B9DR41 Q01213
P08524 O14230 A0A1D8PH78 P18900 P29704 P36596 Q92459 Q4WAG4 P50264 P46882 Q10062 P23623 B8NI10 Q7S860 A0A0D1E6L2 A0A384JIP2 Q4W9H1 G2QIL8 Q52008 P54115 P46367 P08157 Q04458 P38067 P10614 P54781 P23295 P21595 Q9Y8G7 Q12181 P39692 O94613 P21158 P38715 Q12458 P14065 A0A0D1E9J4""".split()
for a in accs:
    e = uniprot_entry(a)
    d = raw("uniprot_full", a, f"https://rest.uniprot.org/uniprotkb/{a}.json")
    pf = sorted(x["id"] for x in d.get("uniProtKBCrossReferences", []) if x["database"] in ("Pfam", "NCBIfam"))
    pdn = d["proteinDescription"]
    nm = (pdn.get("recommendedName") or (pdn.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value")
    au = d["entryAudit"]
    print(a, "|", e["organism"], "|", e["gene_names"][:3], e["length"], e["ec"], "rev" if e["reviewed"] else "unrev", "|", nm, "|", pf, f"| v{au['entryVersion']}/s{au['sequenceVersion']}", "| pm", e["pmids"][:6])
kg = raw("kegg", "C08607", "https://rest.kegg.jp/get/C08607", as_json=False); print(kg[:200])
