from h import *
import json, time
from pathway_census.kb import formula_mass, net_charge
ids = "16584 17407 48861 27787 15948 62474 27740 17579 9638 17898 53157 15825 7540 4746 3362 80215 80214 80216 40968 53160 53175 211238".split()
out = {}
for i in ids:
    try:
        c = chebi_compound(i)
    except KbFetchError as x:
        print(i, "CHEBI ERROR", x); continue
    r = raw("chebi_raw", i, f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{i}/")
    kegg = []
    da = r.get("database_accessions") or {}
    for grp, lst in (da.items() if isinstance(da, dict) else []):
        for x in lst:
            if "KEGG" in str(x.get("source_name","")).upper() or "KEGG" in str(x.get("type","")).upper():
                kegg.append(x.get("accession_number"))
    try:
        p = pubchem_by_inchikey(c["inchikey"])
    except KbFetchError as x:
        p = {"error": str(x)[:120]}
    ok = p.get("inchikey","")[:14] == c["inchikey"][:14]
    print(i, "|", c["name"], "|", c["formula"], c["charge"], c["monoisotopic_mass"], f"{formula_mass(c['formula']):.4f}", "|", c["inchikey"], "| pubchem", p.get("cid"), p.get("formula"), p.get("exact_mass"), "block1 ok" if ok else p, "| kegg", sorted(set(kegg)), "| netcharge", net_charge(c["smiles"] or ""))
    out[i] = {"chebi": c, "pubchem": p, "kegg": sorted(set(kegg))}
    time.sleep(0.3)
json.dump(out, open("compounds.json","w"), indent=1)
r = raw("chebi_raw", "40968", "https://www.ebi.ac.uk/chebi/backend/api/public/compound/40968/")
print(list(r.keys())); print(json.dumps(r.get("database_accessions"))[:1500])
