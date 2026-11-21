"""Exploration: for each compound name, ask PubChem for the compounds of that name (CID, InChIKey,
formula, charge), then search ChEBI by the InChIKey and by its first block. A name with no ChEBI
hit is a USER DECISION (no id is invented). Usage: find_compounds.py 'name' ['name' ...]"""
from helpers import *
import urllib.parse, time

def pubchem_name(name):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(name) +
           "/property/InChIKey,MolecularFormula,ExactMass,Charge,IUPACName/JSON")
    return raw("pubchem_name", name, url)["PropertyTable"]["Properties"]

def chebi_by_text(term):
    url = "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/?term=" + urllib.parse.quote(term) + "&page=1&size=8"
    try:
        return raw("chebi_search", term, url).get("results", [])
    except KbFetchError as ex:
        if "404" in str(ex):   # the ChEBI search answers 404 when it has no hit
            return []
        raise

if __name__ == "__main__":
    for name in sys.argv[1:]:
        print(f"## {name}")
        try:
            props = pubchem_name(name)
        except KbFetchError as ex:
            print("   PubChem: no compound of this name", str(ex)[-40:]); continue
        for p in props[:3]:
            hits = chebi_by_text(p["InChIKey"]) or chebi_by_text(p["InChIKey"][:14])
            hs = [(h["_source"]["chebi_accession"], h["_source"].get("name"), h["_source"].get("formula"), h["_source"].get("charge")) for h in hits]
            print(f"   CID {p['CID']} {p['InChIKey']} {p['MolecularFormula']} charge {p.get('Charge')} | ChEBI: {hs or 'NONE'}")
            time.sleep(0.4)
