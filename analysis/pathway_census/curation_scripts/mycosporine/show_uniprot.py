"""Exploration: print the checked fields of UniProt entries. Usage: show_uniprot.py ACC [ACC ...]"""
from helpers import *

def full(acc):
    return raw("uniprot_full", acc, f"https://rest.uniprot.org/uniprotkb/{acc}.json")

def show(acc):
    d = full(acc)
    pdn = d["proteinDescription"]
    nm = (pdn.get("recommendedName") or (pdn.get("submissionNames") or [{}])[0])["fullName"]["value"]
    e = uniprot_entry(acc)
    print(f"== {acc} [{'reviewed' if e['reviewed'] else 'unreviewed'}] {d['uniProtkbId']} | {nm}")
    print(f"   organism: {e['organism']} (taxid {d['organism']['taxonId']}) | genes: {e['gene_names']} | len {e['length']} | EC {e['ec']}")
    au = d["entryAudit"]
    print(f"   entry v{au['entryVersion']} seq v{au['sequenceVersion']} | pmids: {e['pmids']}")
    for c in d.get("comments", []):
        if c["commentType"] == "CATALYTIC ACTIVITY":
            r = c["reaction"]
            xr = [x["id"] for x in r.get("reactionCrossReferences", []) if x["database"] == "Rhea"]
            print(f"   CA: {r['name']} | EC {r.get('ecNumber')} | {xr[:2]}")
        elif c["commentType"] in ("FUNCTION", "PATHWAY", "DISRUPTION PHENOTYPE"):
            for t in c.get("texts", []):
                print(f"   {c['commentType']}: {t['value'][:450]}")
    pf = [(x["id"], [p["value"] for p in x.get("properties", []) if p["key"] == "EntryName"]) for x in d.get("uniProtKBCrossReferences", []) if x["database"] in ("Pfam", "NCBIfam", "PROSITE")]
    print(f"   Pfam/NCBIfam/PROSITE: {pf}")
    # kbsources.uniprot_entry does not read ordered locus names; they are printed from the full record
    print(f"   ordered locus names: {[o['value'] for g in d.get('genes', []) for o in g.get('orderedLocusNames', [])]}")
    for r in d.get("references", []):
        c = r["citation"]
        pm = [x["id"] for x in c.get("citationCrossReferences", []) or [] if x["database"] == "PubMed"]
        print(f"   ref: {pm} {c.get("title","")[:90]} | {r.get("referencePositions")}") if len(d.get("references", [])) <= 8 or any(k in " ".join(r.get("referencePositions") or []) for k in ("CATALYTIC", "BIOPHYS")) else None

if __name__ == "__main__":
    for a in sys.argv[1:]:
        try:
            show(a)
        except Exception as ex:
            print(f"== {a} ERROR {ex!r}")
