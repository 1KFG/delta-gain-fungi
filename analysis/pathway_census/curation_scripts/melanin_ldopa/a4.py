"""Exploration: ChEBI name search and Rhea search.
Usage: a4.py chebi 'term' ... | a4.py rhea 'query' ... | a4.py id CHEBI_NUM ..."""
from h import *
import urllib.parse

def chebi_search(term):
    url = "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/?term=" + urllib.parse.quote(term) + "&page=1&size=8"
    return raw("chebi_search", term, url)

def rhea(q, limit=10):
    url = "https://www.rhea-db.org/rhea/?query=" + urllib.parse.quote(q) + f"&columns=rhea-id,equation,chebi-id,ec&format=tsv&limit={limit}"
    return raw("rhea", q, url, as_json=False)

if __name__ == "__main__":
    mode = sys.argv[1]
    for t in sys.argv[2:]:
        print(f"## {mode} {t}")
        try:
            if mode == "chebi":
                for h in chebi_search(t).get("results", []):
                    s = h.get("_source", h)
                    print("  ", s.get("chebi_accession"), "|", s.get("name"), "|", s.get("formula"), "|", s.get("charge"), "| stars", s.get("stars"))
            elif mode == "rhea":
                print(rhea(t))
            else:
                c = chebi_compound(t)
                print("  ", c)
        except Exception as ex:
            print("  ERROR", repr(ex))
