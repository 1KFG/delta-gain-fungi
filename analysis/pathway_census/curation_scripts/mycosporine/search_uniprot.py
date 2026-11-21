"""Exploration: UniProt search. Usage: search_uniprot.py 'query' ['query' ...]"""
from helpers import *
import urllib.parse

def search(q, size=15):
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) +
           f"&fields=accession,reviewed,id,protein_name,gene_names,organism_name,length,ec,lit_pubmed_id&format=tsv&size={size}")
    return raw("uniprot_search", q + f"|{size}", url, as_json=False)

if __name__ == "__main__":
    for q in sys.argv[1:]:
        print(f"## {q}")
        try:
            print(search(q))
        except Exception as ex:
            print("ERROR", repr(ex))
