"""Exploration: fetch PubMed abstracts (paced at 2 per second or less) and print them.
Usage: show_pubmed.py PMID [PMID ...]"""
from helpers import *
import time
from pathway_census.kbsources import CACHE_DIR

if __name__ == "__main__":
    for p in sys.argv[1:]:
        if not (CACHE_DIR / "pubmed" / f"{p}.json").exists():
            time.sleep(0.6)
        try:
            a = pubmed_abstract(p)
            print(f"## PMID {p} ({a['year']}) {a['title']}\n{a['abstract']}\n")
        except Exception as ex:
            print(f"## PMID {p} ERROR {ex!r}\n")
