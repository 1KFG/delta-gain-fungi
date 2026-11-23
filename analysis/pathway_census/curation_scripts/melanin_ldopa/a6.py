"""Exploration: PubMed search (esearch + esummary titles), paced at 2 requests per second or less.
Usage: a6.py 'query' ['query' ...]"""
from h import *
import time, urllib.parse

def pm_search(q, n=12):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=%d&term=%s" % (n, urllib.parse.quote(q))
    ids = raw("pubmed_search", q + f"|{n}", url)["esearchresult"]["idlist"]
    time.sleep(0.6)
    if not ids:
        return []
    s = raw("pubmed_summary", ",".join(ids), "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(ids))["result"]
    time.sleep(0.6)
    return [(i, s[i].get("pubdate", "")[:4], s[i].get("title", "")) for i in ids]

if __name__ == "__main__":
    for q in sys.argv[1:]:
        print(f"## {q}")
        for i, y, t in pm_search(q):
            print(f"  {i} ({y}) {t[:170]}")
