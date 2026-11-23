"""Exploration: fetch the open-access full text of a paper from Europe PMC (XML) through the
cached fetch, strip the tags, and print the passages that contain a search word.
Usage: fetch_fulltext.py PMID [word ...]. Prints 'NO FULL TEXT' when Europe PMC has none."""
from helpers import *
import re

def pmcid(pmid):
    d = raw("europepmc_search", f"ext_id_{pmid}",
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=lite")
    res = d["resultList"]["result"]
    return (res[0].get("pmcid"), res[0].get("isOpenAccess")) if res else (None, None)

def fulltext(pmid):
    pmc, oa = pmcid(pmid)
    if not pmc:
        return None, None
    x = raw("europepmc_fulltext", pmc, f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML", as_json=False)
    x = re.sub(r"<ref-list.*?</ref-list>", " ", x, flags=re.S)
    t = re.sub(r"<[^>]+>", " ", x)
    return pmc, " ".join(t.split())

if __name__ == "__main__":
    pmid, words = sys.argv[1], sys.argv[2:]
    try:
        pmc, t = fulltext(pmid)
    except KbFetchError as ex:
        print("NO FULL TEXT (fetch failed):", ex); sys.exit(0)
    if not t:
        print("NO FULL TEXT"); sys.exit(0)
    print(f"## PMID {pmid} = {pmc}, {len(t)} characters of text")
    for w in words:
        for m in re.finditer(re.escape(w), t, flags=re.I):
            print(f"[{w}] ...{t[max(0, m.start()-300):m.end()+400]}...\n")
