"""NCBI Taxonomy taxid of the organism names used by the mycosporine tables.

For each name: esearch (db=taxonomy, field All Names) gives the taxid(s); efetch gives the
scientific name and rank. Prints the rows. With --write it appends the rows that are not yet in
tests/data/kb/ncbi_taxonomy_names.tsv (existing rows are kept as they are).

- "Trichormus variabilis" (UniProt name of Q3M6C3) and "Anabaena variabilis" (the name in the
  papers) have one taxid (264691), so "anabaena variabilis" is an alias of the key
  "trichormus variabilis" in kb.SPECIES_ALIASES.
- Names in CHECK_ONLY are printed and never written: "Phaffia rhodozyma" and "Xanthophyllomyces
  dendrorhous" are in the table already (alias from the carotenoid task); "Naganishia
  friedmannii" and "Cryptococcus friedmannii" occur in the lineage table only, not as a key.
"""
from helpers import *
import csv, time, urllib.parse
import xml.etree.ElementTree as ET
NAMES = ["Nostoc punctiforme", "Trichormus variabilis", "Anabaena variabilis"]
CHECK_ONLY = ["Phaffia rhodozyma", "Xanthophyllomyces dendrorhous", "Naganishia friedmannii", "Cryptococcus friedmannii"]
TSV = PKG / "tests" / "data" / "kb" / "ncbi_taxonomy_names.tsv"
rows = []
for n in NAMES + CHECK_ONLY:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmode=json&term=" + urllib.parse.quote(f'"{n}"[All Names]')
    d = raw("ncbi_taxonomy_search", n, url)
    ids = d["esearchresult"]["idlist"]
    time.sleep(0.6)
    if len(ids) != 1:
        print("NOT ONE TAXID", n, ids); continue
    x = raw("ncbi_taxonomy", ids[0], f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ids[0]}&retmode=xml", as_json=False)
    t = ET.fromstring(x).find("Taxon")
    row = [n, ids[0], t.findtext("ScientificName"), t.findtext("Rank"), "2026-11-21"]
    print(("check only: " if n in CHECK_ONLY else "") + str(row)); time.sleep(0.6)
    if n not in CHECK_ONLY:
        rows.append(row)
if "--write" in sys.argv:
    have = {ln.split("\t")[0] for ln in TSV.read_text(encoding="utf-8").splitlines()[1:]}
    with open(TSV, "a", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerows([r for r in rows if r[0] not in have])
