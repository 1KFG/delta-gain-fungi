"""NCBI Taxonomy taxid of every organism name used by the DHN-melanin tables.

Not in the list: "Glomerella lagenarium" and "Colletotrichum lagenarium" (the UniProt name of
P79068 and the name in the papers). esearch returned no taxid for either name on 2026-11-21.
UniProt gives P79068 the taxid 5462, whose NCBI scientific name is Colletotrichum lagenaria.

For each name: esearch (db=taxonomy, field All Names) gives the taxid(s); efetch gives the
scientific name and rank. Prints the rows. With --write it appends the rows that are not yet in
tests/data/kb/ncbi_taxonomy_names.tsv (existing rows are kept as they are).
"""
from h import *
import csv, time, urllib.parse
import xml.etree.ElementTree as ET
NAMES = ["Colletotrichum lagenaria", "Colletotrichum orbiculare",
         "Exophiala dermatitidis", "Wangiella dermatitidis", "Aspergillus fumigatus",
         "Pyricularia oryzae", "Magnaporthe oryzae", "Pyricularia grisea", "Magnaporthe grisea",
         "Pestalotiopsis fici", "Sordaria macrospora", "Aspergillus niger"]
TSV = PKG / "tests" / "data" / "kb" / "ncbi_taxonomy_names.tsv"
rows = []
for n in NAMES:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmode=json&term=" + urllib.parse.quote(f'"{n}"[All Names]')
    d = raw("ncbi_taxonomy_search", n, url)
    ids = d["esearchresult"]["idlist"]
    time.sleep(0.6)
    if len(ids) != 1:
        print("NOT ONE TAXID", n, ids); continue
    x = raw("ncbi_taxonomy", ids[0], f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ids[0]}&retmode=xml", as_json=False)
    t = ET.fromstring(x).find("Taxon")
    rows.append([n, ids[0], t.findtext("ScientificName"), t.findtext("Rank"), "2026-11-21"])
    print(rows[-1]); time.sleep(0.6)
if "--write" in sys.argv:
    have = {ln.split("\t")[0] for ln in TSV.read_text(encoding="utf-8").splitlines()[1:]}
    with open(TSV, "a", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerows([r for r in rows if r[0] not in have])
