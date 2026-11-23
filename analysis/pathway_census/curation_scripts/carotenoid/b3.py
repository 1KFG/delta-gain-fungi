"""Fix round 3: NCBI Taxonomy taxid of every organism name used by kb.species_key.

For each name: esearch (db=taxonomy, field All Names) gives the taxid(s); efetch gives the
scientific name and rank of that taxid. Writes tests/data/kb/ncbi_taxonomy_names.tsv.
"""
from h import *
import csv, time, urllib.parse
import xml.etree.ElementTree as ET
NAMES = ["Gibberella fujikuroi", "Fusarium fujikuroi", "Phaffia rhodozyma", "Xanthophyllomyces dendrorhous",
         "Mucor circinelloides f. lusitanicus", "Mucor racemosus var. lusitanicus", "Mucor lusitanicus",
         "Mucor circinelloides", "Mucor racemosus", "Mycosarcoma maydis", "Ustilago maydis",
         "Neurospora crassa", "Phycomyces blakesleeanus", "Blakeslea trispora", "Mucor mucedo"]
rows = []
for n in NAMES:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmode=json&term=" + urllib.parse.quote(f'"{n}"[All Names]')
    d = raw("ncbi_taxonomy_search", n, url)
    ids = d["esearchresult"]["idlist"]
    time.sleep(1)
    assert len(ids) == 1, (n, ids)
    x = raw("ncbi_taxonomy", ids[0], f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ids[0]}&retmode=xml", as_json=False)
    t = ET.fromstring(x).find("Taxon")
    rows.append([n, ids[0], t.findtext("ScientificName"), t.findtext("Rank"), "2026-11-21"])
    print(rows[-1]); time.sleep(1)
with open(PKG / "tests" / "data" / "kb" / "ncbi_taxonomy_names.tsv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["name", "taxid", "ncbi_scientific_name", "rank", "fetched_on"]); w.writerows(rows)
