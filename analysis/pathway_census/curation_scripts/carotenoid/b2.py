"""Fix round 2 lookups: NCBI Taxonomy names for the organism aliases; GenBank OAD75744 against A0ABR3AZC2."""
from h import *
import re, time
import xml.etree.ElementTree as ET
for tid in ["5127", "1279085", "264483", "29924", "5270", "367110", "763407", "4842", "101091", "29922"]:
    x = raw("ncbi_taxonomy", tid, f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tid}&retmode=xml", as_json=False)
    t = ET.fromstring(x).find("Taxon")
    other = t.find("OtherNames")
    names = [(e.tag, e.text) for e in other if e.text and e.text.strip()] if other is not None else []
    names += [(n.findtext("ClassCDE"), n.findtext("DispName")) for n in (other.findall("Name") if other is not None else [])]
    print(tid, "|", t.findtext("ScientificName"), "|", t.findtext("Rank"), "| parent", t.findtext("ParentTaxId"), "|", names)
    time.sleep(1)
gp = raw("ncbi_protein", "OAD75744_gp", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=OAD75744&rettype=gp&retmode=text", as_json=False)
print(gp.split("REFERENCE")[0]); print([l for l in gp.splitlines() if "locus_tag" in l or "PUBMED" in l])
seq = re.sub(r"[^a-z]", "", gp.split("ORIGIN")[1].split("//")[0]).upper()
u = uniprot_entry("A0ABR3AZC2")
print("len", len(seq), u["length"], "identical:", seq == u["sequence"])
ab = pubmed_abstract("18931040")["abstract"]; print("activity sentence present:", "At the enzyme level, activity is found exclusively in the (-) mating type." in ab)
