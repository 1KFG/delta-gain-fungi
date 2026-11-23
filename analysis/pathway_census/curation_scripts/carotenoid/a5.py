from h import *
import urllib.parse
F = "accession,reviewed,gene_names,organism_name,length,ec,protein_name"
def s(key, q, size=25):
    url = "https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) + f"&fields={F}&format=tsv&size={size}"
    t = raw("uniprot_search", key, url, as_json=False)
    print(f"--- {key}: {q}"); print(t)
s("um_pf03055", "(xref:pfam-PF03055) AND (taxonomy_id:5270)")
s("p450_sc_nc", "(xref:pfam-PF00067) AND (reviewed:true) AND ((taxonomy_id:559292) OR (taxonomy_id:367110) OR (taxonomy_id:5507) OR (taxonomy_id:5346))", 40)
s("tigr03462_rev", "(xref:ncbifam-TIGR03462) AND (reviewed:true)", 40)
s("pf00667_fungi", "(xref:pfam-PF00667) AND (reviewed:true) AND ((taxonomy_id:559292) OR (taxonomy_id:284812))", 40)
s("csga", "(gene:csgA) AND (taxonomy_id:34) AND (reviewed:true)")
s("akr_sc", "(xref:pfam-PF00248) AND (reviewed:true) AND (taxonomy_id:559292)", 40)
s("coq1", "(gene:COQ1) AND (taxonomy_id:559292) AND (reviewed:true)")
