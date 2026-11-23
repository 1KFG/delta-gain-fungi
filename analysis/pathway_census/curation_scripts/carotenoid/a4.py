from h import *
import urllib.parse, time
F = "accession,reviewed,gene_names,organism_name,length,ec,protein_name"
def s(key, q, size=25):
    url = "https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) + f"&fields={F}&format=tsv&size={size}"
    t = raw("uniprot_search", key, url, as_json=False)
    print(f"--- {key}: {q}"); print(t)
s("cco1_um", "(gene:cco1) AND (taxonomy_id:5270)")
s("fpps_fungi", "(ec:2.5.1.10) AND (taxonomy_id:4751) AND (reviewed:true)")
s("sqs_fungi", "(ec:2.5.1.21) AND (taxonomy_id:4751) AND (reviewed:true)")
s("pf01593_fungi", "(xref:pfam-PF01593) AND (taxonomy_id:4751) AND (reviewed:true)", 60)
s("pf03055_rev", "(xref:pfam-PF03055) AND (reviewed:true) AND ((taxonomy_id:4751) OR (taxonomy_id:2))", 60)
s("aldh_nc_sc", "(xref:pfam-PF00171) AND (reviewed:true) AND ((taxonomy_id:559292) OR (taxonomy_id:367110) OR (taxonomy_id:227321))", 60)
