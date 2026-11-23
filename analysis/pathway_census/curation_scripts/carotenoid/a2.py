from h import *
import sys
for a in ["Q7Z859"]:
    e = uniprot_entry(a)
    print(a, e['accession'], '|', e['organism'], '|', e['gene_names'], e['length'], e['ec'], 'rev' if e['reviewed'] else 'unrev', e['pmids'])
pm = "1826006 12695851 8163509 12172798 11172012 10951210 10589832 11862485 10397832 11017770 9079885 17293483 19584000 17493127 17610084 18627463 21749649 18769767 16416328 18837978 25504221 21854466 22445959 17822929 18931040 8828234 23893079 15345447".split()
with open("abs.txt","w") as fh:
    for p in pm:
        try:
            d = pubmed_abstract(p)
            fh.write(f"=== {p} ({d['year']}) {d['title']}\n{d['abstract']}\n\n")
        except KbFetchError as x:
            fh.write(f"=== {p} ERROR {x}\n\n"); print(p, "ERROR", x)
