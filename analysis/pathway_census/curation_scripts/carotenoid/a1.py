from h import *
refs = "P24322 Q9P885 P37295 Q8X0Z1 Q9P854 Q9UUQ6 Q7Z859 P21334 Q8X0Z0 P54982 O13506 S0ENQ2 A1KQY3 A7UXI1 Q1K615 F6IBC7 B5LEC8 Q3HR17 B6D970 G8Z9M3 A0ABR3AZC2 A5HXP5 B9DR41 Q01213 Q7S860 Q67GH9".split()
for a in refs:
    try:
        e = uniprot_entry(a)
        print(a, e['accession'], '|', e['organism'], '|', e['gene_names'], e['length'], e['ec'], 'rev' if e['reviewed'] else 'unrev', e['pmids'])
    except KbFetchError as x:
        print(a, 'ERROR', x)
