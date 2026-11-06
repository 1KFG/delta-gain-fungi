import csv
bad=[]
with open("prior_annotation_status.tsv",newline="") as fh:
    for r in csv.DictReader(fh,delimiter="\t"):
        if "\r" in r["folder"]: bad.append((r["accession"],repr(r["folder"]),r["status"]))
print(len(bad)); [print(*b) for b in bad[:6]]
import collections; print(collections.Counter(b[2] for b in bad))
