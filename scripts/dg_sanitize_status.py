import csv,sys
p=sys.argv[1]
with open(p,newline="") as fh: rows=list(csv.reader(fh,delimiter="\t"))
rows=[[c.strip() for c in r] for r in rows]
with open(p,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerows(rows)
print(len(rows)-1,"rows")
