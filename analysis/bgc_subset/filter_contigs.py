#!/usr/bin/env python3
"""Keep contigs of at least MINLEN bp in a FASTA and in the matching GFF3.
Usage: filter_contigs.py MINLEN IN.fa IN.gff3 OUT.fa OUT.gff3 [FAI]
Contig lengths come from the .fai (name, length) if given, else from the FASTA. Prints counts.
The same contig set is used for both files, so antiSMASH sees no GFF3 feature without a sequence."""
import sys
minlen = int(sys.argv[1]); fa, gff, ofa, ogff = sys.argv[2:6]
fai = sys.argv[6] if len(sys.argv) > 6 else None
keep = set(); total = 0
if fai:
    for l in open(fai):
        a = l.split("\t"); total += 1
        if int(a[1]) >= minlen: keep.add(a[0])
else:
    name = None; ln = 0; L = {}
    for l in open(fa):
        if l.startswith(">"):
            if name: L[name] = ln
            name = l[1:].split()[0]; ln = 0
        else: ln += len(l.strip())
    if name: L[name] = ln
    total = len(L); keep = {k for k, v in L.items() if v >= minlen}
w = False; nout = 0
with open(fa) as fin, open(ofa, "w") as fo:
    for l in fin:
        if l.startswith(">"):
            w = l[1:].split()[0] in keep; nout += w
        if w: fo.write(l)
ng = 0
with open(gff) as fin, open(ogff, "w") as fo:
    for l in fin:
        if l.startswith("##sequence-region"):
            if l.split()[1] in keep: fo.write(l)
        elif l.startswith("#"): fo.write(l)
        else:
            if l.split("\t", 1)[0] in keep: fo.write(l); ng += 1
print(f"minlen {minlen}: kept {len(keep)} of {total} contigs; wrote {nout} fasta records, {ng} gff3 lines", file=sys.stderr)
