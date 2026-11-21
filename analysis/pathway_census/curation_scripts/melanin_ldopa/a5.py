"""Exploration: pairwise global identity (BLOSUM62, gap open -10, extend -0.5) between
references and decoys. Identity = identical columns / alignment columns, in percent.
Usage: a5.py ACC ACC [ACC ...]. Needs Biopython."""
from h import *
from Bio import Align
from Bio.Align import substitution_matrices
al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = "global"
def ident(a, b):
    x = al.align(a, b)[0]
    s1, s2 = x[0], x[1]
    same = sum(1 for p, q in zip(s1, s2) if p == q and p != "-")
    return 100.0 * same / len(s1)
ACC = sys.argv[1:] or ["Q2UP46", "A0A1S9DK56", "Q00234", "P00440", "C7FF05", "C7FF04", "A7BHQ9", "Q2TL94", "Q0CRX0", "Q5BGU9", "Q5AUW8", "A0A8J9R8H5", "Q2UNF9", "B1VTI5"]
seq = {a: uniprot_entry(a)["sequence"] for a in ACC}
print("\t" + "\t".join(ACC))
for a in ACC:
    print(a + "\t" + "\t".join(f"{ident(seq[a], seq[b]):.0f}" for b in ACC))
