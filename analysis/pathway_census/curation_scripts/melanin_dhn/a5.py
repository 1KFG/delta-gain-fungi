"""Exploration: pairwise global identity (BLOSUM62, gap open -10, extend -0.5) between
reductase references and decoys. Identity = identical columns / alignment columns.
Also prints KEGG entry names for the cross-referenced KOs and compounds (names only)."""
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
ACC = sys.argv[1:] or ["Q9HFV6", "E9QUT3", "W3XAA0", "F8WTI4", "A7DWM4", "Q12634", "P87025", "A7DWM6", "P50161", "Q00791", "Q8TFD5"]
seq = {a: uniprot_entry(a)["sequence"] for a in ACC}
print("\t" + "\t".join(ACC))
for a in ACC:
    print(a + "\t" + "\t".join(f"{ident(seq[a], seq[b]):.0f}" for b in ACC))
