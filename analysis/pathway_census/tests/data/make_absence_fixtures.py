"""Generate the synthetic fixtures for tests/test_absence.py (seeded, deterministic).

Run: /usr/bin/python3.12 tests/data/make_absence_fixtures.py tests/data
Writes:
  mini_step_refs.faa           2 random proteins: refA (150 aa), refB (150 aa)
  mini_genome.fna              2 contigs; refA back-translated (full length) planted mid-contig
  mini_genome_fragment.fna     refB back-translated, first 45 aa only (30% of the query), planted
  mini_genome_nohit.fna        random sequence only
  mini_genome_decoy.fna        refA back-translated after 55% of positions are substituted (diverged ORF)
  mini_genome_overlap.fna      refB[0:70] and refB[50:100] planted at two loci (overlapping query ranges)
  mini_genome_stop.fna         refA back-translated with a TAA stop codon inserted after residue 75
Proteins and flanks are random, not real. Back-translation uses one fixed codon per amino acid.
"""
import random
import sys
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
         "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
         "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}


def rprot(rng, n):
    return "".join(rng.choice(AA) for _ in range(n))


def rdna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def backtranslate(p):
    return "".join(CODON[a] for a in p)


def wrap(s, w=70):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def fa(records):
    return "".join(f">{k}\n{wrap(v)}\n" for k, v in records)


def main(d):
    d = Path(d)
    rng = random.Random(20261120)
    ref_a, ref_b = rprot(rng, 150), rprot(rng, 150)
    (d / "mini_step_refs.faa").write_text(fa([("refA", ref_a), ("refB", ref_b)]))
    planted = rdna(rng, 1500) + backtranslate(ref_a) + rdna(rng, 1500)
    (d / "mini_genome.fna").write_text(fa([("contig1", rdna(rng, 3000)), ("contig2", planted)]))
    frag = rdna(rng, 1500) + backtranslate(ref_b[:45]) + rdna(rng, 1500)
    (d / "mini_genome_fragment.fna").write_text(fa([("contig1", frag)]))
    (d / "mini_genome_nohit.fna").write_text(fa([("contig1", rdna(rng, 3000)), ("contig2", rdna(rng, 3000))]))

    # Added after the fixtures above so their random draws are unchanged.
    div = [rng.choice(AA) if rng.random() < 0.55 else a for a in ref_a]
    (d / "mini_genome_decoy.fna").write_text(fa([("contig1", rdna(rng, 1500) + backtranslate("".join(div)) + rdna(rng, 1500))]))
    ov = (rdna(rng, 800) + backtranslate(ref_b[0:70]) + rdna(rng, 1500)
          + backtranslate(ref_b[50:100]) + rdna(rng, 800))
    (d / "mini_genome_overlap.fna").write_text(fa([("contig1", ov)]))
    stop = rdna(rng, 1500) + backtranslate(ref_a[:75]) + "TAA" + backtranslate(ref_a[75:]) + rdna(rng, 1500)
    (d / "mini_genome_stop.fna").write_text(fa([("contig1", stop)]))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
