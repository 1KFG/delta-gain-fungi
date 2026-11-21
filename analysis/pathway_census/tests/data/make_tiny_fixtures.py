"""Generate the synthetic fixtures for tests/test_hits.py (seeded, deterministic).

Run: /usr/bin/python3.12 tests/data/make_tiny_fixtures.py tests/data
Sequences are random; they are not real proteins.
"""
import random
import sys
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"


def rnd(rng, n):
    return "".join(rng.choice(AA) for _ in range(n))


def mutate(rng, s, frac):
    out = list(s)
    for i in rng.sample(range(len(s)), int(len(s) * frac)):
        out[i] = rng.choice(AA)
    return "".join(out)


def fa(records):
    return "".join(f">{k}\n{v}\n" for k, v in records)


def main(d):
    d = Path(d)
    rng = random.Random(20261120)
    base = rnd(rng, 150)
    refs = [(f"ref{i}", mutate(rng, base, 0.15)) for i in range(1, 6)]
    (d / "tiny_refs.faa").write_text(fa(refs))
    prot = [("homologue", mutate(rng, base, 0.2)), ("decoy", rnd(rng, 150))]
    prot += [(f"filler{i}", rnd(rng, 120)) for i in range(5)]
    (d / "tiny_proteome.faa").write_text(fa(prot))
    rng2 = random.Random(7)
    (d / "unrelated_proteome.faa").write_text(fa([(f"unrel{i}", rnd(rng2, 130)) for i in range(8)]))
    # dereplication fixture: three near-identical proteins plus one distinct
    b2 = rnd(rng, 200)
    (d / "derep_in.faa").write_text(fa([("a", b2), ("b", mutate(rng, b2, 0.01)),
                                        ("c", mutate(rng, b2, 0.01)), ("d", rnd(rng, 200))]))


if __name__ == "__main__":
    main(sys.argv[1])
