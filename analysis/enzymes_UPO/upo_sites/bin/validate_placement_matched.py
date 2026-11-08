#!/usr/bin/env python3.12
"""
Validate the placement clade caller AT MATCHED IDENTITY.

The problem this addresses, stated in PLAN.md section 3.1 item 3: ordinary leave-one-out tests the
caller only in the regime where references are close to each other (52.7-96.2 % identity), while
BFD proteins sit at 25-35 % identity to anything characterized. An accuracy figure from the easy
regime says nothing about the hard one.

Method. For each labelled reference q and each identity ceiling T:
  1. remove q from the reference set;
  2. ALSO remove every reference whose identity to q exceeds T;
  3. require at least one reference of q's own clade to survive, otherwise the case is impossible
     and is excluded (and counted separately);
  4. rebuild the tree on what is left, place q with EPA, and ask whether its clade is recovered.

Step 2 is the point: it forces the caller to work at the divergence BFD actually presents.
Accuracy is then reported as a function of the identity to the NEAREST SURVIVING reference, which
is the quantity that can be compared directly against BFD queries.
"""
from __future__ import annotations
import argparse, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--ref-aln", required=True, help="aligned references (ref_aln.filtered.fa)")
ap.add_argument("--labels", required=True)
ap.add_argument("--thresholds", default="95,70,55,45,40,35")
ap.add_argument("--out", required=True)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--query-list", default=None,
                help="file of tip names to use as QUERIES. Everything else stays in the reference "
                     "set. Used to test whether an EXPANDED reference set recovers the ORIGINAL "
                     "characterized enzymes at low identity.")
a = ap.parse_args()


def rd(p):
    d, n = {}, None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = ""
        else:
            d[n] += l.strip()
    return d


aln = rd(a.ref_aln)
lab = {}
with open(a.labels) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            lab[f[0]] = f[1]
names = list(aln)
print(f"{len(names)} references, {len(lab)} labelled", file=sys.stderr)


def ident(x, y):
    n = s = 0
    for p, q in zip(x, y):
        if p == "-" and q == "-":
            continue
        n += 1
        if p == q and p != "-":
            s += 1
    return 100.0 * s / n if n else 0.0


ID = {}
for i, x in enumerate(names):
    for y in names[i + 1:]:
        v = ident(aln[x], aln[y]); ID[(x, y)] = v; ID[(y, x)] = v

thresholds = [float(t) for t in a.thresholds.split(",")]
rows = []
work = tempfile.mkdtemp(prefix="matched.")

qlist = None
if a.query_list:
    qlist = {l.strip() for l in open(a.query_list) if l.strip()}
    print(f"restricting queries to {len(qlist)} tips", file=sys.stderr)

for q in sorted(lab):
    if qlist is not None and q not in qlist:
        continue
    tc = lab[q]
    for T in thresholds:
        keep = [n for n in names if n != q and ID[(q, n)] <= T]
        same = [n for n in keep if lab.get(n) == tc]
        if len(keep) < 4 or not same:
            rows.append((q, tc, T, None, None, None, "impossible_no_same_clade_ref"))
            continue
        nearest = max(ID[(q, n)] for n in keep)
        d = os.path.join(work, f"{q}_{int(T)}")
        os.makedirs(d, exist_ok=True)
        refp = os.path.join(d, "ref.fa")
        with open(refp, "w") as fh:
            for n in keep:
                fh.write(f">{n}\n{aln[n]}\n")
        allp = os.path.join(d, "all.fa")
        shutil.copy(refp, allp)
        with open(allp, "a") as fh:
            fh.write(f">QUERY\n{aln[q]}\n")
        labp = os.path.join(d, "labels.tsv")
        with open(labp, "w") as fh:
            fh.write("tip\tclade\n")
            for n in keep:
                if n in lab:
                    fh.write(f"{n}\t{lab[n]}\n")
        r = subprocess.run(["raxml-ng", "--search1", "--msa", refp, "--model", "LG+G4",
                            "--prefix", os.path.join(d, "t"), "--threads", "1",
                            "--seed", "1337", "--force", "perf_threads"],
                           capture_output=True, text=True)
        tre = os.path.join(d, "t.raxml.bestTree")
        if not os.path.exists(tre):
            rows.append((q, tc, T, nearest, None, None, "tree_failed")); continue
        cwd = os.getcwd(); os.chdir(d)
        subprocess.run(["raxmlHPC-PTHREADS", "-T", str(a.threads), "-f", "v", "-s", "all.fa",
                        "-t", "t.raxml.bestTree", "-m", "PROTGAMMALG", "-n", "epa", "-p", "1337"],
                       capture_output=True, text=True)
        os.chdir(cwd)
        jp = os.path.join(d, "RAxML_portableTree.epa.jplace")
        if not os.path.exists(jp):
            rows.append((q, tc, T, nearest, None, None, "epa_failed")); continue
        outp = os.path.join(d, "call.tsv")
        subprocess.run(["/usr/bin/python3.12",
                        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "jplace_to_clade.py"),
                        "--jplace", jp, "--labels", labp, "--out", outp],
                       capture_output=True, text=True)
        call, wt = "unassigned", 0.0
        if os.path.exists(outp):
            with open(outp) as fh:
                fh.readline()
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f[0] == "QUERY":
                        call, wt = f[1], float(f[2])
        rows.append((q, tc, T, nearest, call, wt,
                     "correct" if call == tc else ("declined" if call == "unassigned" else "WRONG")))
        print(f"  {q:12} T={T:5.1f} nearest={nearest:5.1f} true={tc:3} call={call:11} "
              f"w={wt:.2f} {rows[-1][6]}", file=sys.stderr)

with open(a.out, "w") as fh:
    fh.write("query\ttrue_clade\tthreshold\tnearest_surviving_identity\tcall\tweight\toutcome\n")
    for r in rows:
        fh.write("\t".join("" if x is None else (f"{x:.1f}" if isinstance(x, float) else str(x))
                           for x in r) + "\n")

print("\n=== accuracy by identity to nearest surviving reference ===", file=sys.stderr)
bins = [(0, 35), (35, 45), (45, 55), (55, 70), (70, 101)]
for lo, hi in bins:
    sel = [r for r in rows if r[3] is not None and lo <= r[3] < hi and r[6] in
           ("correct", "declined", "WRONG")]
    if not sel:
        continue
    c = sum(1 for r in sel if r[6] == "correct")
    d = sum(1 for r in sel if r[6] == "declined")
    w = sum(1 for r in sel if r[6] == "WRONG")
    print(f"  nearest {lo:3}-{hi:3}% id  n={len(sel):3}  correct={c:3}  declined={d:3}  WRONG={w:3}"
          f"   precision(when it commits)="
          f"{(100*c/(c+w)) if (c+w) else float('nan'):.0f}%", file=sys.stderr)
shutil.rmtree(work, ignore_errors=True)
