#!/usr/bin/env python3
"""Summarise and check the full BFD search for melanin_dhn, melanin_ldopa, mycosporine.

Usage: melanin_myco_step4_summary.py <search_dir> <manifests_dir> <out_dir>
Reads candidates.chunk*.zst and searched.chunk*.tsv. Writes small tables to <out_dir>.
Python 3.12 stdlib only. zstd is read by shelling out to `zstd -dc`.
"""
import csv, subprocess, sys, glob, os, collections, statistics

PROJ = {  # step -> pathway
    **{s: "melanin_dhn" for s in "ayg1 heptaketide_pks hexaketide_pks laccase_abr2 mco_abr1 scytalone_dehydratase t3hn_reductase t4hn_reductase thn_synthase".split()},
    **{s: "melanin_ldopa" for s in "asp_melanin_tyrosinase aspulvinone_e_synthase hmga hppd laccase_diphenol tyrosinase".split()},
    **{s: "mycosporine" for s in "atp_grasp_ligase ddgs omt".split()},
}
PROJECTION = {"tyrosinase": 51310, "hppd": 18000, "omt": 14610, "ayg1": 10169, "atp_grasp_ligase": 10052,
              "ddgs": 9467, "t4hn_reductase": 8298, "scytalone_dehydratase": 7480, "hexaketide_pks": 7363,
              "laccase_abr2": 7363, "thn_synthase": 6195, "t3hn_reductase": 2454, "heptaketide_pks": 1753,
              "hmga": 1753, "asp_melanin_tyrosinase": 1519, "mco_abr1": 1519, "laccase_diphenol": 468,
              "aspulvinone_e_synthase": 0}
PKS = ("heptaketide_pks", "hexaketide_pks", "thn_synthase")
OMT_HIGH = 120.7


def zlines(path):
    p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
    yield from p.stdout
    p.wait()


def main():
    sdir, mdir, odir = sys.argv[1:4]
    os.makedirs(odir, exist_ok=True)
    seen = collections.Counter()
    nprot = 0
    nbad = 0
    status = collections.Counter()
    per_chunk = {}
    for f in sorted(glob.glob(os.path.join(sdir, "searched.chunk*.tsv"))):
        n = 0
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("locustag"):
                    continue
                c = line.rstrip("\n").split("\t")
                seen[c[0]] += 1
                status[c[2]] += 1
                nprot += int(c[3])
                nbad += int(c[4])
                n += 1
        per_chunk[os.path.basename(f)] = n
    expected = set()
    for line in zlines(os.path.join(mdir, "genomes.tsv.zst")):
        c = line.rstrip("\n").split("\t")
        if c[0] != "locustag" and c[14] == "True":
            expected.add(c[0])
    print("searched files", per_chunk)
    print("unique genomes searched", len(seen), "duplicates", sum(1 for v in seen.values() if v > 1))
    print("status", dict(status))
    print("expected (has_proteins) genomes", len(expected), "missing", len(expected - set(seen)), "extra", len(set(seen) - expected))
    print("proteins searched", nprot, "bad id prefixes", nbad)

    rows = collections.Counter()
    prots = collections.defaultdict(set)
    genomes = collections.defaultdict(set)
    maxscore, minscore = {}, {}
    scores = collections.defaultdict(list)
    prot_steps = collections.defaultdict(set)
    omt_low_genomes = collections.Counter()
    omt_low = omt_high = 0
    tyr_per_genome = collections.Counter()
    all_cand_genomes = set()
    for f in sorted(glob.glob(os.path.join(sdir, "candidates.chunk*.zst"))):
        for line in zlines(f):
            if line.startswith("#") or line.startswith("protein_id"):
                continue
            c = line.rstrip("\n").split("\t")
            pid, step, score = c[0], c[1], float(c[2])
            g = pid.split("_")[0]
            rows[step] += 1
            prots[step].add(pid)
            genomes[step].add(g)
            all_cand_genomes.add(g)
            prot_steps[pid].add(step)
            maxscore[step] = max(maxscore.get(step, -1), score)
            minscore[step] = min(minscore.get(step, 1e9), score)
            scores[step].append(score)
            if step == "omt":
                if score < OMT_HIGH:
                    omt_low += 1
                    omt_low_genomes[g] += 1
                else:
                    omt_high += 1
            if step == "tyrosinase":
                tyr_per_genome[g] += 1
    tot = sum(rows.values())
    proj_tot = sum(PROJECTION.values())
    with open(os.path.join(odir, "candidates_per_step.tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["pathway", "step", "rows", "unique_proteins", "genomes_with_candidate", "min_score_seen", "median_score", "max_score", "projection_rows", "ratio_observed_to_projection"])
        for s in sorted(PROJ, key=lambda s: (PROJ[s], -rows[s])):
            r = rows[s]
            w.writerow([PROJ[s], s, r, len(prots[s]), len(genomes[s]),
                        minscore.get(s, ""), statistics.median(scores[s]) if scores[s] else "",
                        maxscore.get(s, ""), PROJECTION[s], f"{r / PROJECTION[s]:.2f}" if PROJECTION[s] else "NA"])
        w.writerow(["ALL", "total", tot, len(prot_steps), len(all_cand_genomes), "", "", "", proj_tot, f"{tot / proj_tot:.2f}"])
    print("total rows", tot, "unique proteins", len(prot_steps), "genomes with >=1 candidate", len(all_cand_genomes), "ratio", round(tot / proj_tot, 3))
    for p in ("melanin_dhn", "melanin_ldopa", "mycosporine"):
        st = [s for s in PROJ if PROJ[s] == p]
        up = set().union(*[prots[s] for s in st])
        print(p, "rows", sum(rows[s] for s in st), "unique proteins", len(up),
              "genomes", len(set().union(*[genomes[s] for s in st])))
    for s in sorted(PROJ, key=lambda s: -rows[s]):
        ratio = rows[s] / PROJECTION[s] if PROJECTION[s] else float("nan")
        print(f"{s:26s} rows {rows[s]:7d} prot {len(prots[s]):7d} genomes {len(genomes[s]):6d} max {maxscore.get(s)} proj {PROJECTION[s]} ratio {ratio:.2f}")
    print("omt rows 100.6<=s<120.7", omt_low, "genomes", len(omt_low_genomes), "; >=120.7", omt_high)
    cnt = collections.Counter()
    for pid, ss in prot_steps.items():
        k = len(ss & set(PKS))
        if k:
            cnt[k] += 1
    print("PKS proteins by number of PKS steps hit", dict(cnt), "total", sum(cnt.values()))
    pair = collections.Counter(tuple(sorted(ss & set(PKS))) for ss in prot_steps.values() if len(ss & set(PKS)) >= 2)
    print("PKS combos", dict(pair))
    for grp, name in ((("t3hn_reductase", "t4hn_reductase"), "SDR"),
                      (("laccase_abr2", "laccase_diphenol", "mco_abr1"), "laccase/MCO"),
                      (("tyrosinase", "asp_melanin_tyrosinase"), "tyrosinase pair"),
                      (("tyrosinase", "laccase_diphenol", "mco_abr1", "laccase_abr2"), "tyr+laccase/MCO")):
        n = sum(1 for ss in prot_steps.values() if len(ss & set(grp)) >= 2)
        print("proteins hitting >=2 of", name, grp, n)
    d = collections.Counter(tyr_per_genome.values())
    print("tyrosinase rows per genome (genomes with >=1):", sorted(d.items()), "genomes", len(tyr_per_genome),
          "max", max(tyr_per_genome.values()) if tyr_per_genome else 0)
    with open(os.path.join(odir, "tyrosinase_rows_per_genome.tsv"), "w") as o:
        o.write("rows_per_genome\tgenomes\n")
        for k in sorted(d):
            o.write(f"{k}\t{d[k]}\n")
    with open(os.path.join(odir, "omt_low_band_genomes.tsv"), "w") as o:
        o.write("locustag\tomt_rows_100.6_to_120.7\n")
        for g, n in sorted(omt_low_genomes.items()):
            o.write(f"{g}\t{n}\n")


main()
