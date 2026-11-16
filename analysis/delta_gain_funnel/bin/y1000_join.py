#!/usr/bin/env python3
"""How much of the enriched-Pfam novelty is Y1000+ deposition rather than new biology?

The delta-gain funnel calls a protein cluster "novel" when it has no member from a genome
that carried a prior NCBI annotation. That is a *deposition* boundary, not a knowledge
boundary, and the funnel's own limitations say so. This script measures the size of the gap
for the one result where it matters most.

The result in question: Pfam families selected for being enriched in new genomes show a
+14.6 to +14.9 point excess of novel clusters over an order-matched baseline null, against
+2.2 for size-matched random Pfam controls. The funnel notes the enriched list "is a
Saccharomycotina yeast list". The open question was whether those yeasts are new biology or
new *venue* - genomes whose assemblies went to GenBank without annotation and whose
proteins were published somewhere nr does not index.

The test is a join on NCBI BioProject. Y1000+ (PRJNA736342, "Y1000+ Saccharomycotina
yeasts genome sequencing and assembly", UW-Madison) and its 2018 precursor (PRJNA429441,
196 yeast species, also UW-Madison) published proteomes through the project's own releases
rather than through NCBI's annotation pipeline, so every one of their genomes lands on the
funnel's "new" side by construction while the proteins have been public for years.

For each family and clustering tier this reports:

  clusters            clusters at that tier
  novel_now           clusters with no baseline member      <- the funnel's definition
  novel_if_y_known    clusters with no baseline AND no Y1000+ member
  y_only              clusters that are novel now and stop being novel if Y1000+ counts
                      as known - i.e. novelty that is purely a deposition artefact

Run it over the enriched families and over the random-control families; the control is what
says whether any shrinkage is specific to the enriched list or applies everywhere.

Memory: families are processed one at a time and only the three needed columns are read.
Doing all 204 at once OOM-killed a login node, which is also why this has a slurm driver.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

TIERS = ("c90", "c70", "c50", "c30")


def family_names(null_dir, prefix="pfam_"):
    """Family ids, taken from the null run's outputs.

    The enriched list is not stored as a list anywhere; the families that were *run* are the
    list. `<fam>.ordnull.tsv` is one file per family, so their basenames are the set.
    """
    out = set()
    for p in glob.glob(os.path.join(null_dir, f"{prefix}*.ordnull.tsv")):
        out.add(os.path.basename(p).split(".ordnull")[0])
    return sorted(out)


def y1000_asmids(con, samples_parquet, bioprojects):
    """ASMIDs whose BIOPROJECT field names any of `bioprojects`.

    BIOPROJECT holds one or more accessions joined by ';' (a genome can belong to several
    umbrella projects), so this matches as substring on accession boundaries rather than by
    equality - `= 'PRJNA736342'` silently misses 'PRJEB65317;PRJNA736342'.
    """
    if not bioprojects:
        return set()
    terms = " or ".join(
        f"';' || BIOPROJECT || ';' like '%;{b};%'" for b in bioprojects
    )
    rows = con.execute(
        f"select ASMID from read_parquet('{samples_parquet}') where {terms}"
    ).fetchall()
    return {r[0] for r in rows}


def family_rows(con, fam, regions_dir, clusters_dir, y_table="y"):
    """-> [{family, tier, clusters, novel_now, novel_if_y_known, y_only, ...}] for one family.

    Returns [] when either input is missing: the wide run and the control run do not cover
    exactly the same families, and a family with a null result but no regions file is not an
    error, it is a family that was dropped upstream.
    """
    mem = os.path.join(regions_dir, f"{fam}.members.tsv")
    clu = os.path.join(clusters_dir, f"{fam}.clusters.tsv")
    if not (os.path.exists(mem) and os.path.exists(clu)):
        return []

    con.execute("drop table if exists m")
    con.execute("drop table if exists cl")
    con.execute(
        f"create table m as select region_id, asmid, class from "
        f"read_csv('{mem}', delim='\t', header=true, all_varchar=true)"
    )
    con.execute(
        f"create table cl as select * from "
        f"read_csv('{clu}', delim='\t', header=true, all_varchar=true)"
    )
    have = {r[0] for r in con.execute("describe cl").fetchall()}

    out = []
    for tier in TIERS:
        if tier not in have:
            continue
        r = con.execute(f"""
            with j as (
              select cl.{tier} as cid, m.class,
                     (yy.ASMID is not null) as is_y
              from cl join m using(region_id)
              left join {y_table} yy on m.asmid = yy.ASMID
            ), g as (
              select cid,
                     max(case when class='baseline' then 1 else 0 end) as has_base,
                     max(case when class='baseline' or is_y then 1 else 0 end) as has_base_or_y
              from j group by cid
            )
            select count(*), sum(1-has_base), sum(1-has_base_or_y) from g
        """).fetchone()
        clusters, novel_now, novel_if = int(r[0]), int(r[1] or 0), int(r[2] or 0)
        out.append({
            "family": fam, "tier": tier, "clusters": clusters,
            "novel_now": novel_now, "novel_if_y_known": novel_if,
            "y_only": novel_now - novel_if,
            "novel_now_pct": round(100.0 * novel_now / clusters, 4) if clusters else "",
            "novel_if_y_known_pct": round(100.0 * novel_if / clusters, 4) if clusters else "",
        })
    return out


def summarize(rows, tier):
    """Pooled totals for one tier: clusters, the two novel counts, and the share lost."""
    sel = [r for r in rows if r["tier"] == tier]
    tot = sum(r["clusters"] for r in sel)
    now = sum(r["novel_now"] for r in sel)
    aft = sum(r["novel_if_y_known"] for r in sel)
    return {
        "families": len(sel), "clusters": tot, "novel_now": now,
        "novel_if_y_known": aft, "y_only": now - aft,
        "novel_now_pct": (100.0 * now / tot) if tot else 0.0,
        "novel_if_y_known_pct": (100.0 * aft / tot) if tot else 0.0,
        "pct_of_novel_lost": (100.0 * (now - aft) / now) if now else 0.0,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--null-dir", required=True, help="dir of <fam>.ordnull.tsv; defines the family set")
    ap.add_argument("--regions-dir", required=True)
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--samples", required=True, help="BFD samples.parquet (ASMID, BIOPROJECT)")
    ap.add_argument("--bioproject", action="append", default=[],
                    help="BioProject accession to treat as already-public; repeatable")
    ap.add_argument("--prefix", default="pfam_")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="24GB")
    ap.add_argument("--label", default="run")
    ap.add_argument("--out", required=True, help="per-family TSV")
    a = ap.parse_args(argv)

    import duckdb
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={a.threads}")
    con.execute(f"PRAGMA memory_limit='{a.memory}'")

    asmids = y1000_asmids(con, a.samples, a.bioproject)
    print(f"{len(asmids)} genomes match {a.bioproject}", file=sys.stderr)
    con.execute("create table y (ASMID VARCHAR)")
    if asmids:
        con.executemany("insert into y values (?)", [(x,) for x in sorted(asmids)])

    fams = family_names(a.null_dir, a.prefix)
    print(f"{len(fams)} families from {a.null_dir}", file=sys.stderr)

    rows = []
    for i, fam in enumerate(fams, 1):
        got = family_rows(con, fam, a.regions_dir, a.clusters_dir)
        if not got:
            print(f"  skip {fam}: no regions/clusters file", file=sys.stderr)
        rows.extend(got)
        if i % 25 == 0:
            print(f"  {i}/{len(fams)}", file=sys.stderr)

    cols = ["family", "tier", "clusters", "novel_now", "novel_if_y_known", "y_only",
            "novel_now_pct", "novel_if_y_known_pct"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\n=== {a.label} ===")
    for tier in TIERS:
        s = summarize(rows, tier)
        if not s["families"]:
            continue
        print(f"{tier}: {s['families']} families, {s['clusters']:,} clusters | "
              f"novel now {s['novel_now']:,} ({s['novel_now_pct']:.1f}%) | "
              f"if Y1000+ known {s['novel_if_y_known']:,} ({s['novel_if_y_known_pct']:.1f}%) | "
              f"lost {s['y_only']:,} ({s['pct_of_novel_lost']:.1f}% of novel)")
    print(f"{len(rows)} rows -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
