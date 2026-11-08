#!/usr/bin/env python3.12
"""
G3 Phase A step 1a: compute promoter intervals for target UPO genes and a matched background.

Promoter = --upstream bp immediately 5' of the CDS start on the gene's own strand, truncated at
the contig edge and at the nearest neighbouring gene so no promoter contains coding sequence.
Genes whose usable promoter is shorter than --min-len are dropped and counted, never padded.

Background = --bg-per-genome non-UPO protein-coding genes per genome, drawn with a fixed seed
from the SAME genomes, and length-matched: each background promoter is truncated to the length
of a randomly paired target promoter from that genome. Matching to the same genomes controls GC
content and gene-prediction style, which differ enormously across 12,000 fungal assemblies.

Runs in LOCUSTAG batches so the window function never touches the whole 245M-row table.
"""
import argparse, os, random, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--targets", required=True, help="file of target gene_id, one per line")
ap.add_argument("--gene-info", required=True)
ap.add_argument("--upstream", type=int, default=1000)
ap.add_argument("--min-len", type=int, default=300)
ap.add_argument("--bg-per-genome", type=int, default=5)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--batch", type=int, default=800, help="LOCUSTAGs per duckdb query")
ap.add_argument("--out", required=True, help="output TSV of promoter intervals")
ap.add_argument("--memory", default="24GB")
a = ap.parse_args()

targets = [l.strip() for l in open(a.targets) if l.strip()]
tset = set(targets)
lts = sorted({t[:8] for t in targets})
print(f"{len(targets)} target genes in {len(lts)} genomes", file=sys.stderr)

rng = random.Random(a.seed)
out = open(a.out, "w")
out.write("gene_id\tLOCUSTAG\tchrom\tprom_start\tprom_end\tstrand\tset\tlength\n")
n_t = n_b = 0
dropped_short = 0

for bi in range(0, len(lts), a.batch):
    chunk = lts[bi:bi + a.batch]
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("\n".join(chunk) + "\n"); ltfile = fh.name
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        res = fh.name
    q = f"""
    set memory_limit='{a.memory}';
    create table lt as select column0 as LOCUSTAG from read_csv('{ltfile}', header=false);
    create table g as select gene_id, LOCUSTAG, chrom, start, "end", strand
      from '{a.gene_info}'
      where gene_type='protein_coding' and LOCUSTAG in (select LOCUSTAG from lt);
    create table w as select *,
      lag("end") over (partition by LOCUSTAG, chrom order by start) prev_end,
      lead(start) over (partition by LOCUSTAG, chrom order by start) next_start
      from g;
    copy (select gene_id, LOCUSTAG, chrom, start, "end", strand, prev_end, next_start from w)
      to '{res}' (header, delimiter '\t');
    """
    r = subprocess.run(["duckdb", "-c", q], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"duckdb failed on batch {bi}: {r.stderr[:500]}")

    per_genome_targets = {}
    bg_pool = {}
    with open(res) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid, lt_, chrom = f[ix["gene_id"]], f[ix["LOCUSTAG"]], f[ix["chrom"]]
            st, en = int(f[ix["start"]]), int(f[ix["end"]])
            strand = int(f[ix["strand"]])
            pe = f[ix["prev_end"]]; ns = f[ix["next_start"]]
            pe = int(pe) if pe not in ("", "NULL") else None
            ns = int(ns) if ns not in ("", "NULL") else None
            if strand >= 0:
                lo = max(1, st - a.upstream, (pe + 1) if pe else 1)
                hi = st - 1
            else:
                lo = en + 1
                hi = en + a.upstream
                if ns:
                    hi = min(hi, ns - 1)
            length = hi - lo + 1
            rec = (gid, lt_, chrom, lo, hi, strand, length)
            if gid in tset:
                per_genome_targets.setdefault(lt_, []).append(rec)
            else:
                bg_pool.setdefault(lt_, []).append(rec)

    for lt_, recs in per_genome_targets.items():
        kept = [r_ for r_ in recs if r_[6] >= a.min_len]
        dropped_short += len(recs) - len(kept)
        for r_ in kept:
            out.write(f"{r_[0]}\t{r_[1]}\t{r_[2]}\t{r_[3]}\t{r_[4]}\t{r_[5]}\ttarget\t{r_[6]}\n")
            n_t += 1
        pool = [r_ for r_ in bg_pool.get(lt_, []) if r_[6] >= a.min_len]
        if not pool or not kept:
            continue
        rng.shuffle(pool)
        for r_ in pool[:a.bg_per_genome]:
            # length-match: truncate the background promoter to a paired target's length
            tgt_len = rng.choice(kept)[6]
            L = min(r_[6], tgt_len)
            if r_[5] >= 0:
                lo, hi = r_[4] - L + 1, r_[4]
            else:
                lo, hi = r_[3], r_[3] + L - 1
            out.write(f"{r_[0]}\t{r_[1]}\t{r_[2]}\t{lo}\t{hi}\t{r_[5]}\tbackground\t{L}\n")
            n_b += 1
    os.unlink(ltfile); os.unlink(res)
    print(f"  batch {bi//a.batch+1}/{(len(lts)+a.batch-1)//a.batch}: "
          f"{n_t} targets, {n_b} background", file=sys.stderr)

out.close()
print(f"\ntargets written: {n_t}", file=sys.stderr)
print(f"background written: {n_b}", file=sys.stderr)
print(f"targets dropped (< {a.min_len} bp usable): {dropped_short}", file=sys.stderr)
