#!/usr/bin/env python3.12
"""
G3 Phase A, step 1: pull promoter sequences for UPO genes, plus a matched background.

Scope is set by the gene list you pass, so this runs unchanged under either answer to decision D6
(all-BFD, or only clades where RNA-seq follow-up is possible).

Promoter = --upstream bp immediately 5' of the CDS start, on the gene's own strand, truncated at
the contig edge and at the 3' end of the nearest upstream gene so the sequence does not run into
a neighbouring coding region. Genes whose usable promoter is shorter than --min-len are reported
and dropped rather than padded.

Background = promoters of --bg-per-genome randomly chosen non-UPO genes from the SAME genomes,
drawn with a fixed seed. Matching the background to the same genomes is what makes any motif
enrichment interpretable: it controls for genome GC content and for gene-prediction style, which
differ enormously across 15,000 fungal assemblies.
"""
import argparse, gzip, hashlib, os, random, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--genes", required=True,
                help="TSV with a protein_id or gene_id column (upo_sites output works)")
ap.add_argument("--id-col", default="protein_id")
ap.add_argument("--tables", required=True, help="BFD tables dir holding gene_info.parquet")
ap.add_argument("--genomes-dir", required=True,
                help=".../genome_annotation  (per-genome funannotate output)")
ap.add_argument("--seqdb-genomes", default=None,
                help="optional seqdb genomes.tsv mapping LOCUSTAG to output folder name")
ap.add_argument("--upstream", type=int, default=1000)
ap.add_argument("--min-len", type=int, default=200)
ap.add_argument("--bg-per-genome", type=int, default=5)
ap.add_argument("--seed", type=int, default=1337)
ap.add_argument("--out-fg", required=True)
ap.add_argument("--out-bg", required=True)
ap.add_argument("--report", required=True)
a = ap.parse_args()

sys.exit("NOT YET RUNNABLE: blocked on decision D6 (scope) and on confirming the genome FASTA "
         "layout for contig sequence retrieval. The argument contract above is settled; the "
         "body is not written. See PLAN.md section 4, G3.")
