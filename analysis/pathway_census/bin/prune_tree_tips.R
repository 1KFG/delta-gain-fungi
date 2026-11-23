#!/usr/bin/env Rscript
# Prune tips listed in a file from a tree. Usage: Rscript bin/prune_tree_tips.R <tree> <tips.txt> <out.nwk>
suppressPackageStartupMessages(library(ape))
a <- commandArgs(TRUE); t <- read.tree(a[1]); drop <- readLines(a[2]); drop <- drop[nzchar(drop)]
cat("tips before:", Ntip(t), " listed:", length(drop), " listed and in tree:", sum(drop %in% t$tip.label), "\n")
t2 <- drop.tip(t, intersect(drop, t$tip.label)); write.tree(t2, a[3]); cat("tips after:", Ntip(t2), "\n")
