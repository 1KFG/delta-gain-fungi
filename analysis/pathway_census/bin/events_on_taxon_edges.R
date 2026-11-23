#!/usr/bin/env Rscript
# Count the events of a run that sit on edges whose tips all belong to one phylum (default Microsporidia), or only partly.
# Usage: Rscript bin/events_on_taxon_edges.R <run dir with primary.rds, edge_events_primary.tsv, calls.tsv> <phylum> <out tsv>
suppressPackageStartupMessages(library(ape))
a <- commandArgs(TRUE); prim <- readRDS(file.path(a[1], "primary.rds")); tree <- prim$common_tree
ee <- read.delim(file.path(a[1], "edge_events_primary.tsv"), stringsAsFactors = FALSE)
cal <- read.delim(file.path(a[1], "calls.tsv"), stringsAsFactors = FALSE); ph <- setNames(cal$PHYLUM, cal$locustag)
pp <- prop.part(tree, check.labels = TRUE); tl <- attr(pp, "labels"); n <- Ntip(tree)
below <- function(nd) if (nd <= n) tree$tip.label[nd] else tl[pp[[nd - n]]]
frac <- vapply(ee$child, function(nd) mean(ph[below(nd)] == a[2]), 0)
cat("tips of", a[2], "in the common tree:", sum(ph[tree$tip.label] == a[2]), "\n")
rows <- list()
for (s in prim$steps) for (ev in c("gain", "loss")) for (et in c("tip", "internal")) {
  sel <- ee[[paste0(s, "_event")]] == ev & (ee$n_tips == 1) == (et == "tip")
  rows[[length(rows) + 1]] <- data.frame(step = s, event = ev, edge_type = et, n_events = sum(sel), all_tips_in_phylum = sum(sel & frac == 1),
                                         some_tips_in_phylum = sum(sel & frac > 0 & frac < 1))
}
res <- do.call(rbind, rows); write.table(res, a[3], sep = "\t", quote = FALSE, row.names = FALSE); print(res[res$all_tips_in_phylum > 0 | res$some_tips_in_phylum > 0, ])
