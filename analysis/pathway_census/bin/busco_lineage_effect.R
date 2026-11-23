#!/usr/bin/env Rscript
# How the BUSCO filter changes events by lineage. Usage: Rscript bin/busco_lineage_effect.R <unfiltered dir> <busco dir> <reps tsv> <min busco> <outdir>
# For each event edge of the unfiltered run (ARD, 0.8 rule), report whether all, some or none of the tips below it are
# removed by the filter, by step, event type and edge type (tip or internal).
suppressPackageStartupMessages({library(ape)})
a <- commandArgs(TRUE)
u <- readRDS(file.path(a[1], "primary.rds")); tree <- u$common_tree
reps <- read.delim(a[3], stringsAsFactors = FALSE); removed <- reps$LOCUSTAG[reps$busco_complete < as.numeric(a[4])]
ee <- read.delim(file.path(a[1], "edge_events_primary.tsv"), stringsAsFactors = FALSE)
pp <- prop.part(tree, check.labels = TRUE); tl <- attr(pp, "labels"); n <- Ntip(tree)
below <- function(node) if (node <= n) tree$tip.label[node] else tl[pp[[node - n]]]
ee$n_removed_below <- vapply(ee$child, function(nd) sum(below(nd) %in% removed), 0)
ee$frac_removed_below <- ee$n_removed_below / ee$n_tips
ee$edge_type <- ifelse(ee$n_tips == 1, "tip", "internal")
rows <- list()
for (s in u$steps) for (ev in c("gain", "loss")) {
  x <- ee[ee[[paste0(s, "_event")]] == ev, ]
  for (et in c("tip", "internal")) { y <- x[x$edge_type == et, ]
    rows[[length(rows) + 1]] <- data.frame(step = s, event = ev, edge_type = et, n_events = nrow(y), all_tips_removed = sum(y$frac_removed_below == 1),
      some_tips_removed = sum(y$frac_removed_below > 0 & y$frac_removed_below < 1), no_tips_removed = sum(y$frac_removed_below == 0)) }
}
write.table(do.call(rbind, rows), file.path(a[5], "busco_events_by_removal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# removed tips per class among tree tips
cal <- read.delim(file.path(a[1], "calls.tsv"), stringsAsFactors = FALSE)
cl <- setNames(ifelse(cal$CLASS == "" | is.na(cal$CLASS), "(no class)", cal$CLASS), cal$locustag)[tree$tip.label]
tb <- data.frame(class = names(table(cl)), n_tips_in_tree = as.integer(table(cl)), stringsAsFactors = FALSE)
tb$n_removed <- sapply(tb$class, function(k) sum(tree$tip.label[cl == k] %in% removed))
tb$fraction_removed <- round(tb$n_removed / tb$n_tips_in_tree, 3)
# presence among removed vs kept tips per step (calls of the unfiltered run)
for (s in u$steps) { v <- suppressWarnings(as.numeric(cal[[s]])); names(v) <- cal$locustag; v <- v[tree$tip.label]
  tb[[paste0(s, "_present_frac_removed")]] <- sapply(tb$class, function(k) { t <- tree$tip.label[cl == k & tree$tip.label %in% removed]; if (length(t)) round(mean(v[t]), 3) else NA }) }
tb <- tb[order(-tb$n_tips_in_tree), ]
write.table(tb, file.path(a[5], "busco_removed_by_class_in_tree.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
v <- tb[tb$n_tips_in_tree >= 8, c("class", "n_tips_in_tree", "n_removed", "fraction_removed")]; print(v, row.names = FALSE)
cat("tips in tree:", n, " removed among them:", sum(tree$tip.label %in% removed), "\n")
