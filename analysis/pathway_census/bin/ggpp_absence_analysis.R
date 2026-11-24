#!/usr/bin/env Rscript
# Where is a near-universal step absent, and how are the absences clustered on the tree? Also node-call counts for every step.
# Usage: Rscript bin/ggpp_absence_analysis.R <run dir with primary.rds, calls.tsv, edge_events_primary.tsv> <step> <outdir>
suppressPackageStartupMessages(library(ape))
a <- commandArgs(TRUE); rd <- a[1]; step <- a[2]; od <- a[3]
prim <- readRDS(file.path(rd, "primary.rds")); cal <- read.delim(file.path(rd, "calls.tsv"), stringsAsFactors = FALSE)
tree <- prim$step_trees[[step]]; n <- Ntip(tree)
v <- setNames(suppressWarnings(as.numeric(cal[[step]])), cal$locustag)[tree$tip.label]
cls <- setNames(ifelse(cal$CLASS == "" | is.na(cal$CLASS), "(no class)", cal$CLASS), cal$locustag)[tree$tip.label]
# absent tips by class
tb <- data.frame(class = names(table(cls)), n_tips = as.integer(table(cls)), stringsAsFactors = FALSE)
tb$n_absent <- sapply(tb$class, function(k) sum(v[cls == k] == 0)); tb$frac_absent <- round(tb$n_absent / tb$n_tips, 3)
tb <- tb[order(-tb$n_absent, -tb$n_tips), ]; write.table(tb, file.path(od, paste0(step, "_absent_by_class.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
# maximal all-absent clades (blocks)
pp <- prop.part(tree, check.labels = TRUE); tl <- attr(pp, "labels")
node_all_absent <- vapply(seq_along(pp), function(i) all(v[tl[pp[[i]]]] == 0), TRUE)
cand <- which(node_all_absent)
covered <- rep(FALSE, n); names(covered) <- tree$tip.label; blocks <- list()
# take larger clades first
ord <- cand[order(-lengths(pp[cand]))]
for (i in ord) { tips <- tl[pp[[i]]]; if (!any(covered[tips])) { covered[tips] <- TRUE; blocks[[length(blocks) + 1]] <- tips } }
single <- names(v)[v == 0 & !covered]
rows <- lapply(blocks, function(t) { tb2 <- sort(table(cls[t]), decreasing = TRUE); data.frame(block_size = length(t), top_class = names(tb2)[1], n_top_class = as.integer(tb2[1]), classes = paste(names(tb2), tb2, sep = ":", collapse = "; ")) })
rows <- c(rows, lapply(single, function(t) data.frame(block_size = 1, top_class = cls[t], n_top_class = 1, classes = paste0(cls[t], ":1"))))
bl <- do.call(rbind, rows); bl <- bl[order(-bl$block_size), ]
write.table(bl, file.path(od, paste0(step, "_absent_blocks.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("absent tips:", sum(v == 0), " in blocks of size >= 2:", sum(covered), " singletons:", length(single), " blocks (size>=2):", length(blocks), "\n")
print(head(bl, 15), row.names = FALSE); print(tb[tb$n_absent > 0, ], row.names = FALSE)
# node calls for all steps (internal nodes only)
rows <- lapply(prim$steps, function(s) { p <- prim$p1[[s]]; tr <- prim$step_trees[[s]]; ni <- Ntip(tr)
  pn <- p[(ni + 1):length(p)]; data.frame(step = s, internal_nodes = length(pn), called_present = sum(pn >= 0.8), called_absent = sum(pn <= 0.2), uncertain = sum(pn > 0.2 & pn < 0.8)) })
nc <- do.call(rbind, rows); write.table(nc, file.path(od, "internal_node_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(nc, row.names = FALSE)
# loss events of this step by clade class and edge type
ee <- read.delim(file.path(rd, "edge_events_primary.tsv"), stringsAsFactors = FALSE); ev <- ee[[paste0(step, "_event")]]
lo <- ee[ev == "loss", ]; lo$edge_type <- ifelse(lo$n_tips == 1, "tip", "internal")
write.table(lo[, c("edge", "n_tips", "top_label", "top_frac", "edge_type")], file.path(od, paste0(step, "_loss_edges.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(lo$edge_type)); print(sort(table(lo$top_label[lo$edge_type == "internal"]), decreasing = TRUE))
