#!/usr/bin/env Rscript
# Compare the candidate genus-level trees and write tree_comparison.tsv and tree_rootsets.tsv.
# Usage: Rscript bin/gainloss_compare_trees.R <calls.tsv> <outdir>
suppressPackageStartupMessages({library(ape); library(phangorn)})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "gainloss_lib.R"))
a <- commandArgs(TRUE); calls_f <- a[1]; outdir <- a[2]
P <- "/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/"
paths <- c(astral = paste0(P, "results/protein/tree/fungi_odb10/final_tree.nw"),
  fasttree_support = paste0(P, "results/protein/buildtree/fungi_odb10/fasttree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fasttree.support.treefile"),
  fasttree_nosupport = paste0(P, "results/protein/buildtree/fungi_odb10/fasttree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fasttree.nosupport.treefile"))
reps <- read.delim(paste0(P, "genus_representatives.tsv"), stringsAsFactors = FALSE)
cal <- read.delim(calls_f, stringsAsFactors = FALSE)
trees <- lapply(paths, read.tree)
rows <- list()
for (n in names(trees)) {
  t <- trees[[n]]; el <- t$edge.length
  sup <- suppressWarnings(as.numeric(t$node.label)); sup <- sup[!is.na(sup)]
  rows[[n]] <- data.frame(tree = n, tips = Ntip(t), internal_nodes = t$Nnode, rooted = is.rooted(t), binary = is.binary(t),
    tips_in_reps = sum(t$tip.label %in% reps$LOCUSTAG), reps_missing_from_tree = sum(!reps$LOCUSTAG %in% t$tip.label),
    duplicate_tips = sum(duplicated(t$tip.label)), edge_length_min = min(el), edge_length_median = median(el), edge_length_max = max(el),
    edges_zero = sum(el == 0), edges_lt_1e6 = sum(el < 1e-6), edges_NA = sum(is.na(el)),
    root_to_tip_min = min(node.depth.edgelength(t)[seq_len(Ntip(t))]), root_to_tip_max = max(node.depth.edgelength(t)[seq_len(Ntip(t))]),
    n_support_values = length(sup), support_ge_0.8 = if (length(sup)) mean(sup >= 0.8) else NA_real_,
    tips_with_data_all_steps = sum(t$tip.label %in% cal$locustag[cal$searched == 1]))
}
cmp <- do.call(rbind, rows)
write.table(cmp, file.path(outdir, "tree_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(t(cmp))
# Topology distance (unrooted Robinson-Foulds) between trees on common tips
common <- Reduce(intersect, lapply(trees, `[[`, "tip.label"))
u <- lapply(trees, function(t) unroot(keep.tip(t, common)))
rf <- function(x, y) phangorn::RF.dist(u[[x]], u[[y]])
cat("common tips:", length(common), " max RF:", 2 * (length(common) - 3), "\n")
cat("RF astral vs fasttree_support:", rf("astral", "fasttree_support"), "\n")
cat("RF fasttree support vs nosupport:", rf("fasttree_support", "fasttree_nosupport"), "\n")
# Edge-length identity of the two FastTree files
cat("support/nosupport edge lengths identical:", isTRUE(all.equal(trees$fasttree_support$edge.length, trees$fasttree_nosupport$edge.length)), "\n")
# Mk log-likelihood (ER, flat root prior) of each step on each tree, same tip set (tips with data)
steps <- c("phytoene_synthase", "phytoene_desaturase", "lycopene_cyclase", "cart")
ll <- list()
for (n in c("astral", "fasttree_support")) for (s in steps) {
  v <- setNames(cal[[s]], cal$locustag); v <- suppressWarnings(as.numeric(v)); names(v) <- cal$locustag
  pr <- prune_to_data(trees[[n]], v); f <- asr_binary(pr$tree, pr$states, "ER", Ntrials = 3)
  ll[[length(ll) + 1]] <- data.frame(tree = n, step = s, n_tips = Ntip(pr$tree), logLik_ER = f$loglik, rate_ER = f$rate_01)
}
ll <- do.call(rbind, ll); print(ll)
write.table(ll, file.path(outdir, "tree_mk_loglik.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Root candidates on the FastTree support tree (tips of each group that form one side of a split of the unrooted tree)
ph <- setNames(cal$PHYLUM, cal$locustag)
ft <- trees$fasttree_support
groups <- list(Opisthosporidia = c("Rozellomycota", "Aphelidiomycota", "Microsporidia"),
               Blastocladiomycota = "Blastocladiomycota", Chytridiomycota = "Chytridiomycota")
rs <- lapply(names(groups), function(g) {
  tips <- largest_pure_split(ft, ft$tip.label[ph[ft$tip.label] %in% groups[[g]]])
  all_g <- ft$tip.label[ph[ft$tip.label] %in% groups[[g]]]
  data.frame(root = g, phyla = paste(groups[[g]], collapse = "+"), tips_of_phyla_in_tree = length(all_g), tips_in_split = length(tips),
             split_equals_all_phylum_tips = setequal(tips, all_g))
})
rs <- do.call(rbind, rs); print(rs)
write.table(rs, file.path(outdir, "tree_rootsets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
