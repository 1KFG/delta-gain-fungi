#!/usr/bin/env Rscript
# Compare event counts on the ASTRAL consensus tree and on the FastTree primary tree, same tips, same rules.
# Usage: Rscript bin/run_astral_compare.R <calls.tsv> <primary dir (tree.rooted.primary.nwk)> <outdir>
# Variants: fasttree_primary (substitutions per site), astral_asis (ASTRAL lengths, coalescent units for internal branches),
# astral_unit (every branch length set to 1: topology only). The ASTRAL tree keeps its own root (a single Blastocladiomycota tip).
suppressPackageStartupMessages({library(ape)})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "gainloss_lib.R")); source(file.path(here, "gainloss_null_lib.R"))
a <- commandArgs(TRUE); cal <- read.delim(a[1], stringsAsFactors = FALSE); pdir <- a[2]; od <- a[3]
P <- "/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/"
src_cols <- grep("_source$", names(cal), value = TRUE); steps <- sub("_source$", "", src_cols)
st_of <- function(s) { v <- suppressWarnings(as.numeric(cal[[s]])); names(v) <- cal$locustag; v }
common <- Reduce(intersect, lapply(steps, function(s) names(which(!is.na(st_of(s))))))
astral <- read.tree(paste0(P, "results/protein/tree/fungi_odb10/final_tree.nw"))
ft <- read.tree(file.path(pdir, "tree.rooted.primary.nwk"))
trees <- list(fasttree_primary = ft, astral_asis = astral, astral_unit = astral)
rows <- list(); co <- list(); lens <- list()
for (tn in names(trees)) {
  t <- trees[[tn]]; common_t <- intersect(common, t$tip.label)
  tc <- keep.tip(t, common_t)
  if (tn == "astral_unit") tc$edge.length <- rep(1, nrow(tc$edge))
  tc$edge.length[tc$edge.length < 1e-8] <- 1e-8
  internal <- tc$edge[, 2] > Ntip(tc); evs <- list(); best <- list()
  for (s in steps) {
    v <- st_of(s)[tc$tip.label]
    fits <- lapply(c(ER = "ER", ARD = "ARD"), function(m) asr_binary(tc, v, m, Ntrials = 5))
    best[[s]] <- if (fits$ARD$AIC < fits$ER$AIC) "ARD" else "ER"
    for (m in names(fits)) {
      ev <- call_events(tc, fits[[m]]$p1); if (m == "ARD") evs[[s]] <- list(gain = ev$class == "gain", loss = ev$class == "loss")
      rows[[length(rows) + 1]] <- data.frame(tree = tn, step = s, model = m, n_tips = Ntip(tc), logLik = fits[[m]]$loglik, AIC = fits[[m]]$AIC,
        rate_01 = fits[[m]]$rate_01, rate_10 = fits[[m]]$rate_10, gains = sum(ev$class == "gain"), losses = sum(ev$class == "loss"),
        uncertain = sum(ev$class == "uncertain"), gains_tip_edge = sum(ev$class == "gain" & !internal), gains_internal_edge = sum(ev$class == "gain" & internal),
        losses_tip_edge = sum(ev$class == "loss" & !internal), losses_internal_edge = sum(ev$class == "loss" & internal), best_by_AIC = (m == best[[s]]))
    }
  }
  ct <- co_count_table(tc, evs); ct$tree <- tn; co[[tn]] <- ct
}
w <- function(x, f) write.table(x, file.path(od, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(do.call(rbind, rows), "astral_vs_fasttree_counts.tsv"); w(do.call(rbind, co), "astral_vs_fasttree_co_events.tsv")
cat("done; tips used:", length(common), "\n")
