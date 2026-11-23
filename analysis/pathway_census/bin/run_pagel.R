#!/usr/bin/env Rscript
# Pagel (1994) test of correlated evolution of two binary traits with phytools::fitPagel (fitMk back end).
# Usage: Rscript bin/run_pagel.R <primary.rds> <step_a> <step_b> <out.tsv> [max_tips]
# Uses the common tree of the primary run (tips with data for all steps). If max_tips is given, a random subset of
# tips of that size is kept (seed 1), to measure run time on a smaller tree.
suppressPackageStartupMessages({library(ape); library(phytools)})
`%||%` <- function(a, b) if (is.null(a)) b else a
a <- commandArgs(TRUE)
prim <- readRDS(a[1]); sa <- a[2]; sb <- a[3]; out <- a[4]; maxn <- if (length(a) >= 5) as.integer(a[5]) else NA
calls <- read.delim(file.path(dirname(a[1]), "calls.tsv"), stringsAsFactors = FALSE)
tr <- prim$common_tree
if (!is.na(maxn) && maxn < Ntip(tr)) { set.seed(1); tr <- keep.tip(tr, sample(tr$tip.label, maxn)) }
tr$edge.length[tr$edge.length < 1e-8] <- 1e-8
x <- factor(setNames(calls[[sa]], calls$locustag)[tr$tip.label], levels = c(0, 1))
y <- factor(setNames(calls[[sb]], calls$locustag)[tr$tip.label], levels = c(0, 1))
t0 <- proc.time()[["elapsed"]]
fit <- tryCatch(fitPagel(tr, x, y), error = function(e) e)
secs <- proc.time()[["elapsed"]] - t0
if (inherits(fit, "error")) {
  res <- data.frame(step_a = sa, step_b = sb, n_tips = Ntip(tr), ok = FALSE, msg = conditionMessage(fit), secs = secs)
} else {
  res <- data.frame(step_a = sa, step_b = sb, n_tips = Ntip(tr), ok = TRUE, msg = "",
    logL_indep = as.numeric(fit$independent.logL), logL_dep = as.numeric(fit$dependent.logL),
    lik_ratio = as.numeric(fit$lik.ratio), df = attr(fit$dependent.logL, "df") - attr(fit$independent.logL, "df"), p = as.numeric(fit$P),
    AIC_indep = as.numeric(fit$independent.AIC), AIC_dep = as.numeric(fit$dependent.AIC), secs = secs)
}
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
