#!/usr/bin/env Rscript
# Summarise the simulation null: empirical p, null mean and 95% interval next to the observed count, per variant.
# Usage: Rscript bin/summarise_null.R <primary dir> <null dir>
# Observed counts come from the same code as the replicates (ARD refit on the common tree, 0.8 rule).
suppressPackageStartupMessages({library(ape)})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "gainloss_lib.R")); source(file.path(here, "gainloss_null_lib.R"))
a <- commandArgs(TRUE); pd <- a[1]; nd <- a[2]
prim <- readRDS(file.path(pd, "primary.rds")); tree <- prim$common_tree; steps <- prim$steps
cal <- read.delim(file.path(pd, "calls.tsv"), stringsAsFactors = FALSE)
internal <- tree$edge[, 2] > Ntip(tree)
flags <- list(); obs_steps <- list()
for (s in steps) {
  v <- setNames(suppressWarnings(as.numeric(cal[[s]])), cal$locustag)[tree$tip.label]
  f <- asr_binary(tree, v, "ARD", Ntrials = 5); ev <- call_events(tree, f$p1)
  flags[[s]] <- list(gain = ev$class == "gain", loss = ev$class == "loss")
  obs_steps[[s]] <- data.frame(step = s, n_present = sum(v), gains = sum(flags[[s]]$gain), losses = sum(flags[[s]]$loss),
                               internal_gains = sum(flags[[s]]$gain & internal), internal_losses = sum(flags[[s]]$loss & internal))
}
obs <- co_count_table(tree, flags); obs_steps <- do.call(rbind, obs_steps)
write.table(obs, file.path(nd, "observed_co_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
files <- list.files(nd, pattern = "^sim_counts\\..*\\.task[0-9]+\\.tsv$", full.names = TRUE)
variants <- unique(sub("^sim_counts\\.(.*)\\.task[0-9]+\\.tsv$", "\\1", basename(files)))
all_sum <- list(); all_step <- list()
for (vr in variants) {
  cs <- do.call(rbind, lapply(list.files(nd, pattern = paste0("^sim_counts\\.", vr, "\\.task[0-9]+\\.tsv$"), full.names = TRUE), read.delim, stringsAsFactors = FALSE))
  ss <- do.call(rbind, lapply(list.files(nd, pattern = paste0("^sim_steps\\.", vr, "\\.task[0-9]+\\.tsv$"), full.names = TRUE), read.delim, stringsAsFactors = FALSE))
  N <- length(unique(cs$rep)); rows <- list()
  for (i in seq_len(nrow(obs))) for (sc in c("n_all", "n_internal")) {
    sub <- cs[cs$step_a == obs$step_a[i] & cs$step_b == obs$step_b[i] & cs$type == obs$type[i], ]
    nv <- sub[[sc]]; ov <- obs[[sc]][i]
    rows[[length(rows) + 1]] <- data.frame(variant = vr, step_a = obs$step_a[i], step_b = obs$step_b[i], type = obs$type[i],
      scope = ifelse(sc == "n_all", "all_edges", "internal_edges"), observed = ov, null_mean = mean(nv), null_lo = unname(quantile(nv, 0.025)),
      null_hi = unname(quantile(nv, 0.975)), n_rep = length(nv), n_ge = sum(nv >= ov), p = empirical_p(ov, nv))
  }
  sm <- do.call(rbind, rows); sm$q <- NA_real_
  for (sc in unique(sm$scope)) { e <- sm$scope == sc; sm$q[e] <- p.adjust(sm$p[e], "BH") }
  all_sum[[vr]] <- sm
  dg <- do.call(rbind, lapply(steps, function(s) { x <- ss[ss$step == s, ]; o <- obs_steps[obs_steps$step == s, ]
    data.frame(variant = vr, step = s, n_rep = nrow(x), obs_present = o$n_present, null_present_mean = mean(x$n_present), null_present_lo = quantile(x$n_present, 0.025), null_present_hi = quantile(x$n_present, 0.975),
      frac_null_present_ge_obs = mean(x$n_present >= o$n_present), obs_gains = o$gains, null_gains_mean = mean(x$gains), null_gains_lo = quantile(x$gains, 0.025), null_gains_hi = quantile(x$gains, 0.975),
      obs_losses = o$losses, null_losses_mean = mean(x$losses), null_losses_lo = quantile(x$losses, 0.025), null_losses_hi = quantile(x$losses, 0.975),
      obs_internal_gains = o$internal_gains, null_internal_gains_mean = mean(x$internal_gains), obs_internal_losses = o$internal_losses, null_internal_losses_mean = mean(x$internal_losses),
      frac_degenerate = mean(x$degenerate), mean_draws_per_accept = if ("draws" %in% names(x)) mean(x$draws) else NA_real_) }))
  all_step[[vr]] <- dg
  cat(sprintf("variant %s: N=%d, smallest possible p = %.5f\n", vr, N, 1 / (N + 1)))
}
w <- function(x, f) write.table(x, file.path(nd, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(do.call(rbind, all_sum), "null_summary.tsv"); w(do.call(rbind, all_step), "null_step_diagnostics.tsv")
