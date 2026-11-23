#!/usr/bin/env Rscript
# Parametric-bootstrap null for co-events. One task = `nrep` replicates.
# Usage: Rscript bin/run_null_sim.R --primary <dir with primary.rds> --outdir <dir> --nrep 100 --task 1 [--root stationary|asr|flat] [--match_tol 0.1] [--seed 20261122]
# Steps are simulated independently of each other on the common tree of the primary run, with the ARD rates fitted
# to the real tip states of each step. Each replicate is refitted and reconstructed with the same code as the real data.
suppressPackageStartupMessages({library(ape)})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "gainloss_lib.R")); source(file.path(here, "gainloss_null_lib.R"))
a <- commandArgs(TRUE); opt <- list(root = "stationary", seed = 20261122, nrep = 100, task = 1, match_tol = 0)
for (i in seq(1, length(a), 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
opt$nrep <- as.integer(opt$nrep); opt$task <- as.integer(opt$task); opt$seed <- as.integer(opt$seed); opt$match_tol <- as.numeric(opt$match_tol)
prim <- readRDS(file.path(opt$primary, "primary.rds")); tree <- prim$common_tree
cal <- read.delim(file.path(opt$primary, "calls.tsv"), stringsAsFactors = FALSE)
params <- do.call(rbind, lapply(prim$steps, function(s) {
  v <- setNames(suppressWarnings(as.numeric(cal[[s]])), cal$locustag)[tree$tip.label]
  f <- asr_binary(tree, v, "ARD", Ntrials = 5)
  data.frame(step = s, r01 = f$rate_01, r10 = f$rate_10, root_p1 = f$p1[Ntip(tree) + 1], n_present = sum(v))
}))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
if (opt$task == 1) write.table(params, file.path(opt$outdir, "sim_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
set.seed(opt$seed + 1000 * opt$task)
t0 <- proc.time()[["elapsed"]]; cs <- list(); ss <- list()
for (r in seq_len(opt$nrep)) {
  x <- null_replicate(tree, params, opt$root, opt$match_tol)
  cs[[r]] <- cbind(rep = (opt$task - 1) * opt$nrep + r, x$counts); ss[[r]] <- cbind(rep = (opt$task - 1) * opt$nrep + r, x$steps)
}
secs <- proc.time()[["elapsed"]] - t0
tag <- sprintf("%s%s.task%03d", opt$root, if (opt$match_tol > 0) sprintf("_match%g", opt$match_tol) else "", opt$task)
write.table(do.call(rbind, cs), file.path(opt$outdir, paste0("sim_counts.", tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, ss), file.path(opt$outdir, paste0("sim_steps.", tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("task %d: %d replicates in %.1f s (%.2f s per replicate), root=%s match_tol=%g\n", opt$task, opt$nrep, secs, secs / opt$nrep, opt$root, opt$match_tol))
