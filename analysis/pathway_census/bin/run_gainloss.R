#!/usr/bin/env Rscript
# Ancestral state reconstruction, event counts, parsimony and co-event permutation tests for the carotenoid steps.
# Usage: Rscript bin/run_gainloss.R --calls calls.tsv --tree <fasttree treefile> --outdir <dir> [--nperm 10000] [--primary Opisthosporidia]
# Steps are the step columns of calls.tsv (columns `<step>_source` give hits or placed).
suppressPackageStartupMessages({library(ape)})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "gainloss_lib.R"))
a <- commandArgs(TRUE)
opt <- list(nperm = 10000, primary = "Opisthosporidia")
for (i in seq(1, length(a), 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
opt$nperm <- as.integer(opt$nperm)
cal <- read.delim(opt$calls, stringsAsFactors = FALSE, check.names = FALSE)
src_cols <- grep("_source$", names(cal), value = TRUE); steps <- sub("_source$", "", src_cols)
source_of <- setNames(vapply(src_cols, function(c) cal[[c]][1], ""), steps)
tree0 <- read.tree(opt$tree)
ph <- setNames(cal$PHYLUM, cal$locustag); cls <- setNames(cal$CLASS, cal$locustag)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
st_of <- function(s) { v <- suppressWarnings(as.numeric(cal[[s]])); names(v) <- cal$locustag; v }

# --- roots (PROPOSED): each is a set of tips that forms one side of a split of the unrooted tree
groups <- list(Opisthosporidia = c("Rozellomycota", "Aphelidiomycota", "Microsporidia"),
               Blastocladiomycota = "Blastocladiomycota", Chytridiomycota = "Chytridiomycota")
root_tips <- lapply(groups, function(g) largest_pure_split(tree0, tree0$tip.label[ph[tree0$tip.label] %in% g]))
roots <- c(opt$primary, setdiff(names(groups), opt$primary))
cat("root tips:", paste(names(root_tips), lengths(root_tips), collapse = "; "), "\n")

# common tip set: tips with data for every step
common <- Reduce(intersect, lapply(steps, function(s) names(which(!is.na(st_of(s))))))
common <- intersect(common, tree0$tip.label)
cat("tips in tree:", Ntip(tree0), " tips with data for all steps:", length(common), "\n")

counts <- list(); thr <- list(); primary <- NULL; pairs_all <- list(); edge_tabs <- list()
for (r in roots) {
  rt <- root_on_split(tree0, root_tips[[r]])
  fits <- list(); evs <- list(); trees_s <- list(); best <- list()
  for (s in steps) {
    pr <- prune_to_data(rt, st_of(s)); trees_s[[s]] <- pr$tree
    par <- parsimony_counts(pr$tree, pr$states)
    for (m in c("ER", "ARD")) {
      t0 <- proc.time()[["elapsed"]]
      f <- asr_binary(pr$tree, pr$states, m, Ntrials = 5)
      ev <- call_events(pr$tree, f$p1); sm <- summarise_events(ev)
      fits[[paste(s, m)]] <- f; evs[[paste(s, m)]] <- ev
      counts[[length(counts) + 1]] <- data.frame(root = r, step = s, source = source_of[[s]], model = m, n_tips = Ntip(pr$tree),
        n_present = sum(pr$states == 1), n_absent = sum(pr$states == 0), logLik = f$loglik, npar = f$npar, AIC = f$AIC,
        rate_01 = f$rate_01, rate_10 = f$rate_10, gains = sm$gains, losses = sm$losses, uncertain = sm$uncertain, same = sm$same,
        n_edges = sm$n_edges, pars_changes = par$changes, pars_changes_castor = par$changes_castor,
        pars_gains_late = par$gains_late, pars_losses_late = par$losses_late,
        pars_gains_early = par$gains_early, pars_losses_early = par$losses_early, secs = round(proc.time()[["elapsed"]] - t0, 1))
      if (r == roots[1]) {
        for (cut in c(0.7, 0.8, 0.9)) { e2 <- summarise_events(call_events(pr$tree, f$p1, cut))
          thr[[length(thr) + 1]] <- data.frame(root = r, step = s, model = m, cut = cut, gains = e2$gains, losses = e2$losses, uncertain = e2$uncertain) }
      }
    }
    aic <- c(ER = fits[[paste(s, "ER")]]$AIC, ARD = fits[[paste(s, "ARD")]]$AIC); best[[s]] <- names(which.min(aic))
  }
  # --- pair tests on the common tree (same as per-step trees when the NA sets are identical)
  tc <- keep.tip(rt, common); evc <- list()
  for (s in steps) { v <- st_of(s)[common]; f <- asr_binary(tc, v, best[[s]], Ntrials = 5); evc[[s]] <- call_events(tc, f$p1) }
  tests <- list(); k <- 0
  for (i in seq_along(steps)) for (j in seq_along(steps)) if (i < j) {
    sa <- steps[i]; sb <- steps[j]
    for (ty in list(c("gain", "gain"), c("loss", "loss"), c("gain", "loss"), c("loss", "gain"))) {
      k <- k + 1
      a_ev <- evc[[sa]]$class == ty[1]; b_ev <- evc[[sb]]$class == ty[2]
      po <- perm_overlap(a_ev, b_ev, eligible_edges(evc[[sb]], ty[2]), nperm = opt$nperm, seed = 1000 + k)
      tests[[k]] <- data.frame(root = r, step_a = sa, step_b = sb, event_a = ty[1], event_b = ty[2], model_a = best[[sa]], model_b = best[[sb]],
        n_events_a = po$n_a, n_events_b = po$n_b, n_eligible_b = po$n_elig_b, observed = po$observed, expected = po$expected,
        perm_mean = po$perm_mean, p = po$p)
    }
  }
  tests <- do.call(rbind, tests); tests$q <- p.adjust(tests$p, method = "BH")
  pairs_all[[r]] <- tests
  if (r == roots[1]) {
    clade <- edge_clade_summary(tc, cls)
    et <- data.frame(edge = seq_len(nrow(tc$edge)), parent = tc$edge[, 1], child = tc$edge[, 2], length = tc$edge.length, clade)
    for (s in steps) et[[paste0(s, "_event")]] <- evc[[s]]$class
    edge_tabs <- et
    primary <- list(root = r, tree = rt, common_tree = tc, step_trees = trees_s, best_model = best, source = source_of,
                    p1 = lapply(steps, function(s) fits[[paste(s, best[[s]])]]$p1), events = lapply(steps, function(s) evs[[paste(s, best[[s]])]]), steps = steps)
    names(primary$p1) <- steps; names(primary$events) <- steps
  }
}
w <- function(x, f) write.table(x, file.path(opt$outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
cnt <- do.call(rbind, counts)
cnt$best_by_AIC <- ave(cnt$AIC, cnt$root, cnt$step, FUN = function(x) x == min(x))
cnt$delta_AIC_ARD_minus_ER <- NA_real_
for (r in unique(cnt$root)) for (s in steps) {
  e <- cnt$root == r & cnt$step == s
  cnt$delta_AIC_ARD_minus_ER[e] <- cnt$AIC[e & cnt$model == "ARD"] - cnt$AIC[e & cnt$model == "ER"]
}
w(cnt, "gainloss_counts.tsv")
w(do.call(rbind, thr), "threshold_sensitivity.tsv")
pa <- do.call(rbind, pairs_all); w(pa, "co_events.tsv")
w(edge_tabs, "edge_events_primary.tsv")
primary$best_fits <- NULL
saveRDS(primary, file.path(opt$outdir, "primary.rds"))
write.tree(primary$tree, file.path(opt$outdir, "tree.rooted.primary.nwk"))
cat("done\n")
