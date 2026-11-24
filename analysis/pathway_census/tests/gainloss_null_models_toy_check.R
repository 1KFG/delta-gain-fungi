# Toy checks for bin/gainloss_null_models.R (model-aware simulation and fit). Prints key=value lines.
args <- commandArgs(FALSE); f <- sub("--file=", "", args[grep("--file=", args)]); bin <- file.path(dirname(normalizePath(f)), "..", "bin")
source(file.path(bin, "gainloss_lib.R")); source(file.path(bin, "gainloss_null_lib.R")); source(file.path(bin, "gainloss_null_models.R"))
set.seed(4); tr <- rcoal(120); tr$tip.label <- paste0("T", seq_len(120)); tr$edge.length <- tr$edge.length * 4
calls <- data.frame(locustag = tr$tip.label)
calls$X <- simulate_binary(tr, 0.4, 0.4, 1)[[1]][tr$tip.label]
calls$Y <- simulate_binary(tr, 0.03, 1.2, 1, root = c(1, 0))[[1]][tr$tip.label]
calls$Z <- simulate_binary(tr, 0.4, 0.4, 1)[[1]][tr$tip.label]
pr <- fit_params_by_aic(tr, calls, c("X", "Y", "Z"))
cat("params_rows=", nrow(pr), "\nparams_cols_ok=", as.integer(all(c("step", "model", "r01", "r10", "root_p1", "n_present", "AIC_ER", "AIC_ARD") %in% names(pr))), "\n", sep = "")
cat("model_matches_aic=", as.integer(all(pr$model == ifelse(pr$AIC_ARD < pr$AIC_ER, "ARD", "ER"))), "\n", sep = "")
cat("n_present_ok=", as.integer(all(pr$n_present == c(sum(calls$X), sum(calls$Y), sum(calls$Z)))), "\n", sep = "")
# ER fit has equal rates
er <- pr[pr$model == "ER", ]; cat("er_equal_rates=", as.integer(all(er$r01 == er$r10)), "\n", sep = "")
# model-aware flags: ER and ARD on the same data both return edge-length vectors; constant data degenerate
fe <- sim_edge_flags(tr, calls$X, "ER"); fa <- sim_edge_flags(tr, calls$X, "ARD")
cat("flags_len=", length(fe$gain), "\nflags_same_len=", as.integer(length(fe$gain) == length(fa$gain)), "\n", sep = "")
cat("const_degenerate=", as.integer(sim_edge_flags(tr, setNames(rep(1L, 120), tr$tip.label), "ER")$degenerate), "\n", sep = "")
# replicate: 3 steps -> 3 pairs x 4 types; the model column is used (forcing ER on all still returns rows)
pr$model <- "ER"; set.seed(2); r <- null_replicate(tr, pr); cat("rep_rows=", nrow(r$counts), "\nrep_steps=", nrow(r$steps), "\n", sep = "")
pr$model <- c("ER", "ARD", "ER"); set.seed(2); r2 <- null_replicate(tr, pr, match_tol = 0.2); cat("matched_rows=", nrow(r2$counts), "\n", sep = "")
