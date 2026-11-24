# Model-aware versions of sim_edge_flags and null_replicate: each step is simulated and refitted under its own model
# (ER or ARD, chosen by AIC on the real data). Source AFTER gainloss_null_lib.R; the definitions below replace the ARD-only ones.
# Needed for a near-universal trait: its ARD gain rate is about 0 and an ARD simulation could not gain, while ER fits as well.

sim_edge_flags <- function(tree, states, model = "ARD") {
  ne <- nrow(tree$edge)
  none <- list(gain = rep(FALSE, ne), loss = rep(FALSE, ne), degenerate = TRUE)
  if (length(unique(states)) < 2) return(none)
  f <- tryCatch(asr_binary(tree, states, model, Ntrials = 3), error = function(e) NULL)
  if (is.null(f)) return(none)
  ev <- call_events(tree, f$p1)
  list(gain = ev$class == "gain", loss = ev$class == "loss", degenerate = FALSE)
}

# params: step, model, r01, r10, root_p1, n_present
null_replicate <- function(tree, params, root = "stationary", match_tol = 0) {
  flags <- list(); stp <- list(); internal <- tree$edge[, 2] > Ntip(tree)
  for (k in seq_len(nrow(params))) {
    rt <- if (identical(root, "asr")) c(1 - params$root_p1[k], params$root_p1[k]) else if (identical(root, "flat")) c(0.5, 0.5) else "stationary"
    draws <- 1
    if (match_tol > 0) { m <- simulate_binary_matched(tree, params$r01[k], params$r10[k], params$n_present[k], match_tol, rt); st <- m$states; draws <- m$draws
    } else st <- simulate_binary(tree, params$r01[k], params$r10[k], 1, rt)[[1]]
    fl <- sim_edge_flags(tree, st, params$model[k]); flags[[params$step[k]]] <- fl
    stp[[k]] <- data.frame(step = params$step[k], n_present = sum(st), draws = draws, gains = sum(fl$gain), losses = sum(fl$loss), degenerate = fl$degenerate,
                           internal_gains = sum(fl$gain & internal), internal_losses = sum(fl$loss & internal))
  }
  list(counts = co_count_table(tree, flags), steps = do.call(rbind, stp))
}

# Fit ER and ARD to the real tip states of every step; return params (best model by AIC) for the simulation.
fit_params_by_aic <- function(tree, calls, steps) {
  do.call(rbind, lapply(steps, function(s) {
    v <- setNames(suppressWarnings(as.numeric(calls[[s]])), calls$locustag)[tree$tip.label]
    fe <- asr_binary(tree, v, "ER", Ntrials = 5); fa <- asr_binary(tree, v, "ARD", Ntrials = 5)
    best <- if (fa$AIC < fe$AIC) fa else fe; m <- if (fa$AIC < fe$AIC) "ARD" else "ER"
    data.frame(step = s, model = m, r01 = best$rate_01, r10 = best$rate_10, root_p1 = best$p1[Ntip(tree) + 1], n_present = sum(v), AIC_ER = fe$AIC, AIC_ARD = fa$AIC)
  }))
}
