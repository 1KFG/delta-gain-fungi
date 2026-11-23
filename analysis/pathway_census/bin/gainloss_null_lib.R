# Tree-aware simulation null for co-events. Sourced by run_null_sim.R, summarise_null.R and the toy-tree test.
# Needs bin/gainloss_lib.R (asr_binary, call_events) to be sourced first.

SIM_PAIR_TYPES <- list(gain_gain = c("gain", "gain"), loss_loss = c("loss", "loss"), gain_loss = c("gain", "loss"), loss_gain = c("loss", "gain"))

# Q matrix of a 2-state model (state 1 = absent, 2 = present) from rates r01 (0->1) and r10 (1->0).
q_from_rates <- function(r01, r10) rbind(c(-r01, r01), c(r10, -r10))

# Stationary probabilities (P(absent), P(present)). If both rates are 0, use (0.5, 0.5).
stationary_probs <- function(r01, r10) { s <- r01 + r10; if (s == 0) c(0.5, 0.5) else c(r10 / s, r01 / s) }

# Simulate `n` tip-state vectors (0/1, named by tip) on `tree`. root: "stationary" or a numeric vector c(P(absent), P(present)).
simulate_binary <- function(tree, r01, r10, n = 1, root = "stationary") {
  rp <- if (identical(root, "stationary")) stationary_probs(r01, r10) else root
  s <- castor::simulate_mk_model(tree, q_from_rates(r01, r10), root_probabilities = rp, include_tips = TRUE,
                                 include_nodes = FALSE, Nsimulations = n, drop_dims = FALSE)$tip_states
  out <- lapply(seq_len(n), function(i) setNames(as.integer(s[i, ]) - 1L, tree$tip.label))
  out
}

# Simulate one tip-state vector whose number of present tips is within `tol` (fraction) of `n_obs`, by rejection sampling in
# batches. Returns the vector and the number of draws used. Stops with an error after `max_draws` draws.
simulate_binary_matched <- function(tree, r01, r10, n_obs, tol = 0.1, root = "stationary", batch = 2000, max_draws = 2e6) {
  rp <- if (identical(root, "stationary")) stationary_probs(r01, r10) else root
  Q <- q_from_rates(r01, r10); used <- 0
  while (used < max_draws) {
    s <- castor::simulate_mk_model(tree, Q, root_probabilities = rp, include_tips = TRUE, include_nodes = FALSE,
                                   Nsimulations = batch, drop_dims = FALSE)$tip_states
    np <- rowSums(s == 2L); ok <- which(abs(np - n_obs) <= tol * n_obs)
    if (length(ok)) { i <- ok[1]; return(list(states = setNames(as.integer(s[i, ]) - 1L, tree$tip.label), draws = used + i)) }
    used <- used + batch
  }
  stop("no simulated vector within the prevalence window after ", used, " draws")
}

# Refit ARD on simulated tip states, reconstruct, and return edge flags. Constant tip states cannot be fitted:
# no event is called and `degenerate` is TRUE. Any other fit failure also counts as degenerate.
sim_edge_flags <- function(tree, states) {
  ne <- nrow(tree$edge)
  none <- list(gain = rep(FALSE, ne), loss = rep(FALSE, ne), degenerate = TRUE)
  if (length(unique(states)) < 2) return(none)
  f <- tryCatch(asr_binary(tree, states, "ARD", Ntrials = 3), error = function(e) NULL)
  if (is.null(f)) return(none)
  ev <- call_events(tree, f$p1)
  list(gain = ev$class == "gain", loss = ev$class == "loss", degenerate = FALSE)
}

# Counts of branches that carry event A in step a and event B in step b.
# flags: named list (by step) of lists with logical `gain` and `loss`. Returns a data.frame, one row per pair x type x scope.
co_count_table <- function(tree, flags) {
  internal <- tree$edge[, 2] > Ntip(tree)
  steps <- names(flags); rows <- list()
  for (i in seq_along(steps)) for (j in seq_along(steps)) if (i < j) for (ty in names(SIM_PAIR_TYPES)) {
    e <- SIM_PAIR_TYPES[[ty]]
    both <- flags[[steps[i]]][[e[1]]] & flags[[steps[j]]][[e[2]]]
    rows[[length(rows) + 1]] <- data.frame(step_a = steps[i], step_b = steps[j], type = ty, n_all = sum(both), n_internal = sum(both & internal),
                                           stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Empirical one-sided p: (1 + number of null values >= observed) / (N + 1). Ties count as exceedances.
empirical_p <- function(observed, null) (1 + sum(null >= observed)) / (length(null) + 1)

# One replicate: simulate every step independently on `tree`, refit, call events, count co-events.
# params: data.frame with columns step, r01, r10 (and root_p1 for root = "asr"). Returns list(counts, steps).
null_replicate <- function(tree, params, root = "stationary", match_tol = 0) {
  flags <- list(); stp <- list()
  for (k in seq_len(nrow(params))) {
    rt <- if (identical(root, "asr")) c(1 - params$root_p1[k], params$root_p1[k]) else if (identical(root, "flat")) c(0.5, 0.5) else "stationary"
    draws <- 1
    if (match_tol > 0) { m <- simulate_binary_matched(tree, params$r01[k], params$r10[k], params$n_present[k], match_tol, rt); st <- m$states; draws <- m$draws
    } else st <- simulate_binary(tree, params$r01[k], params$r10[k], 1, rt)[[1]]
    fl <- sim_edge_flags(tree, st); flags[[params$step[k]]] <- fl
    stp[[k]] <- data.frame(step = params$step[k], n_present = sum(st), draws = draws, gains = sum(fl$gain), losses = sum(fl$loss), degenerate = fl$degenerate,
                           internal_gains = sum(fl$gain & tree$edge[, 2] > Ntip(tree)), internal_losses = sum(fl$loss & tree$edge[, 2] > Ntip(tree)))
  }
  list(counts = co_count_table(tree, flags), steps = do.call(rbind, stp))
}
