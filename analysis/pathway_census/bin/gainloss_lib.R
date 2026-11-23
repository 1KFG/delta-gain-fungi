# Gain and loss helpers for binary traits on a tree. Sourced by run_gainloss.R and by the toy-tree test.
# Dependencies: ape, castor.
# State coding in this file: tip and node probability `p1` = P(state 1 = present). States in castor are 1 (absent) and 2 (present).

suppressPackageStartupMessages({library(ape); library(castor)})

PROB_CUT <- 0.8

# Keep tips with a non-NA state. `states` is a named 0/1/NA vector. Returns list(tree, states).
prune_to_data <- function(tree, states) {
  states <- states[tree$tip.label]
  keep <- names(states)[!is.na(states)]
  if (length(keep) < 3) stop("fewer than 3 tips with data")
  list(tree = if (length(keep) < Ntip(tree)) keep.tip(tree, keep) else tree, states = states[keep])
}

# Root on a set of tips that form one side of a split of the unrooted tree. The root goes in the middle of the edge.
# Stops if the set is not a clade of the unrooted tree.
root_on_split <- function(tree, tips) {
  tips <- intersect(tips, tree$tip.label)
  if (!length(tips)) stop("no outgroup tips in tree")
  r <- ape::root(tree, outgroup = tips, resolve.root = TRUE)
  if (!is.rooted(r)) stop("rooting failed")
  r
}

# Largest split of the unrooted tree whose side holds only tips from `allowed`. Returns the tip vector (or NULL).
largest_pure_split <- function(tree, allowed) {
  u <- unroot(tree); pp <- prop.part(u, check.labels = TRUE); labs <- attr(pp, "labels"); best <- NULL
  for (i in seq_along(pp)) {
    for (s in list(labs[pp[[i]]], setdiff(labs, labs[pp[[i]]]))) {
      if (length(s) >= 2 && all(s %in% allowed) && (is.null(best) || length(s) > length(best))) best <- s
    }
  }
  best
}

# Internal node ids ordered so that every node comes after all its descendants (postorder of the edge table).
children_first <- function(tree) unique(reorder(tree, "postorder")$edge[, 1])

# Two-state transition probabilities for rates r01 (0->1) and r10 (1->0) over length t (vectorised over t).
# Returns a list of the four vectors p00, p01, p10, p11.
tp_2state <- function(r01, r10, t) {
  s <- r01 + r10
  if (s == 0) return(list(p00 = rep(1, length(t)), p01 = rep(0, length(t)), p10 = rep(0, length(t)), p11 = rep(1, length(t))))
  e <- exp(-s * t)
  list(p00 = (r10 + r01 * e) / s, p01 = r01 * (1 - e) / s, p10 = r10 * (1 - e) / s, p11 = (r01 + r10 * e) / s)
}

# Full marginal ancestral probabilities for a 2-state continuous-time Markov model (up and down pass).
# This is used because castor's re-rooting option gives wrong marginals for a model that is not time-reversible
# (checked on a toy tree: all node probabilities were 0 for ARD), and the non-rerooted option returns subtree-conditional
# likelihoods only. Root prior is flat (0.5, 0.5). Returns P(state 1) for tips then nodes (ape node ids).
marginal_2state <- function(tree, tip_p1, r01, r10) {
  n <- Ntip(tree); N <- n + tree$Nnode; ed <- tree$edge; len <- tree$edge.length
  P <- tp_2state(r01, r10, len)
  Lin <- matrix(1, N, 2)                   # inside likelihood, normalised per node
  Lin[1:n, 1] <- 1 - tip_p1; Lin[1:n, 2] <- tip_p1
  kids <- split(seq_len(nrow(ed)), ed[, 1])  # edge indices by parent
  # message from child edge e up to its parent: m_i = sum_j P_ij L_c(j)
  msg <- matrix(1, nrow(ed), 2)
  nodes <- children_first(tree)
  for (u in nodes) {
    es <- kids[[as.character(u)]]; acc <- c(1, 1)
    for (e in es) {
      c <- ed[e, 2]
      m <- c(P$p00[e] * Lin[c, 1] + P$p01[e] * Lin[c, 2], P$p10[e] * Lin[c, 1] + P$p11[e] * Lin[c, 2])
      msg[e, ] <- m; acc <- acc * m
    }
    Lin[u, ] <- acc / sum(acc)
  }
  Out <- matrix(1, N, 2); root <- n + 1L; Out[root, ] <- c(0.5, 0.5)
  for (u in rev(nodes)) {                  # parents before children
    es <- kids[[as.character(u)]]
    for (e in es) {
      oth <- Out[u, ]
      for (e2 in setdiff(es, e)) oth <- oth * msg[e2, ]
      c <- ed[e, 2]
      o <- c(oth[1] * P$p00[e] + oth[2] * P$p10[e], oth[1] * P$p01[e] + oth[2] * P$p11[e])
      Out[c, ] <- o / sum(o)
    }
  }
  post <- Lin * Out; post <- post / rowSums(post)
  post[, 2]
}

# Fit a binary Mk model with castor (rates and log-likelihood), then compute full marginals with marginal_2state.
# tip_states: named 0/1 vector without NA. Root prior flat in both the fit and the marginals.
asr_binary <- function(tree, states, rate_model = c("ER", "ARD"), Ntrials = 5) {
  rate_model <- match.arg(rate_model)
  ts <- as.integer(states[tree$tip.label]) + 1L
  fit <- castor::asr_mk_model(tree, tip_states = ts, Nstates = 2, rate_model = rate_model, root_prior = "flat",
                              Ntrials = Ntrials, include_ancestral_likelihoods = FALSE)
  if (!isTRUE(fit$success)) stop("asr_mk_model failed: ", fit$error)
  npar <- if (rate_model == "ER") 1 else 2
  r01 <- fit$transition_matrix[1, 2]; r10 <- fit$transition_matrix[2, 1]
  tp <- as.numeric(states[tree$tip.label])
  p1 <- c(tp, marginal_2state(tree, tp, r01, r10)[-seq_along(tp)])
  list(p1 = p1, loglik = fit$loglikelihood, AIC = 2 * npar - 2 * fit$loglikelihood, npar = npar, rate_01 = r01, rate_10 = r10)
}

# Call events per edge from node probabilities p1 (tips then nodes).
# state called 1 if p1 >= cut, 0 if p1 <= 1 - cut, otherwise NA.
# Edge classes: gain (0->1), loss (1->0), same (called, no change), uncertain (parent or child not called).
call_events <- function(tree, p1, cut = PROB_CUT) {
  called <- ifelse(p1 >= cut, 1L, ifelse(p1 <= 1 - cut, 0L, NA_integer_))
  par <- called[tree$edge[, 1]]; chi <- called[tree$edge[, 2]]
  cls <- ifelse(is.na(par) | is.na(chi), "uncertain",
         ifelse(par == 0L & chi == 1L, "gain", ifelse(par == 1L & chi == 0L, "loss", "same")))
  list(class = cls, parent_called = par, child_called = chi)
}

# Sankoff minimum cost for a binary trait. states: named 0/1 vector. cost_gain = cost of 0->1, cost_loss = cost of 1->0.
# root_state: NA (free), 0 or 1 (forced). Returns the minimum total cost.
sankoff_cost <- function(tree, states, cost_gain = 1, cost_loss = 1, root_state = NA) {
  n <- Ntip(tree); N <- n + tree$Nnode; ed <- tree$edge
  S <- matrix(0, N, 2); st <- states[tree$tip.label]
  S[1:n, 1] <- ifelse(st == 0, 0, Inf); S[1:n, 2] <- ifelse(st == 1, 0, Inf)
  kids <- split(ed[, 2], ed[, 1]); nodes <- children_first(tree)
  for (u in nodes) {
    cs <- kids[[as.character(u)]]
    S[u, 1] <- sum(vapply(cs, function(c) min(S[c, 1], S[c, 2] + cost_gain), 0))
    S[u, 2] <- sum(vapply(cs, function(c) min(S[c, 2], S[c, 1] + cost_loss), 0))
  }
  r <- n + 1L
  if (is.na(root_state)) min(S[r, ]) else S[r, root_state + 1]
}

# One optimal Sankoff reconstruction (equal costs, root free) with a tie rule; returns the number of gains and losses.
#   tie = "late": when both states are optimal for a child, keep the parent's state (changes move toward the tips).
#   tie = "early": when both are optimal, take the state that differs from the parent (changes move toward the root).
# Root ties go to state 0. The sum gains + losses is the minimum number of changes.
sankoff_gains_losses <- function(tree, states, tie = c("late", "early")) {
  tie <- match.arg(tie)
  n <- Ntip(tree); N <- n + tree$Nnode; ed <- tree$edge
  S <- matrix(0, N, 2); st <- states[tree$tip.label]
  S[1:n, 1] <- ifelse(st == 0, 0, Inf); S[1:n, 2] <- ifelse(st == 1, 0, Inf)
  kids <- split(ed[, 2], ed[, 1]); nodes <- children_first(tree)
  for (u in nodes) {
    cs <- kids[[as.character(u)]]
    S[u, 1] <- sum(vapply(cs, function(c) min(S[c, 1], S[c, 2] + 1), 0))
    S[u, 2] <- sum(vapply(cs, function(c) min(S[c, 2], S[c, 1] + 1), 0))
  }
  state <- integer(N); r <- n + 1L; state[r] <- if (S[r, 2] < S[r, 1]) 1L else 0L
  gains <- 0L; losses <- 0L
  for (u in rev(nodes)) for (c in kids[[as.character(u)]]) {
    p <- state[u]
    cost <- c(S[c, 1] + (p == 1L), S[c, 2] + (p == 0L))   # total cost if child is 0 or 1
    if (cost[1] == cost[2]) { cs <- if (tie == "late") p else 1L - p } else cs <- as.integer(which.min(cost) - 1L)
    state[c] <- cs
    if (p == 0L && cs == 1L) gains <- gains + 1L
    if (p == 1L && cs == 0L) losses <- losses + 1L
  }
  list(gains = gains, losses = losses)
}

# Parsimony counts. Sankoff is implemented here; the minimum number of changes is cross-checked against castor::asr_max_parsimony.
#   changes, changes_castor: minimum number of changes (equal costs, root free).
#   gains_late, losses_late, gains_early, losses_early: gains and losses in one optimal reconstruction under two tie rules.
# The origins of a trait that can be lost are not defined by one number: the two tie rules bracket them (not all optimal
# reconstructions are covered).
parsimony_counts <- function(tree, states) {
  ts <- as.integer(states[tree$tip.label]) + 1L
  cc <- castor::asr_max_parsimony(tree, ts, Nstates = 2, transition_costs = "all_equal")$total_cost
  l <- sankoff_gains_losses(tree, states, "late"); e <- sankoff_gains_losses(tree, states, "early")
  list(changes = sankoff_cost(tree, states, 1, 1), changes_castor = cc,
       gains_late = l$gains, losses_late = l$losses, gains_early = e$gains, losses_early = e$losses)
}

# Per-step summary row from an ASR fit.
summarise_events <- function(ev) {
  tb <- table(factor(ev$class, levels = c("gain", "loss", "same", "uncertain")))
  list(gains = as.integer(tb["gain"]), losses = as.integer(tb["loss"]), same = as.integer(tb["same"]),
       uncertain = as.integer(tb["uncertain"]), n_edges = length(ev$class))
}

# Permutation test for the overlap of two event sets among edges.
#   a_ev, b_ev: logical vectors over edges (the observed event edges of each step).
#   b_elig: logical vector, edges where step B could have that event (parent called in the state that allows it).
# Null: A is fixed. B's events are redrawn as a uniform sample, without replacement, of size sum(b_ev) from the eligible edges of B.
# Returns observed overlap, expected overlap (hypergeometric mean), mean of permutations, one-sided p (overlap >= observed).
perm_overlap <- function(a_ev, b_ev, b_elig, nperm = 10000, seed = 1) {
  stopifnot(all(b_ev <= b_elig))
  nA <- sum(a_ev); nB <- sum(b_ev); nE <- sum(b_elig)
  obs <- sum(a_ev & b_ev)
  out <- list(n_a = nA, n_b = nB, n_elig_b = nE, observed = obs, expected = NA_real_, perm_mean = NA_real_, p = NA_real_)
  if (nA == 0 || nB == 0) return(out)
  set.seed(seed)
  idx <- which(b_elig)
  a_in <- a_ev[idx]  # A's events among B-eligible edges
  out$expected <- sum(a_in) * nB / nE
  sims <- vapply(seq_len(nperm), function(i) sum(a_in[sample.int(nE, nB)]), numeric(1))
  out$perm_mean <- mean(sims)
  out$p <- (1 + sum(sims >= obs)) / (1 + nperm)
  out
}

# Eligible edges for an event type: gain needs parent called 0; loss needs parent called 1.
eligible_edges <- function(ev, type = c("gain", "loss")) {
  type <- match.arg(type)
  if (type == "gain") ev$parent_called %in% 0L else ev$parent_called %in% 1L
}

# Number of tips below each edge's child node, and the majority label of those tips.
edge_clade_summary <- function(tree, labels) {
  n <- Ntip(tree); ch <- tree$edge[, 2]
  pp <- prop.part(tree, check.labels = TRUE)
  tl <- attr(pp, "labels")
  desc <- function(node) if (node <= n) tree$tip.label[node] else tl[pp[[node - n]]]
  res <- lapply(ch, function(node) {
    tips <- desc(node); lab <- labels[tips]; lab[is.na(lab) | lab == ""] <- "unlabelled"
    tb <- sort(table(lab), decreasing = TRUE)
    list(n_tips = length(tips), top_label = names(tb)[1], top_frac = as.numeric(tb[1]) / length(tips))
  })
  data.frame(n_tips = sapply(res, `[[`, "n_tips"), top_label = sapply(res, `[[`, "top_label"),
             top_frac = sapply(res, `[[`, "top_frac"), stringsAsFactors = FALSE)
}
