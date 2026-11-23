# Toy-tree check of bin/gainloss_lib.R. Prints "key=value" lines; tests/test_gainloss_r.py asserts them.
args <- commandArgs(FALSE); f <- sub("--file=", "", args[grep("--file=", args)])
source(file.path(dirname(normalizePath(f)), "..", "bin", "gainloss_lib.R"))
# Toy tree: root = ((O,G),H). O has 16 tips (stem 0.05). G = (P1 stem 0.05, A stem 0.5). H = (P2 stem 0.05, B stem 0.5).
# P1, P2, A and B have 8 tips each. Inside a clade, branch length is 0.05. The long stems of A and B and the short stems of
# P1 and P2 make the node G and H follow P1 and P2, so a trait that differs between P and A (or B) has a certain gain or loss on the stem.
bal <- function(tips, stem) {
  one <- function(t) if (length(t) == 1) paste0(t, ":0.05") else { h <- length(t) / 2; paste0("(", one(t[1:h]), ",", one(t[(h + 1):length(t)]), "):0.05") }
  sub(":0.05$", paste0(":", stem), one(tips))
}
nw <- paste0("((", bal(paste0("O", 1:16), 0.05), ",(", bal(paste0("P1_", 1:8), 0.05), ",", bal(paste0("A", 1:8), 0.5), "):0.05):0.05,(",
             bal(paste0("P2_", 1:8), 0.05), ",", bal(paste0("B", 1:8), 0.5), "):0.05);")
tr <- read.tree(text = nw)
tips <- tr$tip.label
S <- function(O, P1, A, P2, B) setNames(c(rep(O, 16), rep(P1, 8), rep(A, 8), rep(P2, 8), rep(B, 8))[match(tips, c(paste0("O", 1:16), paste0("P1_", 1:8), paste0("A", 1:8), paste0("P2_", 1:8), paste0("B", 1:8)))], tips)
cat("ntip=", Ntip(tr), "\n", sep = "")
# Check 1: marginal_2state equals castor's re-rooted marginals for a reversible (ER) model.
s <- S(0, 0, 1, 0, 0); s["O1"] <- 1; s["P2_1"] <- 1   # a little noise so that the rate is identifiable
fit <- castor::asr_mk_model(tr, as.integer(s[tr$tip.label]) + 1L, Nstates = 2, rate_model = "ER", root_prior = "flat", reroot = TRUE)
mine <- marginal_2state(tr, as.numeric(s[tr$tip.label]), fit$transition_matrix[1, 2], fit$transition_matrix[2, 1])
cat("marginal_max_abs_diff_vs_castor=", signif(max(abs(mine[-seq_len(Ntip(tr))] - fit$ancestral_likelihoods[, 2])), 3), "\n", sep = "")
ev_of <- function(s, model = "ER") call_events(tr, asr_binary(tr, s, model)$p1)
# add_noise flips tip O1. This makes the rate identifiable and adds exactly one certain event: a tip edge below the absent (or present) node O.
add_noise <- function(s) { s["O1"] <- 1 - s["O1"]; s }
# Trait G1: present in A only. Expected: 2 gains (stem of A, tip O1), 0 losses.
e <- summarise_events(ev_of(add_noise(S(0, 0, 1, 0, 0)))); cat("A_gains=", e$gains, "\nA_losses=", e$losses, "\n", sep = "")
# Trait AB: present in A and B. Expected: 3 gains (stems of A and B, tip O1), 0 losses.
e <- summarise_events(ev_of(add_noise(S(0, 0, 1, 0, 1)))); cat("AB_gains=", e$gains, "\nAB_losses=", e$losses, "\n", sep = "")
# Trait LB: present everywhere except B. Expected: 2 losses (stem of B, tip O1), 0 gains.
e <- summarise_events(ev_of(add_noise(S(1, 1, 1, 1, 0)))); cat("LB_gains=", e$gains, "\nLB_losses=", e$losses, "\n", sep = "")
# Parsimony. AB: 2 changes (2 gains). LB: 1 change (1 loss).
Y <- S(0, 0, 1, 0, 1); L <- S(1, 1, 1, 1, 0)
pY <- parsimony_counts(tr, Y); pL <- parsimony_counts(tr, L)
cat("AB_pars_changes=", pY$changes, "\nAB_pars_changes_castor=", pY$changes_castor, "\nAB_pars_gains_late=", pY$gains_late, "\nAB_pars_gains_early=", pY$gains_early, "\nAB_pars_losses=", pY$losses_late + pY$losses_early, "\n", sep = "")
cat("LB_pars_changes=", pL$changes, "\nLB_pars_losses_late=", pL$losses_late, "\nLB_pars_gains=", pL$gains_late + pL$gains_early, "\n", sep = "")
# A trait with an ambiguous reconstruction: present in P1 and B only (P1 sister to A, B sister to P2). Equal-cost minimum is 2 changes.
T2 <- S(0, 1, 0, 0, 1); p2 <- parsimony_counts(tr, T2)
cat("T2_pars_changes=", p2$changes, "\nT2_sum_late=", p2$gains_late + p2$losses_late, "\nT2_sum_early=", p2$gains_early + p2$losses_early, "\n", sep = "")
# Co-gain. U and V are both present in A and B (3 gain edges each: stem A, stem B, tip O1). W is present in A only (stem A, tip O1): U and W share 2.
gU <- ev_of(add_noise(S(0, 0, 1, 0, 1)))$class == "gain"; gV <- gU; gW <- ev_of(add_noise(S(0, 0, 1, 0, 0)))$class == "gain"
cat("U_V_cogain=", sum(gU & gV), "\nU_W_cogain=", sum(gU & gW), "\n", sep = "")
# Permutation overlap with every gain edge shared: eligible edges are the edges whose parent is called absent.
elig <- eligible_edges(ev_of(add_noise(S(0, 0, 1, 0, 1))), "gain")
po <- perm_overlap(gU, gV, elig, nperm = 2000)
cat("perm_observed=", po$observed, "\nperm_n_b=", po$n_b, "\nperm_expected_below_observed=", as.integer(po$expected < po$observed), "\nperm_p_below_0.05=", as.integer(po$p < 0.05), "\n", sep = "")
pz <- perm_overlap(rep(FALSE, length(gU)), gV, elig, nperm = 10); cat("perm_p_na_when_no_events=", as.integer(is.na(pz$p)), "\n", sep = "")
# NA handling: prune_to_data drops NA tips.
sNA <- Y; sNA[c("A1", "P2_1")] <- NA
pr <- prune_to_data(tr, sNA); cat("pruned_ntip=", Ntip(pr$tree), "\n", sep = "")
