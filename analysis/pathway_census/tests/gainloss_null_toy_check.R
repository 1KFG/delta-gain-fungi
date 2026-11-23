# Toy checks for bin/gainloss_null_lib.R. Prints key=value lines; tests/test_gainloss_null.py asserts them.
args <- commandArgs(FALSE); f <- sub("--file=", "", args[grep("--file=", args)]); bin <- file.path(dirname(normalizePath(f)), "..", "bin")
source(file.path(bin, "gainloss_lib.R")); source(file.path(bin, "gainloss_null_lib.R"))
bal <- function(tips, stem) {
  one <- function(t) if (length(t) == 1) paste0(t, ":0.05") else { h <- length(t) / 2; paste0("(", one(t[1:h]), ",", one(t[(h + 1):length(t)]), "):0.05") }
  sub(":0.05$", paste0(":", stem), one(tips))
}
tr <- read.tree(text = paste0("((", bal(paste0("O", 1:16), 0.05), ",(", bal(paste0("P1_", 1:8), 0.05), ",", bal(paste0("A", 1:8), 0.5), "):0.05):0.05,(",
                              bal(paste0("P2_", 1:8), 0.05), ",", bal(paste0("B", 1:8), 0.5), "):0.05);"))
cat("ntip=", Ntip(tr), "\n", sep = "")
# 1. number of tips and names
s <- simulate_binary(tr, 0.5, 1.0, n = 5, root = "stationary")
cat("sim_n=", length(s), "\nsim_len=", length(s[[1]]), "\nsim_names_ok=", as.integer(identical(names(s[[1]]), tr$tip.label)), "\n", sep = "")
# 2. same seed gives the same draw; zero rates keep the root state
set.seed(11); a1 <- simulate_binary(tr, 0.5, 1.0, 1)[[1]]; set.seed(11); a2 <- simulate_binary(tr, 0.5, 1.0, 1)[[1]]
cat("same_seed_identical=", as.integer(identical(a1, a2)), "\n", sep = "")
z <- simulate_binary(tr, 0, 0, 3, root = c(0, 1))
cat("zero_rate_all_present=", as.integer(all(sapply(z, function(v) all(v == 1)))), "\n", sep = "")
z0 <- simulate_binary(tr, 0, 0, 3, root = c(1, 0)); cat("zero_rate_all_absent=", as.integer(all(sapply(z0, function(v) all(v == 0)))), "\n", sep = "")
# 3. fixed-seed regression: the sum of the simulated states under seed 1 is a fixed number
set.seed(1); r <- simulate_binary(tr, 0.5, 1.0, 1)[[1]]; cat("seed1_sum=", sum(r), "\n", sep = "")
# 4. known event counts: fixed tip states (as in gainloss_toy_check.R: A and B present, noise tip O1) give 3 gain edges; counting function
S <- function(O, P1, A, P2, B) { v <- c(rep(O, 16), rep(P1, 8), rep(A, 8), rep(P2, 8), rep(B, 8)); names(v) <- c(paste0("O", 1:16), paste0("P1_", 1:8), paste0("A", 1:8), paste0("P2_", 1:8), paste0("B", 1:8)); v[tr$tip.label] }
u <- S(0, 0, 1, 0, 1); u["O1"] <- 1; w <- S(0, 0, 1, 0, 0); w["O1"] <- 1
fu <- sim_edge_flags(tr, u); fw <- sim_edge_flags(tr, w)
cat("U_gains=", sum(fu$gain), "\nW_gains=", sum(fw$gain), "\nU_degenerate=", as.integer(fu$degenerate), "\n", sep = "")
ct <- co_count_table(tr, list(U = fu, W = fw)); g <- ct[ct$type == "gain_gain", ]
cat("UW_cogain_all=", g$n_all, "\nUW_cogain_internal=", g$n_internal, "\n", sep = "")   # shared: stem of A (internal) + tip O1 -> all 2, internal 1
# constant states are degenerate with no events
fc <- sim_edge_flags(tr, setNames(rep(0L, Ntip(tr)), tr$tip.label)); cat("const_degenerate=", as.integer(fc$degenerate), "\nconst_events=", sum(fc$gain) + sum(fc$loss), "\n", sep = "")
# 5. empirical p
cat("p_zero_exceed=", empirical_p(10, c(0, 1, 2, 3)), "\n", sep = "")      # 1/5
cat("p_ties=", empirical_p(2, c(0, 2, 2, 3)), "\n", sep = "")              # (1+3)/5
cat("p_all_exceed=", empirical_p(0, c(0, 0, 0, 0)), "\n", sep = "")        # 1
# 6. replicate pipeline runs and returns counts for 1 pair x 4 types
pr <- data.frame(step = c("X", "Y"), r01 = c(0.3, 0.3), r10 = c(1, 1), root_p1 = c(0.2, 0.2))
set.seed(2); rep1 <- null_replicate(tr, pr); cat("rep_rows=", nrow(rep1$counts), "\nrep_steps=", nrow(rep1$steps), "\n", sep = "")
# 7. prevalence-matched simulation stays inside the window
set.seed(3); m <- simulate_binary_matched(tr, 0.3, 1, n_obs = 10, tol = 0.2); cat("matched_in_window=", as.integer(abs(sum(m$states) - 10) <= 2), "\nmatched_len=", length(m$states), "\n", sep = "")
