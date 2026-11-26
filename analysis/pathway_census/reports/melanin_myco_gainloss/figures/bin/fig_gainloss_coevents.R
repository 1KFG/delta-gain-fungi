#!/usr/bin/env Rscript
# Co-gain and co-loss of step pairs against the tree-aware simulation null (matrix per event type; all edges).
# Usage: Rscript fig_gainloss_coevents.R <null dir> <variant tag> <out dir> <file stem> "<dataset label>" [exclude steps, comma separated]
# Each cell: observed branches with the event in both steps; null mean; empirical p; BH q (over all pair-by-type tests of the scope, 4 event types per pair).
# A grey cell means that the null mean and the observed count are both 0 (no test possible).
suppressPackageStartupMessages({library(ggplot2)})
a <- commandArgs(TRUE); nd <- a[1]; vr <- a[2]; od <- a[3]; stem <- a[4]; ds <- a[5]
excl <- if (length(a) >= 6 && nzchar(a[6])) strsplit(a[6], ",")[[1]] else character(0)
BLUE <- "#2a78d6"; SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"
sm <- read.delim(file.path(nd, "null_summary.tsv"), stringsAsFactors = FALSE)
sm <- sm[sm$variant == vr & sm$scope == "all_edges" & sm$type %in% c("gain_gain", "loss_loss"), ]
sm <- sm[!(sm$step_a %in% excl) & !(sm$step_b %in% excl), ]
steps <- unique(c(sm$step_a, sm$step_b))
ord <- unique(c(sm$step_a, sm$step_b))
pl <- function(x) gsub("_", " ", x)
sm$panel <- ifelse(sm$type == "gain_gain", "co-gain (branches with a gain in both steps)", "co-loss (branches with a loss in both steps)")
sm$x <- factor(pl(sm$step_b), levels = pl(ord)); sm$y <- factor(pl(sm$step_a), levels = rev(pl(ord)))
sm$lab <- sprintf("%d\nnull %.2f\np %s\nq %s", sm$observed, sm$null_mean,
                  ifelse(sm$p < 0.001, "< 0.001", sprintf("%.3f", sm$p)), ifelse(sm$q < 0.001, "< 0.001", sprintf("%.3f", sm$q)))
sm$none <- sm$observed == 0 & sm$null_mean < 0.005
sm$lab[sm$none] <- "0\n(no test)"
sm$sig <- sm$q < 0.05
cap <- max(sm$observed, 1)
p <- ggplot(sm, aes(x, y)) + geom_tile(aes(fill = ifelse(none, NA, observed)), colour = SURFACE, linewidth = 1) +
  geom_tile(data = sm[sm$sig, ], fill = NA, colour = INK, linewidth = 0.7) +
  geom_text(aes(label = lab, colour = !is.na(observed) & observed > cap * 0.5 & !none), size = 2.1, lineheight = 0.9) +
  scale_colour_manual(values = c("TRUE" = "#ffffff", "FALSE" = INK), guide = "none") +
  scale_fill_gradient(low = "#e8f0fb", high = BLUE, name = "observed\nbranches", limits = c(0, NA), na.value = "#e4e3de") +
  facet_wrap(~panel) + labs(x = NULL, y = NULL, title = paste0("Events on the same branch against a tree-aware null (", ds, ")"),
    subtitle = paste0("Null variant: ", vr, ". Cell: observed branches, null mean, empirical p, BH q over all pair-by-event-type tests of the scope (4 event types per pair, as in null_summary.tsv). Black outline: q < 0.05.\nThe null simulates each step alone on the same tree; it does not model missing data, a shared technical cause or a link between steps.\nMicrosporidia excluded (30 tips); tree rooted on Blastocladiomycota.")) +
  theme_minimal(base_size = 9) + theme(plot.background = element_rect(fill = SURFACE, colour = NA), panel.background = element_rect(fill = SURFACE, colour = NA),
    text = element_text(colour = INK), panel.grid = element_blank(), axis.text.x = element_text(size = 6.5, angle = 30, hjust = 1), axis.text.y = element_text(size = 6.5), plot.subtitle = element_text(size = 7))
n <- length(ord)
ggsave(file.path(od, paste0(stem, ".png")), p, width = max(8, 2 * (1 + 0.9 * n)), height = max(4, 1.8 + 0.8 * n), dpi = 150, bg = SURFACE)
ggsave(file.path(od, paste0(stem, ".pdf")), p, width = max(8, 2 * (1 + 0.9 * n)), height = max(4, 1.8 + 0.8 * n), bg = SURFACE)
