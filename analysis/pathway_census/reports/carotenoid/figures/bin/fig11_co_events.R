#!/usr/bin/env Rscript
# Figure 11: co-gain and co-loss of step pairs on the same branches (primary root): observed count and BH q value.
# Usage: Rscript fig11_co_events.R <results dir> <out dir> [small_tree steps] [file stem]
suppressPackageStartupMessages({library(ggplot2)})
a <- commandArgs(TRUE); rd <- a[1]; od <- a[2]
small_tree <- if (length(a) >= 3) strsplit(a[3], ",")[[1]] else c("phytoene_desaturase", "cart")
stem <- if (length(a) >= 4) a[4] else "fig11_co_events"
BLUE <- "#2a78d6"; SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"
ce <- read.delim(file.path(rd, "co_events.tsv"), stringsAsFactors = FALSE)
ce <- ce[ce$root == ce$root[1], ]
cnt <- read.delim(file.path(rd, "gainloss_counts.tsv"), stringsAsFactors = FALSE); cnt <- cnt[cnt$root == ce$root[1], ]
steps <- unique(cnt$step); srcs <- setNames(cnt$source[match(steps, cnt$step)], steps)
nm <- function(s) paste0(gsub("_", " ", s), "\n(", ifelse(srcs[s] == "hits", "hits", ifelse(s %in% small_tree, "placed, small tree", "placed, named members")), ")")
d <- ce[(ce$event_a == "gain" & ce$event_b == "gain") | (ce$event_a == "loss" & ce$event_b == "loss"), ]
d$panel <- ifelse(d$event_a == "gain", "co-gain (branches with a gain in both steps)", "co-loss (branches with a loss in both steps)")
d$x <- factor(nm(d$step_b), levels = nm(steps)); d$y <- factor(nm(d$step_a), levels = rev(nm(steps)))
d$lab <- sprintf("%d\nexpected %.2f\nq %s", d$observed, d$expected, ifelse(d$q < 0.001, "< 0.001", sprintf("= %.3f", d$q)))
p <- ggplot(d, aes(x, y, fill = observed)) + geom_tile(colour = SURFACE, linewidth = 1) +
  geom_text(aes(label = lab, colour = observed > max(observed) * 0.5), size = 2.7, lineheight = 0.9) +
  scale_colour_manual(values = c("TRUE" = "#ffffff", "FALSE" = INK), guide = "none") +
  scale_fill_gradient(low = "#e8f0fb", high = BLUE, name = "observed\nbranches", limits = c(0, NA)) +
  facet_wrap(~panel) + labs(x = NULL, y = NULL, title = "Events on the same branch for each pair of steps",
    subtitle = "Observed count, count expected under a random draw of the second step's event branches, and BH q over all 24 pair-by-event-type tests.\nPhytoene synthase and lycopene cyclase are domains of one fused protein in most fungi: their overlap is not independent evidence.") +
  theme_minimal(base_size = 9) + theme(plot.background = element_rect(fill = SURFACE, colour = NA), panel.background = element_rect(fill = SURFACE, colour = NA),
    text = element_text(colour = INK), panel.grid = element_blank(), axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(od, paste0(stem, ".png")), p, width = 10, height = 4.8, dpi = 150, bg = SURFACE)
ggsave(file.path(od, paste0(stem, ".pdf")), p, width = 10, height = 4.8, bg = SURFACE)
