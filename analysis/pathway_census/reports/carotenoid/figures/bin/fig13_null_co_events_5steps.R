# NOTE: five-step copy of fig13_null_co_events.R for the Microsporidia-excluded, Blastocladiomycota-rooted primary run (adds the exclusion and root to the subtitle).
#!/usr/bin/env Rscript
# Figure 13: tree-aware simulation null for co-events, with the observed count marked.
# Usage: Rscript fig13_null_co_events.R <null dir> <variant tag> <out dir> <file stem> "<dataset label>"
# Rows: step pairs. Columns: co-gain and co-loss, each for all edges and for internal edges only.
suppressPackageStartupMessages({library(ggplot2)})
a <- commandArgs(TRUE); nd <- a[1]; vr <- a[2]; od <- a[3]; stem <- a[4]; ds <- a[5]
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"
sm <- read.delim(file.path(nd, "null_summary.tsv"), stringsAsFactors = FALSE); sm <- sm[sm$variant == vr & sm$type %in% c("gain_gain", "loss_loss"), ]
fs <- list.files(nd, pattern = paste0("^sim_counts\\.", vr, "\\.task[0-9]+\\.tsv$"), full.names = TRUE)
cs <- do.call(rbind, lapply(fs, read.delim, stringsAsFactors = FALSE)); cs <- cs[cs$type %in% c("gain_gain", "loss_loss"), ]
N <- length(unique(cs$rep))
dg <- read.delim(file.path(nd, "null_step_diagnostics.tsv"), stringsAsFactors = FALSE); dg <- dg[dg$variant == vr, ]
prev <- sprintf("Present tips per replicate (mean over replicates, 5 steps): %d to %d. Observed: %d to %d.", round(min(dg$null_present_mean)), round(max(dg$null_present_mean)), min(dg$obs_present), max(dg$obs_present))
pl <- function(x) gsub("_", " ", x)
long <- rbind(data.frame(cs[c("step_a", "step_b", "type")], scope = "all_edges", n = cs$n_all), data.frame(cs[c("step_a", "step_b", "type")], scope = "internal_edges", n = cs$n_internal))
key <- function(d) paste(d$step_a, d$step_b, d$type, d$scope)
tab <- aggregate(rep ~ step_a + step_b + type + scope + n, data = cbind(long, rep = 1), FUN = sum); names(tab)[names(tab) == "rep"] <- "freq"
lab <- function(d) {
  d$pair <- paste0(pl(d$step_a), " x\n", pl(d$step_b)); d$pair[d$step_a == "phytoene_synthase" & d$step_b == "lycopene_cyclase"] <- paste0("phytoene synthase x\nlycopene cyclase\n(one fused gene: not a meaningful test)")
  d$col <- paste0(ifelse(d$type == "gain_gain", "co-gain", "co-loss"), ", ", ifelse(d$scope == "all_edges", "all edges", "internal edges only")); d }
tab <- lab(tab); sm <- lab(sm)
ord <- rev(unique(tab$pair)); tab$pair <- factor(tab$pair, levels = ord); sm$pair <- factor(sm$pair, levels = ord)
cols <- c("co-gain, all edges", "co-gain, internal edges only", "co-loss, all edges", "co-loss, internal edges only")
tab$col <- factor(tab$col, levels = cols); sm$col <- factor(sm$col, levels = cols)
sm$txt <- sprintf("observed %d\nnull %.2f (%d-%d)\np = %.4f", sm$observed, sm$null_mean, round(sm$null_lo), round(sm$null_hi), sm$p)
tab$prop <- tab$freq / N
p <- ggplot() + geom_col(data = tab, aes(n, prop), fill = BLUE, width = 0.8) +
  geom_segment(data = sm, aes(x = observed, xend = observed, y = 0, yend = 1.05), colour = ORANGE, linewidth = 0.9) +
  geom_text(data = sm, aes(x = Inf, y = 1.62, label = txt), hjust = 1.05, vjust = 1, size = 2.2, colour = INK, lineheight = 0.9) +
  facet_grid(pair ~ col, scales = "free_x", switch = "y") + scale_y_continuous(limits = c(0, 1.65), breaks = c(0, 0.5, 1), labels = function(x) sprintf("%.1f", x), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "branches with the event in both steps", y = "share of replicates", title = paste0("Co-events against a tree-aware simulation null (", ds, ")"),
       subtitle = paste0("Blue bars: ", N, " replicates, steps simulated independently on the same tree, each refitted and called with the 0.8 rule. Orange line: observed count.\nThe null does not model missing data, a shared technical cause or any link between steps.\n", prev)) +
  theme_minimal(base_size = 8) + theme(plot.background = element_rect(fill = SURFACE, colour = NA), panel.background = element_rect(fill = SURFACE, colour = NA),
    text = element_text(colour = INK), strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.5), strip.text.x = element_text(size = 7.5), panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), plot.subtitle = element_text(size = 7), panel.spacing = unit(0.5, "lines"))
p <- p + labs(subtitle = paste0("Microsporidia excluded (30 tips); tree rooted on Blastocladiomycota.\n", p$labels$subtitle))
ggsave(file.path(od, paste0(stem, ".png")), p, width = 11, height = 13.5, dpi = 150, bg = SURFACE)
ggsave(file.path(od, paste0(stem, ".pdf")), p, width = 11, height = 13.5, bg = SURFACE)
