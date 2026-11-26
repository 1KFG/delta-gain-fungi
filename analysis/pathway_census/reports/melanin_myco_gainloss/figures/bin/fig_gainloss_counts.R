# NOTE: generic copy of fig10_gain_loss_counts_5steps.R for melanin and mycosporine steps.
#!/usr/bin/env Rscript
# Figure 10: gains, losses and uncertain branches per step (primary root), both Mk models, with parsimony gains as markers.
# Usage: Rscript fig_gainloss_counts.R <results dir> <out dir> <file stem> [ring note TSV: step, note]
suppressPackageStartupMessages({library(ggplot2)})
a <- commandArgs(TRUE); rd <- a[1]; od <- a[2]
stem <- a[3]
notes <- if (length(a) >= 4) setNames(read.delim(a[4], stringsAsFactors = FALSE)$note, read.delim(a[4], stringsAsFactors = FALSE)$step) else character(0)
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"; SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"
cnt <- read.delim(file.path(rd, "gainloss_counts.tsv"), stringsAsFactors = FALSE)
root <- cnt$root[1]; cnt <- cnt[cnt$root == root, ]
srcnote <- function(step, source) ifelse(step %in% names(notes), notes[step], "placed")
cnt$step_lab <- paste0(gsub("_", " ", cnt$step), "\n(", srcnote(cnt$step, cnt$source), ")")
cnt$step_lab <- factor(cnt$step_lab, levels = rev(unique(cnt$step_lab)))
cnt$panel <- paste0(cnt$model, " model"); cnt$best <- cnt$best_by_AIC %in% c(1, TRUE, "TRUE")
long <- rbind(data.frame(cnt[c("step_lab", "panel")], kind = "gain (0 to 1)", n = cnt$gains),
              data.frame(cnt[c("step_lab", "panel")], kind = "loss (1 to 0)", n = cnt$losses),
              data.frame(cnt[c("step_lab", "panel")], kind = "uncertain branch", n = cnt$uncertain))
long$kind <- factor(long$kind, levels = c("gain (0 to 1)", "loss (1 to 0)", "uncertain branch"))
nudge <- -0.8 / 3
pars <- rbind(data.frame(cnt[c("step_lab", "panel")], kind = "parsimony gains, late changes", n = cnt$pars_gains_late),
              data.frame(cnt[c("step_lab", "panel")], kind = "parsimony gains, early changes", n = cnt$pars_gains_early))
gl <- cnt[c("step_lab", "panel")]; gl$xmax <- pmax(cnt$gains, cnt$pars_gains_late, cnt$pars_gains_early)
gl$lab <- sprintf("%d (parsimony %d to %d)%s", cnt$gains, pmin(cnt$pars_gains_late, cnt$pars_gains_early), pmax(cnt$pars_gains_late, cnt$pars_gains_early), ifelse(cnt$best, "   [lower AIC]", ""))
p <- ggplot(long, aes(y = step_lab, x = n, fill = kind)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(data = long[long$kind != "gain (0 to 1)", ], aes(y = as.numeric(step_lab) + ifelse(kind == "loss (1 to 0)", 0, -nudge), x = n, label = n), inherit.aes = FALSE, hjust = -0.2, size = 2.6, colour = INK) +
  geom_text(data = gl, aes(y = as.numeric(step_lab) + nudge, x = xmax, label = lab), inherit.aes = FALSE, hjust = -0.25, size = 2.6, colour = INK) +
  geom_point(data = pars, aes(y = as.numeric(step_lab) + nudge, x = n, shape = kind), inherit.aes = FALSE, size = 2.6, colour = INK, fill = SURFACE) +
  scale_fill_manual(values = c(BLUE, ORANGE, AQUA), name = NULL) +
  scale_shape_manual(values = c(18, 5), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  facet_wrap(~panel, nrow = 1) + labs(x = "number of branches (or origins)", y = NULL,
    title = "Gains, losses and uncertain branches per step",
    subtitle = paste0("Mk model on the genus tree (rooted on ", root, "), ", prettyNum(cnt$n_tips[1], big.mark = ","), " tips with data. A branch is an event only when both ends have probability >= 0.8.\nDiamonds: gains in one parsimony reconstruction (tie rule late or early), not a model estimate."))
p <- p + theme_minimal(base_size = 9) + theme(plot.background = element_rect(fill = SURFACE, colour = NA), panel.background = element_rect(fill = SURFACE, colour = NA),
  text = element_text(colour = INK), legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.subtitle = element_text(size = 7.5), panel.spacing.x = unit(1.6, "lines"))
p <- p + labs(subtitle = paste0("Microsporidia excluded (30 tips); tree rooted on Blastocladiomycota.\n", p$labels$subtitle))
ggsave(file.path(od, paste0(stem, ".png")), p, width = 9.5, height = max(5, 1.0 + 0.62 * length(unique(cnt$step))), dpi = 150, bg = SURFACE)
ggsave(file.path(od, paste0(stem, ".pdf")), p, width = 9.5, height = max(5, 1.0 + 0.62 * length(unique(cnt$step))), bg = SURFACE)
