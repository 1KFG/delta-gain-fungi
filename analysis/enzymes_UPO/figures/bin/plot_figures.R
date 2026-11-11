#!/usr/bin/env Rscript
# Figures 1, 3, 4 of the UPO clade summary, from the TSVs written by prep_figures.py.
#
#   Rscript plot_figures.R <analysis/enzymes_UPO> <figures/data> <figures/out>
#
# Fig 1  reference tree (187 leaves), rooted on Mucorales; Zamocky Figure 1 lineages as shaded
#        sectors (11 clades folded into 8 colour groups), our Clades I-V as labelled arcs with
#        SH-aLRT/UFBoot, the 19 characterized enzymes named.
# Fig 3  BFD proteins: old nearest-neighbour call -> placement call, flows coloured by the old
#        call's confidence tier.
# Fig 4  matched-identity validation of both callers: precision when committing (exact 95 % CI)
#        and abstention rate, by identity to the nearest surviving reference.
# Colours: the dataviz reference categorical palette (light), slots in fixed order.
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(ggalluvial)
  library(patchwork); library(dplyr); library(readr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
upo <- args[1]; dat <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

SLOTS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8984"; SURFACE <- "#fcfcfb"
theme_fig <- theme_minimal(base_size = 9) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        text = element_text(colour = INK), axis.text = element_text(colour = INK2),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(colour = INK2, size = 8))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, bg = SURFACE)
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = SURFACE)
}

# ---------------- Figure 1 ----------------
tips <- read_tsv(file.path(dat, "fig1_tips.tsv"), show_col_types = FALSE)
sup  <- read_tsv(file.path(dat, "fig1_support.tsv"), show_col_types = FALSE,
                 col_types = cols(.default = "c"))
tr <- read.tree(file.path(upo, "data/refset_zamocky2022/tree/mfp.treefile"))
muc <- tips$label[which(tips$fig1_clade == "Mucoromycetous basal")]  # which(): empty clade fields read as NA
tr <- root(tr, outgroup = muc, resolve.root = TRUE)
tr <- ladderize(tr)

GROUPS <- c("Mucorales (basal)", "Glomeromycete", "Basidiomycetous 1st", "Ascomycetous 1st",
            "Basidiomycetous 2nd", "Oomycete, chytrid, amoebozoan", "Family II (long)",
            "Family I (short)")
group_col <- setNames(SLOTS, GROUPS)

hl <- tips %>% filter(fig1_clade != "" & !is.na(fig1_clade)) %>%
  group_by(fig1_clade, fig1_group) %>%
  summarise(node = getMRCA(tr, label), .groups = "drop") %>%
  mutate(fig1_group = factor(fig1_group, levels = GROUPS))

ours <- tips %>% filter(!is.na(our_clade) & our_clade != "") %>% group_by(our_clade) %>%
  summarise(node = getMRCA(tr, label), .groups = "drop") %>%
  left_join(sup %>% filter(system == "our") %>% select(clade, sh_alrt, ufboot),
            by = c("our_clade" = "clade")) %>%
  mutate(lab = sprintf("%s  %s/%s", sub("Clade ", "", our_clade), sh_alrt, ufboot))

tips$char <- ifelse(tips$set == "backbone", sub("^B__([^_]+(_[^_]+)*)__clade.*$", "\\1", tips$label), NA)
tips$char <- sub("PaDa-I_parent_AaeUPO", "AaeUPO", tips$char)

p1 <- ggtree(tr, layout = "circular", linewidth = 0.25, colour = INK2) %<+% tips +
  geom_hilight(data = hl, mapping = aes(node = node, fill = fig1_group),
               alpha = 0.28, extend = 0.35, to.bottom = TRUE) +
  scale_fill_manual(values = group_col, name = "Zamocky 2022 lineage", drop = FALSE) +
  geom_tippoint(aes(subset = set == "backbone"), size = 1.6, colour = INK) +
  geom_tiplab(aes(label = char, subset = set == "backbone"), size = 1.8, offset = 0.04,
              colour = INK) +
  geom_cladelab(data = ours, mapping = aes(node = node, label = lab), offset = 0.75, align = TRUE,
                barsize = 0.8, barcolour = INK, textcolour = INK, fontsize = 2.6,
                angle = "auto", hjust = 0.5, offset.text = 0.12) +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.title = element_text(size = 8, face = "bold"),
        plot.background = element_rect(fill = SURFACE, colour = NA)) +
  labs(title = "UPO reference tree: 19 characterized enzymes + 168 curated HTPs (Zamocky 2022)",
       subtitle = paste0("ML tree (IQ-TREE, Q.PFAM+R7), rooted on Mucorales. Sectors: Zamocky lineages.\n",
                         "Arcs: our Clades I-V with SH-aLRT/UFBoot. Named tips: the 19 characterized enzymes."))
save_fig(p1, "fig1_reference_tree", 9, 7.5)

# ---------------- Figure 3 ----------------
fl <- read_tsv(file.path(dat, "fig3_flows.tsv"), show_col_types = FALSE)
lev <- c("Clade I", "Clade II", "Clade III", "Clade IV", "Clade V", "unassigned")
TIERS <- c("validated", "extrapolated", "below floor")
fl <- fl %>% mutate(old_clade = factor(old_clade, levels = lev),
                    new_clade = factor(new_clade, levels = lev),
                    old_tier = factor(old_tier, levels = TIERS))
p3 <- ggplot(fl, aes(axis1 = old_clade, axis2 = new_clade, y = n)) +
  geom_alluvium(aes(fill = old_tier), width = 0.28, alpha = 0.75, knot.pos = 0.4,
                colour = SURFACE, linewidth = 0.2) +
  geom_stratum(width = 0.28, fill = "#f0efec", colour = INK2, linewidth = 0.3) +
  # big strata labelled inside; small ones (< 6,000) beside the bar, outward, so they cannot overprint
  geom_text(stat = "stratum", aes(label = after_stat(ifelse(count >= 6000, sprintf("%s\n%s", stratum,
            format(count, big.mark = ",")), ""))), size = 2.4, colour = INK, lineheight = 0.9) +
  geom_text(stat = "stratum", aes(label = after_stat(ifelse(count < 6000, sprintf("%s  %s", stratum,
            format(count, big.mark = ",")), "")),
            x = after_stat(ifelse(x == 1, x - 0.16, x + 0.16)),
            hjust = after_stat(ifelse(x == 1, 1, 0))), size = 2.3, colour = INK) +
  scale_x_discrete(limits = c("Nearest neighbour\n(19 enzymes)", "Placement\n(187-leaf tree)"),
                   expand = c(0.25, 0.25)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = setNames(SLOTS[1:3], TIERS), name = "Old call's tier") +
  labs(title = "What placement changed: 97,520 BFD PF01328 proteins",
       subtitle = paste("Validated-range calls hold; most extrapolated Clade III and V calls move",
                        "to IV or II, or leave Clades I-V"),
       x = NULL, y = "proteins") +
  theme_fig + theme(panel.grid.major.x = element_blank(), legend.position = "bottom")
save_fig(p3, "fig3_old_vs_placement", 7, 6)

# ---------------- Figure 4 ----------------
pr <- read_tsv(file.path(dat, "fig4_precision.tsv"), show_col_types = FALSE)
BINS <- c("70-100", "55-70", "45-55", "35-45", "<35")
CALLERS <- c("187-leaf reference (Zamocky 2022)", "19-reference (characterized only)")
ci <- function(k, n, which) {
  if (n == 0) return(NA_real_)
  binom.test(k, n)$conf.int[which]
}
pr <- pr %>% mutate(committed = correct + wrong,
                    precision = ifelse(committed > 0, correct / committed, NA),
                    lo = mapply(ci, correct, committed, 1), hi = mapply(ci, correct, committed, 2),
                    declined_frac = declined / tested,
                    bin = factor(bin, levels = BINS), caller = factor(caller, levels = CALLERS))
dodge <- position_dodge(width = 0.45)
caller_col <- setNames(SLOTS[1:2], CALLERS)
p4a <- ggplot(pr, aes(bin, precision, colour = caller, group = caller)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.6, position = dodge) +
  geom_point(size = 2.6, position = dodge, stroke = 0.6, fill = SURFACE, shape = 21) +
  geom_text(aes(label = sprintf("%d/%d", correct, committed), y = pmax(lo, 0) - 0.06),
            size = 2.3, position = dodge, colour = INK2, show.legend = FALSE) +
  scale_colour_manual(values = caller_col, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.1, 1.02), breaks = seq(0, 1, 0.25)) +
  labs(x = "identity to nearest surviving reference (%)", y = "precision when the caller commits",
       title = "Matched-identity validation on the 19 characterized enzymes",
       subtitle = paste0("Every reference above an identity ceiling removed before placing.\n",
                         "Bars: exact binomial 95 % CI; labels: correct / committed calls.")) +
  theme_fig + theme(legend.position = "top")
p4b <- ggplot(pr, aes(bin, declined_frac, colour = caller, group = caller)) +
  geom_point(size = 2.2, position = dodge) +
  scale_colour_manual(values = caller_col, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "identity to nearest surviving reference (%)", y = "declined to call",
       subtitle = "Abstentions: the 187-leaf caller declines where the 19-reference caller guesses") +
  theme_fig
save_fig(p4a / p4b + plot_layout(heights = c(3, 1.3)), "fig4_matched_identity", 6.5, 6.5)

message("wrote figures to ", out)
