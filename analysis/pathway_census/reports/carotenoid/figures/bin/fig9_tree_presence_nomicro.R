# NOTE: copy of fig9_tree_presence.R for the Microsporidia-excluded, Blastocladiomycota-rooted primary run (adds the exclusion and root to the subtitle).
#!/usr/bin/env Rscript
# Figure 9: genus tree with tip presence of each step (rings) and labelled clades.
# Usage: Rscript fig9_tree_presence.R <results/gainloss/carotenoid dir> <out dir>
suppressPackageStartupMessages({library(ape); library(ggplot2); library(ggtree)})
a <- commandArgs(TRUE); rd <- a[1]; od <- a[2]
small_tree <- if (length(a) >= 3) strsplit(a[3], ",")[[1]] else c("phytoene_desaturase", "cart")  # placed on a small reference tree
stem <- if (length(a) >= 4) a[4] else "fig9_tree_presence"
prim <- readRDS(file.path(rd, "primary.rds")); cal <- read.delim(file.path(rd, "calls.tsv"), stringsAsFactors = FALSE)
tree <- prim$tree; steps <- prim$steps; src <- prim$source
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; ABSENT <- "#dcdbd6"
# ring data: one column per step, rows = tips
lab <- function(s) paste0(gsub("_", " ", s), " (", if (src[[s]] == "hits") "search hits" else if (s %in% small_tree) "placed, small tree" else "placed, with named members", ")")
df <- data.frame(row.names = tree$tip.label)
for (s in steps) {
  v <- setNames(cal[[s]], cal$locustag)[tree$tip.label]
  df[[lab(s)]] <- factor(ifelse(is.na(v), "no data", ifelse(v == 1, "present", "absent")), levels = c("present", "absent", "no data"))
}
# clade labels: for each taxon, the internal node whose clade holds the most of its tips at purity >= 0.9
pp <- prop.part(tree, check.labels = TRUE); tl <- attr(pp, "labels"); n <- Ntip(tree)
best_clade <- function(tips) {
  best <- NULL; bn <- 0
  for (i in seq_along(pp)) { d <- tl[pp[[i]]]; k <- sum(d %in% tips)
    if (k >= 5 && k / length(d) >= 0.9 && k > bn) { bn <- k; best <- n + i } }
  list(node = best, n_in = bn, n_clade = if (is.null(best)) NA else length(pp[[best - n]]), n_taxon = length(tips))
}
cls <- setNames(cal$CLASS, cal$locustag)[tree$tip.label]; phy <- setNames(cal$PHYLUM, cal$locustag)[tree$tip.label]
want <- list()
for (g in setdiff(names(which(table(cls) >= 25)), c("", "Mucoromycetes", "Chytridiomycetes"))) want[[g]] <- tree$tip.label[which(cls == g)]
for (g in c("Mucoromycota", "Chytridiomycota", "Microsporidia", "Zoopagomycota", "Kickxellomycota", "Mortierellomycota", "Glomeromycota"))
  want[[g]] <- tree$tip.label[which(phy == g)]
for (g in c("Basidiomycota", "Ascomycota")) want[[g]] <- tree$tip.label[which(phy == g)]
cl <- lapply(names(want), function(g) { b <- best_clade(want[[g]]); data.frame(label = g, node = if (is.null(b$node)) NA else b$node, tips_of_taxon_in_tree = b$n_taxon,
  tips_of_taxon_in_labelled_clade = b$n_in, tips_in_labelled_clade = b$n_clade) })
cl <- do.call(rbind, cl); cl <- cl[!is.na(cl$node), ]
write.table(cl, file.path(rd, "fig9_clade_labels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cl$txt <- paste0(cl$label, " (", cl$tips_of_taxon_in_labelled_clade, ")")

ntip <- n; gap <- ntip * 24 / (360 - 24)
p <- ggtree(tree, layout = "fan", open.angle = 24, size = 0.1, color = "#6b6a66", branch.length = "none") +
  theme(plot.background = element_rect(fill = SURFACE, color = NA), panel.background = element_rect(fill = SURFACE, color = NA),
        legend.position = "none")
depth <- max(p$data$x)          # tree height in cladogram units
ringw <- depth * 0.16           # width of one ring
p2 <- p
colw <- ringw * 0.9
for (k in seq_along(steps)) {
  p2 <- gheatmap(p2, df[, k, drop = FALSE], offset = depth * 0.04 + (k - 1) * ringw, width = colw / depth, colnames = FALSE, color = NA)
}
p2 <- p2 + scale_fill_manual(values = c(present = BLUE, absent = ABSENT, `no data` = ORANGE), name = NULL) +
  ggtitle(paste0("Presence of carotenoid steps on the genus tree (", prettyNum(n, big.mark = ","), " tips; cladogram, branch lengths not drawn)"),
          subtitle = paste0("Rings, inside to outside:\n", paste0(seq_along(steps)[1:2], " = ", names(df)[1:2], collapse = "; "), "\n", paste0(seq_along(steps)[3:4], " = ", names(df)[3:4], collapse = "; ")))
# ring numbers in the gap of the fan
gx <- depth + depth * 0.04 + ringw * (seq_along(steps) - 1) + colw / 2
ring_lab <- data.frame(x = gx, y = ntip + gap / 2, lab = as.character(seq_along(steps)))
p2 <- p2 + geom_text(data = ring_lab, aes(x = x, y = y, label = lab), inherit.aes = FALSE, size = 3.2, colour = INK, fontface = "bold")
# clade bars: classes on the inner tier, phyla on the outer tier; small phyla alternate between two radii
base <- depth + depth * 0.04 + 4 * ringw + depth * 0.03
cl$tier <- ifelse(cl$label %in% c("Basidiomycota", "Ascomycota"), 2, 1); cl$stag <- 0
cl$first <- vapply(cl$node, function(nd) min(which(tree$tip.label %in% tl[pp[[nd - n]]])), 1)
small <- cl$tier == 1 & cl$tips_in_labelled_clade < 100
cl <- cl[order(cl$first), ]; small <- cl$tier == 1 & cl$tips_in_labelled_clade < 100
cl$stag[small] <- ((seq_len(sum(small)) - 1) %% 3) * depth * 0.18
early <- cl$label %in% c("Mucoromycota", "Chytridiomycota", "Microsporidia", "Zoopagomycota", "Kickxellomycota", "Mortierellomycota", "Glomeromycota")
cl$stag[early] <- 0
for (i in seq_len(nrow(cl))) p2 <- if (early[i]) p2 + geom_cladelab(node = cl$node[i], label = cl$txt[i], angle = 0, hjust = 0,
    offset = base - depth + (cl$tier[i] - 1) * depth * 0.5, offset.text = depth * 0.06, barsize = 0.6, fontsize = 2.3,
    textcolor = INK, barcolor = INK, horizontal = FALSE, align = TRUE) else p2 + geom_cladelab(node = cl$node[i], label = cl$txt[i], angle = "auto",
    offset = base - depth + (cl$tier[i] - 1) * depth * 0.5 + cl$stag[i], offset.text = depth * 0.03, barsize = 0.6, fontsize = 2.3,
    textcolor = INK, barcolor = INK, horizontal = FALSE, align = TRUE)
p2 <- p2 + labs(subtitle = paste0("Microsporidia excluded (30 tips); tree rooted on Blastocladiomycota.\n", p2$labels$subtitle))
ggsave(file.path(od, paste0(stem, ".png")), p2 + theme(legend.position = "bottom"), width = 10, height = 10.5, dpi = 150, bg = SURFACE)
ggsave(file.path(od, paste0(stem, ".pdf")), p2 + theme(legend.position = "bottom"), width = 10, height = 10.5, bg = SURFACE)
cat(sprintf("ring order, inside to outside: %s\n", paste(names(df), collapse = " | ")))
