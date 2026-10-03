#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
  library(gridExtra)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

in_file <- "15_tables/cross_SSRI_gene_level/SSRI_top_pair_shared_downregulated_DEGs.csv"

out_tab <- "15_tables/cross_SSRI_gene_level"
out_fig <- "14_figures/cross_SSRI_gene_level"
out_combined <- "14_figures/manuscript_combined_figures"

dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_combined, recursive = TRUE, showWarnings = FALSE)

shared <- read.csv(in_file, check.names = FALSE)

shared$functional_category <- "Other shared downregulated genes"

shared$functional_category[shared$gene_key %in% c(
  "BDNF", "GABRA2", "THRB", "SEMA3D", "LGI1", "CDH13", "SH2D5", "SCN4A"
)] <- "Neurodevelopment, synaptic or excitability-related"

shared$functional_category[shared$gene_key %in% c(
  "ATF5", "TRIB3", "ASNS", "SESN2", "SLC1A5", "ASS1", "PHGDH", "GPT2", "CHAC1", "PCK2"
)] <- "Amino-acid, metabolic or stress-response related"

shared$functional_category[shared$gene_key %in% c(
  "CEBPB", "JDP2", "STC2", "TNFRSF12A", "IL11", "CRYAB", "SERPINB9", "EMP1"
)] <- "Cellular stress, inflammatory or growth-response related"

unique_genes <- unique(shared$gene_key)

summary_list <- lapply(unique_genes, function(g) {
  sub <- shared[shared$gene_key == g, ]

  data.frame(
    gene_key = g,
    functional_category = unique(sub$functional_category)[1],
    n_shared_pair_rows = nrow(sub),
    n_unique_pairs = length(unique(sub$pair)),
    pairs = paste(unique(sub$pair), collapse = "; "),
    mean_logFC_FLX = mean(sub$logFC_FLX, na.rm = TRUE),
    mean_logFC_CIT = mean(sub$logFC_CIT, na.rm = TRUE),
    min_padj_FLX = min(sub$padj_FLX, na.rm = TRUE),
    min_padj_CIT = min(sub$padj_CIT, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

gene_summary <- do.call(rbind, summary_list)

gene_summary <- gene_summary[
  order(-gene_summary$n_unique_pairs, gene_summary$mean_logFC_FLX, gene_summary$mean_logFC_CIT),
]

write.csv(
  gene_summary,
  file.path(out_tab, "SSRI_shared_downregulated_gene_recurrence_summary.csv"),
  row.names = FALSE
)

top_plot <- head(gene_summary, 25)
top_plot$gene_key <- factor(top_plot$gene_key, levels = rev(top_plot$gene_key))

p_gene <- ggplot(top_plot, aes(x = n_unique_pairs, y = gene_key, shape = functional_category)) +
  geom_point(size = 3) +
  theme_bw(base_size = 11) +
  xlab("Number of top SSRI contrast-pairs sharing downregulation") +
  ylab("") +
  ggtitle("Recurrent shared downregulated genes across fluoxetine and citalopram contrasts") +
  theme(legend.position = "bottom")

ggsave(
  file.path(out_fig, "SSRI_shared_downregulated_gene_recurrence_barplot.pdf"),
  p_gene,
  width = 9,
  height = 7
)

ggsave(
  file.path(out_fig, "SSRI_shared_downregulated_gene_recurrence_barplot.png"),
  p_gene,
  width = 9,
  height = 7,
  dpi = 300
)

# Helper for combined figure panels
make_panel <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing figure: ", path)
  }

  img <- png::readPNG(path)
  g_img <- grid::rasterGrob(img, interpolate = TRUE)

  g_label <- grid::textGrob(
    label,
    x = unit(0.02, "npc"),
    y = unit(0.98, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontface = "bold", fontsize = 18)
  )

  grid::grobTree(g_img, g_label)
}

fig5 <- list(
  make_panel("14_figures/cross_SSRI_integration/SSRI_combined_DEG_count_summary.png", "A"),
  make_panel("14_figures/cross_SSRI_integration/SSRI_focused_D10_D13_module_shift_heatmap.png", "B"),
  make_panel("14_figures/cross_SSRI_gene_level/SSRI_gene_level_spearman_correlation_heatmap.png", "C"),
  make_panel("14_figures/cross_SSRI_gene_level/SSRI_shared_significant_DEG_overlap_barplot.png", "D"),
  make_panel("14_figures/cross_SSRI_gene_level/SSRI_shared_downregulated_gene_recurrence_barplot.png", "E")
)

pdf(
  file.path(out_combined, "Figure_5_cross_SSRI_comparative_integration.pdf"),
  width = 12,
  height = 22
)
gridExtra::grid.arrange(grobs = fig5, ncol = 1)
dev.off()

png(
  file.path(out_combined, "Figure_5_cross_SSRI_comparative_integration.png"),
  width = 12,
  height = 22,
  units = "in",
  res = 300
)
gridExtra::grid.arrange(grobs = fig5, ncol = 1)
dev.off()

message("Final cross-SSRI gene summary and Figure 5 created.")
message("Unique shared downregulated genes: ", nrow(gene_summary))
message("Top recurrent genes:")
print(head(gene_summary[, c("gene_key", "functional_category", "n_unique_pairs", "mean_logFC_FLX", "mean_logFC_CIT")], 20))
