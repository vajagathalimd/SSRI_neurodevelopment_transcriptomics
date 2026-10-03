#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(gridExtra)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

out_dir <- "14_figures/manuscript_combined_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_panel <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing figure file: ", path)
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

save_combined <- function(grobs, filename, ncol, width, height) {
  pdf(file.path(out_dir, paste0(filename, ".pdf")), width = width, height = height)
  gridExtra::grid.arrange(grobs = grobs, ncol = ncol)
  dev.off()

  png(file.path(out_dir, paste0(filename, ".png")), width = width, height = height, units = "in", res = 300)
  gridExtra::grid.arrange(grobs = grobs, ncol = ncol)
  dev.off()
}

# Figure 2: DEG overview
fig2 <- list(
  make_panel("14_figures/GSE166297_DESeq2/GSE166297_PCA_Control_FLX.png", "A"),
  make_panel("14_figures/GSE166297_DESeq2/FLX_IC5_vs_Control_D10_volcano.png", "B"),
  make_panel("14_figures/GSE166297_DESeq2/FLX_IC20_100_vs_Control_D10_volcano.png", "C"),
  make_panel("14_figures/GSE166297_DESeq2/GSE166297_DEG_count_barplot.png", "D"),
  make_panel("14_figures/GSE166297_DESeq2/GSE166297_DEG_overlap_barplot.png", "E")
)

save_combined(
  fig2,
  "Figure_2_GSE166297_DEG_overview",
  ncol = 2,
  width = 12,
  height = 14
)

# Figure 3: Enrichment overview
fig3 <- list(
  make_panel("14_figures/GSE166297_enrichment/FLX_IC5_all_significant_DEGs_top20_enrichment.png", "A"),
  make_panel("14_figures/GSE166297_enrichment/FLX_IC20_100_all_significant_DEGs_top20_enrichment.png", "B"),
  make_panel("14_figures/GSE166297_enrichment/Shared_FLX_IC5_and_IC20_100_DEGs_top20_enrichment.png", "C"),
  make_panel("14_figures/GSE166297_enrichment/FLX_IC5_downregulated_DEGs_top20_enrichment.png", "D"),
  make_panel("14_figures/GSE166297_enrichment/FLX_IC20_100_downregulated_DEGs_top20_enrichment.png", "E")
)

save_combined(
  fig3,
  "Figure_3_GSE166297_enrichment_overview",
  ncol = 1,
  width = 10,
  height = 22
)

# Figure 4: Neurodevelopmental module overview
fig4 <- list(
  make_panel("14_figures/GSE166297_module_scores/GSE166297_condition_mean_neurodevelopmental_module_heatmap.png", "A"),
  make_panel("14_figures/GSE166297_module_stats/GSE166297_module_mean_difference_vs_control.png", "B")
)

save_combined(
  fig4,
  "Figure_4_GSE166297_neurodevelopmental_module_shifts",
  ncol = 1,
  width = 9,
  height = 12
)

# Supplementary module sample-level heatmap
figS <- list(
  make_panel("14_figures/GSE166297_module_scores/GSE166297_sample_level_neurodevelopmental_module_heatmap.png", "A")
)

save_combined(
  figS,
  "Supplementary_Figure_GSE166297_sample_level_module_heatmap",
  ncol = 1,
  width = 12,
  height = 7
)

message("Combined manuscript figures created in: ", out_dir)
