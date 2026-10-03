#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

out_dir <- "14_figures/manuscript_combined_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

draw_box <- function(label, x, y, w = 0.22, h = 0.105, fontsize = 10) {
  grid.roundrect(
    x = x, y = y, width = w, height = h,
    r = unit(0.02, "npc"),
    gp = gpar(fill = "grey95", col = "black", lwd = 1.2)
  )
  grid.text(label, x = x, y = y, gp = gpar(fontsize = fontsize))
}

draw_arrow <- function(x1, y1, x2, y2) {
  grid.lines(
    x = unit(c(x1, x2), "npc"),
    y = unit(c(y1, y2), "npc"),
    arrow = arrow(type = "closed", length = unit(0.018, "npc")),
    gp = gpar(lwd = 1.1)
  )
}

plot_workflow <- function() {
  grid.newpage()

  grid.text(
    "GSE166297 fluoxetine transcriptomic reanalysis workflow",
    x = 0.5, y = 0.95,
    gp = gpar(fontsize = 16, fontface = "bold")
  )

  draw_box("Public RNA-seq dataset\nGSE166297", 0.16, 0.78)
  draw_box("Sample selection\nControl_D10, FLX_IC5,\nFLX_IC20_100\nn = 8/group", 0.50, 0.78)
  draw_box("FASTQ processing\nnf-core/rnaseq\nfastp + FastQC", 0.84, 0.78)

  draw_arrow(0.27, 0.78, 0.39, 0.78)
  draw_arrow(0.61, 0.78, 0.73, 0.78)

  draw_box("Transcript quantification\nSalmon pseudoalignment\nGRCh38 reference", 0.16, 0.55)
  draw_box("Gene-level count matrix\nSalmon/tximport output", 0.50, 0.55)
  draw_box("Differential expression\nDESeq2\n~ batch + condition", 0.84, 0.55)

  draw_arrow(0.84, 0.725, 0.16, 0.61)
  draw_arrow(0.27, 0.55, 0.39, 0.55)
  draw_arrow(0.61, 0.55, 0.73, 0.55)

  draw_box("DEG filtering\nFDR < 0.05\n|log2FC| ≥ 0.58", 0.16, 0.32)
  draw_box("DEG overlap\nshared and condition-specific\nfluoxetine-responsive genes", 0.50, 0.32)
  draw_box("Functional enrichment\ng:Profiler\nGO, Reactome, KEGG, WP", 0.84, 0.32)

  draw_arrow(0.84, 0.495, 0.16, 0.38)
  draw_arrow(0.27, 0.32, 0.39, 0.32)
  draw_arrow(0.61, 0.32, 0.73, 0.32)

  draw_box("Curated neurodevelopmental\nmodule scoring\nVST expression + z-scores", 0.32, 0.11, w = 0.28)
  draw_box("Manuscript figures\nDEGs, enrichment,\nmodule shifts", 0.68, 0.11, w = 0.28)

  draw_arrow(0.84, 0.265, 0.68, 0.17)
  draw_arrow(0.50, 0.265, 0.32, 0.17)
  draw_arrow(0.46, 0.11, 0.54, 0.11)
}

pdf(file.path(out_dir, "Figure_1_GSE166297_workflow.pdf"), width = 11, height = 7)
plot_workflow()
dev.off()

png(file.path(out_dir, "Figure_1_GSE166297_workflow.png"), width = 11, height = 7, units = "in", res = 300)
plot_workflow()
dev.off()

message("Figure 1 workflow created.")
