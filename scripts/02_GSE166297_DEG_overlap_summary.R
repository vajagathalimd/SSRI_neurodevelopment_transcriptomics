#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

de_dir  <- "09_differential_expression/GSE166297"
fig_dir <- "14_figures/GSE166297_DESeq2"
tab_dir <- "15_tables/GSE166297_DESeq2"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

ic5_all <- read.csv(file.path(de_dir, "FLX_IC5_vs_Control_D10_all_genes.csv"), check.names = FALSE)
ic20_all <- read.csv(file.path(de_dir, "FLX_IC20_100_vs_Control_D10_all_genes.csv"), check.names = FALSE)

ic5_sig <- read.csv(file.path(de_dir, "FLX_IC5_vs_Control_D10_significant_FDR05_log2FC058.csv"), check.names = FALSE)
ic20_sig <- read.csv(file.path(de_dir, "FLX_IC20_100_vs_Control_D10_significant_FDR05_log2FC058.csv"), check.names = FALSE)

get_gene_label <- function(df) {
  if ("gene_name" %in% colnames(df) && any(!is.na(df$gene_name) & df$gene_name != "")) {
    return(ifelse(is.na(df$gene_name) | df$gene_name == "", df$gene_id, df$gene_name))
  } else {
    return(df$gene_id)
  }
}

ic5_sig$gene_label <- get_gene_label(ic5_sig)
ic20_sig$gene_label <- get_gene_label(ic20_sig)

ic5_genes <- unique(ic5_sig$gene_label)
ic20_genes <- unique(ic20_sig$gene_label)

overlap_genes <- intersect(ic5_genes, ic20_genes)
ic5_only <- setdiff(ic5_genes, ic20_genes)
ic20_only <- setdiff(ic20_genes, ic5_genes)

summary_df <- data.frame(
  category = c("FLX_IC5 significant", "FLX_IC20_100 significant", "Shared", "FLX_IC5 only", "FLX_IC20_100 only"),
  count = c(length(ic5_genes), length(ic20_genes), length(overlap_genes), length(ic5_only), length(ic20_only))
)

write.csv(summary_df, file.path(tab_dir, "GSE166297_DEG_overlap_summary.csv"), row.names = FALSE)

write.csv(data.frame(shared_DEGs = overlap_genes),
          file.path(tab_dir, "GSE166297_shared_DEGs_FLX_IC5_and_IC20_100.csv"),
          row.names = FALSE)

write.csv(data.frame(FLX_IC5_only_DEGs = ic5_only),
          file.path(tab_dir, "GSE166297_FLX_IC5_only_DEGs.csv"),
          row.names = FALSE)

write.csv(data.frame(FLX_IC20_100_only_DEGs = ic20_only),
          file.path(tab_dir, "GSE166297_FLX_IC20_100_only_DEGs.csv"),
          row.names = FALSE)

# Direction counts
direction_summary <- function(df, contrast_name) {
  data.frame(
    contrast = contrast_name,
    upregulated = sum(df$log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(df$log2FoldChange < 0, na.rm = TRUE),
    total = nrow(df)
  )
}

direction_df <- rbind(
  direction_summary(ic5_sig, "FLX_IC5_vs_Control_D10"),
  direction_summary(ic20_sig, "FLX_IC20_100_vs_Control_D10")
)

write.csv(direction_df, file.path(tab_dir, "GSE166297_DEG_direction_summary.csv"), row.names = FALSE)

# Top genes
top_table <- function(df, contrast_name) {
  df$contrast <- contrast_name
  df$gene_label <- get_gene_label(df)

  top_up <- df[order(-df$log2FoldChange, df$padj), ]
  top_down <- df[order(df$log2FoldChange, df$padj), ]

  out <- rbind(
    head(top_up, 25),
    head(top_down, 25)
  )

  keep_cols <- intersect(c("contrast", "gene_id", "gene_name", "gene_label", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), colnames(out))
  out[, keep_cols, drop = FALSE]
}

top_genes <- rbind(
  top_table(ic5_sig, "FLX_IC5_vs_Control_D10"),
  top_table(ic20_sig, "FLX_IC20_100_vs_Control_D10")
)

write.csv(top_genes, file.path(tab_dir, "GSE166297_top_25_up_down_DEGs_each_contrast.csv"), row.names = FALSE)

# Bar plot: DEG counts
p1 <- ggplot(direction_df, aes(x = contrast, y = total)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.4, size = 4) +
  theme_bw(base_size = 12) +
  xlab("") +
  ylab("Number of significant DEGs") +
  ggtitle("Significant DEGs after fluoxetine exposure") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(fig_dir, "GSE166297_DEG_count_barplot.pdf"), p1, width = 6, height = 4)
ggsave(file.path(fig_dir, "GSE166297_DEG_count_barplot.png"), p1, width = 6, height = 4, dpi = 300)

# Bar plot: overlap categories
p2 <- ggplot(summary_df[3:5, ], aes(x = category, y = count)) +
  geom_col() +
  geom_text(aes(label = count), vjust = -0.4, size = 4) +
  theme_bw(base_size = 12) +
  xlab("") +
  ylab("Number of DEGs") +
  ggtitle("Overlap of fluoxetine-responsive DEGs") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(fig_dir, "GSE166297_DEG_overlap_barplot.pdf"), p2, width = 6, height = 4)
ggsave(file.path(fig_dir, "GSE166297_DEG_overlap_barplot.png"), p2, width = 6, height = 4, dpi = 300)

message("Overlap and DEG summary completed.")
message("Shared DEGs: ", length(overlap_genes))
message("FLX_IC5 only: ", length(ic5_only))
message("FLX_IC20_100 only: ", length(ic20_only))
