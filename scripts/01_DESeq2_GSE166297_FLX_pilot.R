#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

counts_file <- "05_nfcore_results/GSE166297/salmon/salmon.merged.gene_counts.tsv"
meta_file   <- "02_metadata/GSE166297_pilot_metadata.csv"

out_de   <- "09_differential_expression/GSE166297"
out_expr <- "08_processed_expression/GSE166297_DESeq2"
out_fig  <- "14_figures/GSE166297_DESeq2"
out_log  <- "18_software_logs"

dir.create(out_de, recursive = TRUE, showWarnings = FALSE)
dir.create(out_expr, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

message("Reading nf-core Salmon gene counts...")
counts_raw <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)

message("Reading metadata...")
metadata <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)

needed_conditions <- c("Control_D10", "FLX_IC5", "FLX_IC20_100")
metadata <- metadata[metadata$condition_group %in% needed_conditions, , drop = FALSE]

# Detect sample column in metadata
candidate_cols <- c("nfcore_sample", "sample", "Run")
candidate_cols <- candidate_cols[candidate_cols %in% colnames(metadata)]

sample_col <- NULL
for (cc in candidate_cols) {
  if (all(metadata[[cc]] %in% colnames(counts_raw))) {
    sample_col <- cc
    break
  }
}

if (is.null(sample_col)) {
  message("Count matrix columns:")
  print(colnames(counts_raw))
  message("Metadata columns:")
  print(colnames(metadata))
  stop("Could not match metadata sample names with count matrix columns.")
}

metadata$sample_id <- metadata[[sample_col]]
message("Using sample column: ", sample_col)

# Identify gene annotation columns
gene_id_col <- intersect(c("gene_id", "gene", "Geneid", "ensembl_gene_id"), colnames(counts_raw))[1]
if (is.na(gene_id_col)) gene_id_col <- colnames(counts_raw)[1]

gene_name_col <- intersect(c("gene_name", "gene_symbol", "external_gene_name", "symbol"), colnames(counts_raw))[1]

gene_ids <- as.character(counts_raw[[gene_id_col]])
gene_ids_unique <- make.unique(gene_ids)

annotation <- data.frame(
  gene_id = gene_ids,
  stringsAsFactors = FALSE
)

if (!is.na(gene_name_col)) {
  annotation$gene_name <- counts_raw[[gene_name_col]]
} else {
  annotation$gene_name <- NA
}

rownames(annotation) <- gene_ids_unique

# Prepare count matrix
missing_samples <- setdiff(metadata$sample_id, colnames(counts_raw))
if (length(missing_samples) > 0) {
  stop("Missing samples in count matrix: ", paste(missing_samples, collapse = ", "))
}

metadata <- metadata[metadata$sample_id %in% colnames(counts_raw), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

count_mat <- as.matrix(counts_raw[, metadata$sample_id, drop = FALSE])
rownames(count_mat) <- gene_ids_unique
mode(count_mat) <- "numeric"

# Salmon tximport counts can be decimal; DESeq2 needs integer-like counts
count_mat <- round(count_mat)

# Remove all-zero genes
keep_nonzero <- rowSums(count_mat) > 0
count_mat <- count_mat[keep_nonzero, , drop = FALSE]
annotation <- annotation[rownames(count_mat), , drop = FALSE]

metadata$condition_group <- factor(
  metadata$condition_group,
  levels = c("Control_D10", "FLX_IC5", "FLX_IC20_100")
)

# Decide whether batch can be used safely
use_batch <- FALSE
if ("batch" %in% colnames(metadata)) {
  metadata$batch <- as.factor(metadata$batch)
  if (length(unique(metadata$batch)) > 1) {
    mm <- model.matrix(~ batch + condition_group, data = metadata)
    if (qr(mm)$rank == ncol(mm)) {
      use_batch <- TRUE
    }
  }
}

if (use_batch) {
  design_formula <- ~ batch + condition_group
  message("Using design: ~ batch + condition_group")
} else {
  design_formula <- ~ condition_group
  message("Using design: ~ condition_group")
}

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = design_formula
)

# Low-count filtering
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

message("Genes before filtering: ", nrow(count_mat))
message("Genes after filtering: ", nrow(dds))
message("Samples used: ", ncol(dds))

dds <- DESeq(dds)

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_out <- cbind(annotation[rownames(norm_counts), , drop = FALSE], norm_counts)

write.csv(
  norm_out,
  file.path(out_expr, "GSE166297_DESeq2_normalized_counts.csv"),
  row.names = FALSE
)

# VST matrix
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
vst_out <- cbind(annotation[rownames(vst_mat), , drop = FALSE], vst_mat)

write.csv(
  vst_out,
  file.path(out_expr, "GSE166297_DESeq2_vst_matrix.csv"),
  row.names = FALSE
)

# PCA
pcaData <- plotPCA(vsd, intgroup = "condition_group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, shape = condition_group)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.7, size = 2.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 12) +
  ggtitle("GSE166297 PCA: Control vs Fluoxetine")

ggsave(file.path(out_fig, "GSE166297_PCA_Control_FLX.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(out_fig, "GSE166297_PCA_Control_FLX.png"), p_pca, width = 7, height = 5, dpi = 300)

save_deseq_result <- function(dds, contrast_vec, prefix) {
  res <- results(dds, contrast = contrast_vec, alpha = 0.05)
  res <- res[order(res$padj), ]

  res_df <- as.data.frame(res)
  res_df$gene_id_unique <- rownames(res_df)

  res_df <- cbind(
    annotation[rownames(res_df), , drop = FALSE],
    res_df
  )

  res_df$significant_FDR05_log2FC058 <- with(
    res_df,
    ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 0.58, TRUE, FALSE)
  )

  write.csv(
    res_df,
    file.path(out_de, paste0(prefix, "_all_genes.csv")),
    row.names = FALSE
  )

  sig_df <- res_df[res_df$significant_FDR05_log2FC058 == TRUE, , drop = FALSE]

  write.csv(
    sig_df,
    file.path(out_de, paste0(prefix, "_significant_FDR05_log2FC058.csv")),
    row.names = FALSE
  )

  message(prefix, " total tested genes: ", nrow(res_df))
  message(prefix, " significant genes: ", nrow(sig_df))

  return(res_df)
}

res_ic5 <- save_deseq_result(
  dds,
  c("condition_group", "FLX_IC5", "Control_D10"),
  "FLX_IC5_vs_Control_D10"
)

res_ic20 <- save_deseq_result(
  dds,
  c("condition_group", "FLX_IC20_100", "Control_D10"),
  "FLX_IC20_100_vs_Control_D10"
)

plot_volcano <- function(res_df, title, prefix) {
  df <- res_df
  df$padj_plot <- ifelse(is.na(df$padj), 1, df$padj)
  df$neglog10padj <- -log10(pmax(df$padj_plot, .Machine$double.xmin))

  df$status <- "Not significant"
  df$status[!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange >= 0.58] <- "Up"
  df$status[!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange <= -0.58] <- "Down"

  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj, shape = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw(base_size = 12) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    ggtitle(title)

  ggsave(file.path(out_fig, paste0(prefix, "_volcano.pdf")), p, width = 7, height = 5)
  ggsave(file.path(out_fig, paste0(prefix, "_volcano.png")), p, width = 7, height = 5, dpi = 300)
}

plot_volcano(res_ic5, "Fluoxetine IC5 vs Control D10", "FLX_IC5_vs_Control_D10")
plot_volcano(res_ic20, "Fluoxetine IC20_100 vs Control D10", "FLX_IC20_100_vs_Control_D10")

write.csv(
  as.data.frame(colData(dds)),
  file.path(out_expr, "GSE166297_DESeq2_sample_metadata_used.csv"),
  row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(out_log, "GSE166297_DESeq2_sessionInfo.txt"))

message("DESeq2 analysis completed successfully.")
