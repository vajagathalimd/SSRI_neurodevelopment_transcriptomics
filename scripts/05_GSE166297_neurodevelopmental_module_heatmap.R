#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

vst_file  <- "08_processed_expression/GSE166297_DESeq2/GSE166297_DESeq2_vst_matrix.csv"
meta_file <- "02_metadata/GSE166297_pilot_metadata.csv"

fig_dir <- "14_figures/GSE166297_module_scores"
tab_dir <- "15_tables/GSE166297_module_scores"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

vst <- read.csv(vst_file, check.names = FALSE)
meta <- read.csv(meta_file, check.names = FALSE)

modules <- list(
  Neural_progenitor_cell_cycle = c("SOX2", "NES", "PAX6", "HES1", "HES5", "MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1"),
  Neuronal_differentiation = c("DCX", "TUBB3", "MAP2", "NEUROD1", "NEUROG2", "RBFOX3", "ELAVL3", "STMN2"),
  Axon_neurite_outgrowth = c("GAP43", "ROBO1", "ROBO2", "SLIT2", "SEMA3A", "DPYSL2", "NEFL", "NEFM", "NRCAM"),
  Synaptic_function = c("SYN1", "SYP", "SNAP25", "STXBP1", "DLG4", "NRXN1", "NLGN1", "NLGN3", "HOMER1"),
  Glutamatergic_signalling = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2A", "GRIN2B", "GRIA1", "GRIA2", "CAMK2A"),
  GABAergic_signalling = c("GAD1", "GAD2", "SLC6A1", "GABRA1", "GABRB2", "GABRG2", "DLX1", "DLX2"),
  Serotonergic_signalling = c("SLC6A4", "HTR1A", "HTR2A", "HTR2C", "HTR5A", "MAOA", "MAOB", "TPH1", "TPH2"),
  Astrocyte_glial_identity = c("GFAP", "S100B", "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "SOX9", "VIM"),
  Oligodendrocyte_myelin = c("OLIG1", "OLIG2", "SOX10", "MBP", "PLP1", "MAG", "MOG", "CNP"),
  Oxidative_stress_apoptosis = c("HMOX1", "NQO1", "SOD1", "SOD2", "DDIT3", "ATF4", "BAX", "BCL2", "CASP3"),
  Inflammatory_chemokine_signalling = c("IL6", "CXCL8", "CCL2", "NFKBIA", "STAT1", "IFITM3", "TNF", "IRF1")
)

# Detect gene label column
if ("gene_name" %in% colnames(vst) && any(!is.na(vst$gene_name) & vst$gene_name != "")) {
  vst$gene_label <- vst$gene_name
} else {
  vst$gene_label <- vst$gene_id
}

vst$gene_label <- toupper(vst$gene_label)

# Detect sample column in metadata
candidate_cols <- c("nfcore_sample", "sample", "Run")
candidate_cols <- candidate_cols[candidate_cols %in% colnames(meta)]

sample_col <- NULL
for (cc in candidate_cols) {
  if (any(meta[[cc]] %in% colnames(vst))) {
    sample_col <- cc
    break
  }
}

if (is.null(sample_col)) {
  stop("Could not match metadata samples to VST matrix columns.")
}

meta$sample_id <- meta[[sample_col]]
meta <- meta[meta$sample_id %in% colnames(vst), ]

sample_cols <- meta$sample_id

expr <- as.matrix(vst[, sample_cols, drop = FALSE])
rownames(expr) <- vst$gene_label
mode(expr) <- "numeric"

# Collapse duplicated gene symbols by mean
expr_df <- as.data.frame(expr)
expr_df$gene_label <- rownames(expr_df)

expr_collapsed <- aggregate(
  expr_df[, sample_cols, drop = FALSE],
  by = list(gene_label = expr_df$gene_label),
  FUN = mean
)

rownames(expr_collapsed) <- expr_collapsed$gene_label
expr_collapsed$gene_label <- NULL
expr <- as.matrix(expr_collapsed)

# Z-score each gene across samples
zexpr <- t(scale(t(expr)))
zexpr[is.na(zexpr)] <- 0

module_scores <- data.frame(sample_id = sample_cols)

matched_genes_long <- data.frame()

for (module_name in names(modules)) {
  genes <- toupper(modules[[module_name]])
  matched <- intersect(genes, rownames(zexpr))

  matched_genes_long <- rbind(
    matched_genes_long,
    data.frame(
      module = module_name,
      input_gene = genes,
      matched = genes %in% rownames(zexpr)
    )
  )

  if (length(matched) >= 2) {
    module_scores[[module_name]] <- colMeans(zexpr[matched, , drop = FALSE])
  } else if (length(matched) == 1) {
    module_scores[[module_name]] <- as.numeric(zexpr[matched, ])
  } else {
    module_scores[[module_name]] <- NA
  }
}

write.csv(
  matched_genes_long,
  file.path(tab_dir, "GSE166297_module_gene_matching.csv"),
  row.names = FALSE
)

module_scores <- merge(
  meta[, c("sample_id", "condition_group", "replicate_number", "batch"), drop = FALSE],
  module_scores,
  by = "sample_id"
)

write.csv(
  module_scores,
  file.path(tab_dir, "GSE166297_sample_level_neurodevelopmental_module_scores.csv"),
  row.names = FALSE
)

score_cols <- setdiff(
  colnames(module_scores),
  c("sample_id", "condition_group", "replicate_number", "batch")
)

condition_scores <- aggregate(
  module_scores[, score_cols, drop = FALSE],
  by = list(condition_group = module_scores$condition_group),
  FUN = mean,
  na.rm = TRUE
)

write.csv(
  condition_scores,
  file.path(tab_dir, "GSE166297_condition_mean_neurodevelopmental_module_scores.csv"),
  row.names = FALSE
)

# Long format for sample-level heatmap
long_sample <- data.frame()

for (module_name in score_cols) {
  long_sample <- rbind(
    long_sample,
    data.frame(
      sample_id = module_scores$sample_id,
      condition_group = module_scores$condition_group,
      module = module_name,
      score = module_scores[[module_name]]
    )
  )
}

long_sample$condition_group <- factor(
  long_sample$condition_group,
  levels = c("Control_D10", "FLX_IC5", "FLX_IC20_100")
)

p1 <- ggplot(long_sample, aes(x = sample_id, y = module, fill = score)) +
  geom_tile() +
  facet_grid(. ~ condition_group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_rect(fill = "white")
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Sample-level neurodevelopmental module scores")

ggsave(
  file.path(fig_dir, "GSE166297_sample_level_neurodevelopmental_module_heatmap.pdf"),
  p1,
  width = 12,
  height = 6
)

ggsave(
  file.path(fig_dir, "GSE166297_sample_level_neurodevelopmental_module_heatmap.png"),
  p1,
  width = 12,
  height = 6,
  dpi = 300
)

# Long format for condition mean heatmap
long_condition <- data.frame()

for (module_name in score_cols) {
  long_condition <- rbind(
    long_condition,
    data.frame(
      condition_group = condition_scores$condition_group,
      module = module_name,
      score = condition_scores[[module_name]]
    )
  )
}

long_condition$condition_group <- factor(
  long_condition$condition_group,
  levels = c("Control_D10", "FLX_IC5", "FLX_IC20_100")
)

p2 <- ggplot(long_condition, aes(x = condition_group, y = module, fill = score)) +
  geom_tile() +
  geom_text(aes(label = round(score, 2)), size = 3) +
  theme_bw(base_size = 11) +
  xlab("") +
  ylab("") +
  ggtitle("Condition-mean neurodevelopmental module scores")

ggsave(
  file.path(fig_dir, "GSE166297_condition_mean_neurodevelopmental_module_heatmap.pdf"),
  p2,
  width = 7,
  height = 6
)

ggsave(
  file.path(fig_dir, "GSE166297_condition_mean_neurodevelopmental_module_heatmap.png"),
  p2,
  width = 7,
  height = 6,
  dpi = 300
)

message("Neurodevelopmental module heatmap completed.")
message("Sample column used: ", sample_col)
message("Modules analysed: ", length(score_cols))
