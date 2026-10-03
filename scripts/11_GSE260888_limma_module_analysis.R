#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

cpm_file <- "03_raw_downloads/GSE260888_processed/GSE260888/GSE260888_cpm.csv.gz"

out_expr <- "08_processed_expression/GSE260888"
out_de   <- "09_differential_expression/GSE260888"
out_fig  <- "14_figures/GSE260888"
out_tab  <- "15_tables/GSE260888"
out_log  <- "18_software_logs/GSE260888"

dir.create(out_expr, recursive = TRUE, showWarnings = FALSE)
dir.create(out_de, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

message("Reading GSE260888 CPM matrix...")
cpm <- read.csv(gzfile(cpm_file), check.names = FALSE, stringsAsFactors = FALSE)

gene_id <- cpm[[1]]
expr <- as.matrix(cpm[, -1, drop = FALSE])
rownames(expr) <- gene_id
mode(expr) <- "numeric"

# If values look like raw CPM, log2-transform; if already log-scale, keep as is.
max_val <- max(expr, na.rm = TRUE)
if (max_val > 100) {
  message("Expression values appear non-log CPM. Applying log2(CPM + 1).")
  expr <- log2(expr + 1)
} else {
  message("Expression values appear log-scaled CPM. Using as provided.")
}

sample_ids <- colnames(expr)

parse_sample <- function(x) {
  # Expected pattern: 257-D10-CTRL-A or 260-D10-CIT50-A
  parts <- strsplit(x, "-", fixed = TRUE)[[1]]

  if (length(parts) != 4) {
    return(data.frame(
      sample_id = x,
      sample_number = NA,
      day = NA,
      day_group = NA,
      treatment = NA,
      dose_nM = NA,
      replicate = NA,
      stringsAsFactors = FALSE
    ))
  }

  sample_number <- parts[1]
  day_group <- parts[2]
  treatment <- parts[3]
  replicate <- parts[4]

  # Header has one likely typo: CIT40 in D13 should be CIT50
  if (treatment == "CIT40") treatment <- "CIT50"

  day <- as.numeric(gsub("D", "", day_group))

  dose_nM <- ifelse(
    treatment == "CTRL",
    0,
    as.numeric(gsub("CIT", "", treatment))
  )

  data.frame(
    sample_id = x,
    sample_number = sample_number,
    day = day,
    day_group = paste0("D", day),
    treatment = treatment,
    dose_nM = dose_nM,
    replicate = replicate,
    stringsAsFactors = FALSE
  )
}

sample_meta <- do.call(rbind, lapply(sample_ids, parse_sample))

sample_meta$treatment <- factor(
  sample_meta$treatment,
  levels = c("CTRL", "CIT50", "CIT100", "CIT200", "CIT400")
)

sample_meta$day_group <- factor(
  sample_meta$day_group,
  levels = c("D0", "D6", "D10", "D13")
)

write.csv(
  sample_meta,
  file.path(out_expr, "GSE260888_parsed_sample_metadata.csv"),
  row.names = FALSE
)

sample_summary <- as.data.frame.matrix(table(sample_meta$day_group, sample_meta$treatment))
sample_summary$day_group <- rownames(sample_summary)
sample_summary <- sample_summary[, c("day_group", setdiff(colnames(sample_summary), "day_group"))]

write.csv(
  sample_summary,
  file.path(out_tab, "GSE260888_sample_count_summary.csv"),
  row.names = FALSE
)

message("Sample summary:")
print(sample_summary)

# Gene annotation
ensg_clean <- sub("\\..*$", "", rownames(expr))

message("Mapping Ensembl IDs to gene symbols...")
gene_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensg_clean),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

symbol_vec <- gene_symbol[ensg_clean]

gene_annotation <- data.frame(
  gene_id = rownames(expr),
  ensembl_id = ensg_clean,
  gene_symbol = as.character(symbol_vec),
  stringsAsFactors = FALSE
)

rownames(gene_annotation) <- rownames(expr)

write.csv(
  gene_annotation,
  file.path(out_expr, "GSE260888_gene_annotation_Ensembl_to_symbol.csv"),
  row.names = FALSE
)

expr_out <- cbind(gene_annotation, expr)

write.csv(
  expr_out,
  file.path(out_expr, "GSE260888_logCPM_expression_with_annotation.csv"),
  row.names = FALSE
)

# PCA using top variable genes
message("Creating PCA...")
vars <- apply(expr, 1, var, na.rm = TRUE)
top_genes <- names(sort(vars, decreasing = TRUE))[1:min(5000, length(vars))]
pca <- prcomp(t(expr[top_genes, , drop = FALSE]), scale. = TRUE)

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  sample_meta[match(rownames(pca$x), sample_meta$sample_id), ],
  stringsAsFactors = FALSE
)

percent <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = treatment)) +
  geom_point(size = 2.8) +
  facet_wrap(~ day_group, scales = "free") +
  theme_bw(base_size = 12) +
  xlab(paste0("PC1: ", percent[1], "% variance")) +
  ylab(paste0("PC2: ", percent[2], "% variance")) +
  ggtitle("GSE260888 PCA across citalopram neuronal differentiation samples")

ggsave(file.path(out_fig, "GSE260888_PCA_all_samples.pdf"), p_pca, width = 10, height = 6)
ggsave(file.path(out_fig, "GSE260888_PCA_all_samples.png"), p_pca, width = 10, height = 6, dpi = 300)

# Limma DE per day and dose
message("Running limma differential expression...")
days_to_test <- c("D6", "D10", "D13")
doses <- c("CIT50", "CIT100", "CIT200", "CIT400")

de_summary <- data.frame()

for (day_i in days_to_test) {
  meta_day <- sample_meta[sample_meta$day_group == day_i & sample_meta$treatment %in% c("CTRL", doses), ]
  meta_day <- meta_day[!is.na(meta_day$treatment), ]

  present_treatments <- as.character(unique(meta_day$treatment))
  present_treatments <- present_treatments[!is.na(present_treatments)]

  if (!("CTRL" %in% present_treatments)) {
    message("Skipping ", day_i, ": no CTRL.")
    next
  }

  expr_day <- expr[, meta_day$sample_id, drop = FALSE]

  meta_day$treatment <- factor(
    as.character(meta_day$treatment),
    levels = c("CTRL", doses)
  )
  meta_day$treatment <- droplevels(meta_day$treatment)

  design <- model.matrix(~ 0 + treatment, data = meta_day)
  colnames(design) <- gsub("^treatment", "", colnames(design))

  fit <- lmFit(expr_day, design)

  available_doses <- intersect(doses, colnames(design))

  for (dose_i in available_doses) {
    contrast_name <- paste0(day_i, "_", dose_i, "_vs_CTRL")

    contrast_matrix <- makeContrasts(
      contrasts = paste0(dose_i, "-CTRL"),
      levels = design
    )

    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    res <- topTable(fit2, number = Inf, sort.by = "P")
    res$gene_id <- rownames(res)
    res$ensembl_id <- sub("\\..*$", "", res$gene_id)
    res$gene_symbol <- gene_annotation[res$gene_id, "gene_symbol"]

    res <- res[, c("gene_id", "ensembl_id", "gene_symbol", setdiff(colnames(res), c("gene_id", "ensembl_id", "gene_symbol")))]

    res$significant_FDR05_log2FC058 <- with(
      res,
      ifelse(!is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) >= 0.58, TRUE, FALSE)
    )

    all_file <- file.path(out_de, paste0(contrast_name, "_all_genes.csv"))
    sig_file <- file.path(out_de, paste0(contrast_name, "_significant_FDR05_log2FC058.csv"))

    write.csv(res, all_file, row.names = FALSE)

    sig <- res[res$significant_FDR05_log2FC058 == TRUE, , drop = FALSE]
    write.csv(sig, sig_file, row.names = FALSE)

    n_up <- sum(sig$logFC > 0, na.rm = TRUE)
    n_down <- sum(sig$logFC < 0, na.rm = TRUE)

    de_summary <- rbind(
      de_summary,
      data.frame(
        contrast = contrast_name,
        day_group = day_i,
        treatment = dose_i,
        total_sig = nrow(sig),
        upregulated = n_up,
        downregulated = n_down,
        stringsAsFactors = FALSE
      )
    )

    # Volcano
    res$status <- "Not significant"
    res$status[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$logFC >= 0.58] <- "Up"
    res$status[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$logFC <= -0.58] <- "Down"

    res$neglog10padj <- -log10(pmax(res$adj.P.Val, .Machine$double.xmin))

    p_vol <- ggplot(res, aes(x = logFC, y = neglog10padj, shape = status)) +
      geom_point(alpha = 0.6, size = 1.1) +
      geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_bw(base_size = 11) +
      xlab("log2 fold change") +
      ylab("-log10 adjusted p-value") +
      ggtitle(contrast_name)

    ggsave(file.path(out_fig, paste0(contrast_name, "_volcano.pdf")), p_vol, width = 6.5, height = 5)
    ggsave(file.path(out_fig, paste0(contrast_name, "_volcano.png")), p_vol, width = 6.5, height = 5, dpi = 300)

    message(contrast_name, ": ", nrow(sig), " significant genes | up=", n_up, " down=", n_down)
  }
}

write.csv(
  de_summary,
  file.path(out_tab, "GSE260888_limma_DEG_summary.csv"),
  row.names = FALSE
)

# DEG summary plot
if (nrow(de_summary) > 0) {
  de_summary$contrast <- factor(de_summary$contrast, levels = de_summary$contrast)

  p_deg <- ggplot(de_summary, aes(x = contrast, y = total_sig)) +
    geom_col() +
    geom_text(aes(label = total_sig), vjust = -0.4, size = 3) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    xlab("") +
    ylab("Number of significant DEGs") +
    ggtitle("GSE260888 citalopram DEG counts")

  ggsave(file.path(out_fig, "GSE260888_DEG_count_summary.pdf"), p_deg, width = 10, height = 5)
  ggsave(file.path(out_fig, "GSE260888_DEG_count_summary.png"), p_deg, width = 10, height = 5, dpi = 300)
}

# Neurodevelopmental module scoring
message("Calculating neurodevelopmental module scores...")

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

has_symbol <- !is.na(gene_annotation$gene_symbol) & gene_annotation$gene_symbol != ""
group_symbol <- toupper(gene_annotation$gene_symbol[has_symbol])

expr_symbol <- expr[has_symbol, , drop = FALSE]
expr_sum <- rowsum(expr_symbol, group = group_symbol, reorder = FALSE)
gene_counts <- as.numeric(table(group_symbol)[rownames(expr_sum)])
expr_symbol_mean <- sweep(expr_sum, 1, gene_counts, "/")

zexpr <- t(scale(t(expr_symbol_mean)))
zexpr[is.na(zexpr)] <- 0

module_scores <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
module_gene_matching <- data.frame()

for (module_name in names(modules)) {
  genes <- toupper(modules[[module_name]])
  matched <- intersect(genes, rownames(zexpr))

  module_gene_matching <- rbind(
    module_gene_matching,
    data.frame(
      module = module_name,
      input_gene = genes,
      matched = genes %in% rownames(zexpr),
      stringsAsFactors = FALSE
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
  module_gene_matching,
  file.path(out_tab, "GSE260888_module_gene_matching.csv"),
  row.names = FALSE
)

module_scores <- merge(sample_meta, module_scores, by = "sample_id")

write.csv(
  module_scores,
  file.path(out_tab, "GSE260888_sample_level_module_scores.csv"),
  row.names = FALSE
)

score_cols <- names(modules)

condition_scores <- aggregate(
  module_scores[, score_cols, drop = FALSE],
  by = list(day_group = module_scores$day_group, treatment = module_scores$treatment),
  FUN = mean,
  na.rm = TRUE
)

write.csv(
  condition_scores,
  file.path(out_tab, "GSE260888_condition_mean_module_scores.csv"),
  row.names = FALSE
)

# Module stats per day/dose vs control
module_stats <- data.frame()

for (day_i in days_to_test) {
  for (dose_i in doses) {
    sub <- module_scores[module_scores$day_group == day_i & module_scores$treatment %in% c("CTRL", dose_i), ]

    if (nrow(sub) == 0 || length(unique(sub$treatment)) < 2) next

    sub$treatment <- factor(as.character(sub$treatment), levels = c("CTRL", dose_i))

    for (m in score_cols) {
      df <- sub[, c("sample_id", "treatment", m), drop = FALSE]
      colnames(df)[colnames(df) == m] <- "score"
      df <- df[!is.na(df$score), ]

      if (nrow(df) < 4) next

      fit <- lm(score ~ treatment, data = df)
      ct <- summary(fit)$coefficients

      term <- paste0("treatment", dose_i)

      mean_ctrl <- mean(df$score[df$treatment == "CTRL"], na.rm = TRUE)
      mean_treat <- mean(df$score[df$treatment == dose_i], na.rm = TRUE)

      if (term %in% rownames(ct)) {
        module_stats <- rbind(
          module_stats,
          data.frame(
            day_group = day_i,
            treatment = dose_i,
            contrast = paste0(day_i, "_", dose_i, "_vs_CTRL"),
            module = m,
            mean_control = mean_ctrl,
            mean_treated = mean_treat,
            mean_difference = mean_treat - mean_ctrl,
            estimate = ct[term, "Estimate"],
            std_error = ct[term, "Std. Error"],
            t_value = ct[term, "t value"],
            p_value = ct[term, "Pr(>|t|)"],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

module_stats$padj_BH <- p.adjust(module_stats$p_value, method = "BH")

write.csv(
  module_stats,
  file.path(out_tab, "GSE260888_module_statistics_vs_day_matched_control.csv"),
  row.names = FALSE
)

# Module heatmaps
long_condition <- data.frame()

for (m in score_cols) {
  long_condition <- rbind(
    long_condition,
    data.frame(
      day_group = condition_scores$day_group,
      treatment = condition_scores$treatment,
      module = m,
      score = condition_scores[[m]],
      stringsAsFactors = FALSE
    )
  )
}

long_condition$day_group <- factor(long_condition$day_group, levels = c("D0", "D6", "D10", "D13"))
long_condition$treatment <- factor(long_condition$treatment, levels = c("CTRL", "CIT50", "CIT100", "CIT200", "CIT400"))

p_mod <- ggplot(long_condition, aes(x = treatment, y = module, fill = score)) +
  geom_tile() +
  geom_text(aes(label = round(score, 2)), size = 2.5) +
  facet_wrap(~ day_group, nrow = 1) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("GSE260888 citalopram neurodevelopmental module scores")

ggsave(file.path(out_fig, "GSE260888_condition_mean_module_heatmap.pdf"), p_mod, width = 14, height = 6)
ggsave(file.path(out_fig, "GSE260888_condition_mean_module_heatmap.png"), p_mod, width = 14, height = 6, dpi = 300)

# Module difference plot for D10 and D13
plot_stats <- module_stats[module_stats$day_group %in% c("D10", "D13"), ]

if (nrow(plot_stats) > 0) {
  p_diff <- ggplot(plot_stats, aes(x = mean_difference, y = module)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.3) +
    facet_grid(day_group ~ treatment) +
    theme_bw(base_size = 9) +
    xlab("Mean module-score difference vs day-matched control") +
    ylab("") +
    ggtitle("GSE260888 citalopram module shifts vs day-matched controls")

  ggsave(file.path(out_fig, "GSE260888_module_difference_vs_control_D10_D13.pdf"), p_diff, width = 14, height = 8)
  ggsave(file.path(out_fig, "GSE260888_module_difference_vs_control_D10_D13.png"), p_diff, width = 14, height = 8, dpi = 300)
}

capture.output(sessionInfo(), file = file.path(out_log, "GSE260888_limma_module_sessionInfo.txt"))

message("GSE260888 limma and module analysis completed.")
