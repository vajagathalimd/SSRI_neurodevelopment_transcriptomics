#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

module_file <- "15_tables/GSE166297_module_scores/GSE166297_sample_level_neurodevelopmental_module_scores.csv"
overlap_file <- "15_tables/GSE166297_DESeq2/GSE166297_DEG_overlap_summary.csv"
direction_file <- "15_tables/GSE166297_DESeq2/GSE166297_DEG_direction_summary.csv"
top_deg_file <- "15_tables/GSE166297_DESeq2/GSE166297_top_25_up_down_DEGs_each_contrast.csv"

enrich_dir <- "11_pathway_enrichment/GSE166297_gprofiler"

out_tab <- "15_tables/GSE166297_final_summary"
out_fig <- "14_figures/GSE166297_module_stats"

dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

module_scores <- read.csv(module_file, check.names = FALSE)

module_cols <- setdiff(
  colnames(module_scores),
  c("sample_id", "condition_group", "replicate_number", "batch")
)

module_scores$condition_group <- factor(
  module_scores$condition_group,
  levels = c("Control_D10", "FLX_IC5", "FLX_IC20_100")
)

if ("batch" %in% colnames(module_scores)) {
  module_scores$batch <- as.factor(module_scores$batch)
}

# Module statistics: control as reference
stats_all <- data.frame()

for (m in module_cols) {
  df <- module_scores[, c("sample_id", "condition_group", "batch", m), drop = FALSE]
  colnames(df)[colnames(df) == m] <- "score"
  df <- df[!is.na(df$score), ]

  use_batch <- FALSE
  if ("batch" %in% colnames(df) && length(unique(df$batch)) > 1) {
    mm <- model.matrix(~ batch + condition_group, data = df)
    if (qr(mm)$rank == ncol(mm)) {
      use_batch <- TRUE
    }
  }

  if (use_batch) {
    fit <- lm(score ~ batch + condition_group, data = df)
    design_used <- "score ~ batch + condition_group"
  } else {
    fit <- lm(score ~ condition_group, data = df)
    design_used <- "score ~ condition_group"
  }

  coef_tab <- summary(fit)$coefficients

  mean_control <- mean(df$score[df$condition_group == "Control_D10"], na.rm = TRUE)
  mean_ic5 <- mean(df$score[df$condition_group == "FLX_IC5"], na.rm = TRUE)
  mean_ic20 <- mean(df$score[df$condition_group == "FLX_IC20_100"], na.rm = TRUE)

  get_row <- function(term, contrast_name, mean_treated, mean_control) {
    if (term %in% rownames(coef_tab)) {
      data.frame(
        module = m,
        contrast = contrast_name,
        mean_control = mean_control,
        mean_treated = mean_treated,
        mean_difference = mean_treated - mean_control,
        estimate_from_model = coef_tab[term, "Estimate"],
        std_error = coef_tab[term, "Std. Error"],
        t_value = coef_tab[term, "t value"],
        p_value = coef_tab[term, "Pr(>|t|)"],
        design_used = design_used,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        module = m,
        contrast = contrast_name,
        mean_control = mean_control,
        mean_treated = mean_treated,
        mean_difference = mean_treated - mean_control,
        estimate_from_model = NA,
        std_error = NA,
        t_value = NA,
        p_value = NA,
        design_used = design_used,
        stringsAsFactors = FALSE
      )
    }
  }

  stats_all <- rbind(
    stats_all,
    get_row("condition_groupFLX_IC5", "FLX_IC5_vs_Control_D10", mean_ic5, mean_control),
    get_row("condition_groupFLX_IC20_100", "FLX_IC20_100_vs_Control_D10", mean_ic20, mean_control)
  )
}

stats_all$padj_BH <- p.adjust(stats_all$p_value, method = "BH")

write.csv(
  stats_all,
  file.path(out_tab, "GSE166297_neurodevelopmental_module_statistics.csv"),
  row.names = FALSE
)

# Plot module mean differences
plot_df <- stats_all
plot_df$contrast <- factor(
  plot_df$contrast,
  levels = c("FLX_IC5_vs_Control_D10", "FLX_IC20_100_vs_Control_D10")
)

p <- ggplot(plot_df, aes(x = mean_difference, y = module)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  facet_wrap(~ contrast, ncol = 1) +
  theme_bw(base_size = 11) +
  xlab("Mean module-score difference vs Control_D10") +
  ylab("") +
  ggtitle("Fluoxetine-associated neurodevelopmental module shifts")

ggsave(
  file.path(out_fig, "GSE166297_module_mean_difference_vs_control.pdf"),
  p,
  width = 8,
  height = 8
)

ggsave(
  file.path(out_fig, "GSE166297_module_mean_difference_vs_control.png"),
  p,
  width = 8,
  height = 8,
  dpi = 300
)

# Enrichment top-term clean summary
enrich_files <- list.files(enrich_dir, pattern = "_gprofiler_enrichment.csv$", full.names = TRUE)

top_enrichment_all <- data.frame()

for (f in enrich_files) {
  df <- read.csv(f, check.names = FALSE)

  if (nrow(df) == 0) next

  df$p_value <- as.numeric(df$p_value)
  df <- df[order(df$p_value), ]

  gene_set <- basename(f)
  gene_set <- sub("_gprofiler_enrichment.csv", "", gene_set)

  top <- head(df, 10)
  top$gene_set <- gene_set

  keep <- intersect(
    c("gene_set", "source", "native", "term_name", "p_value", "term_size", "query_size", "intersection_size", "precision", "recall"),
    colnames(top)
  )

  top_enrichment_all <- rbind(top_enrichment_all, top[, keep, drop = FALSE])
}

write.csv(
  top_enrichment_all,
  file.path(out_tab, "GSE166297_top10_enrichment_terms_all_gene_sets.csv"),
  row.names = FALSE
)

# Copy DEG summaries into final summary folder
if (file.exists(overlap_file)) {
  file.copy(overlap_file, file.path(out_tab, basename(overlap_file)), overwrite = TRUE)
}

if (file.exists(direction_file)) {
  file.copy(direction_file, file.path(out_tab, basename(direction_file)), overwrite = TRUE)
}

if (file.exists(top_deg_file)) {
  file.copy(top_deg_file, file.path(out_tab, basename(top_deg_file)), overwrite = TRUE)
}

# Create compact text summary
overlap <- read.csv(overlap_file)
direction <- read.csv(direction_file)

summary_txt <- file.path(out_tab, "GSE166297_results_summary_for_manuscript.txt")

sink(summary_txt)

cat("GSE166297 fluoxetine transcriptomic analysis summary\n")
cat("====================================================\n\n")

cat("Differential expression analysis:\n")
print(direction)
cat("\n")

cat("DEG overlap summary:\n")
print(overlap)
cat("\n")

cat("Neurodevelopmental module statistics:\n")
print(stats_all[order(stats_all$contrast, stats_all$p_value), ])
cat("\n")

cat("Top enrichment terms per gene set:\n")
print(top_enrichment_all)
cat("\n")

sink()

message("Final summary tables completed.")
message("Outputs written to: ", out_tab)
message("Figures written to: ", out_fig)
