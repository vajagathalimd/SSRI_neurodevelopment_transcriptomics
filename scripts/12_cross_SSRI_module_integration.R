#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

gse166_file <- "15_tables/GSE166297_final_summary/GSE166297_neurodevelopmental_module_statistics.csv"
gse260_file <- "15_tables/GSE260888/GSE260888_module_statistics_vs_day_matched_control.csv"

gse166_deg_file <- "15_tables/GSE166297_final_summary/GSE166297_DEG_direction_summary.csv"
gse260_deg_file <- "15_tables/GSE260888/GSE260888_limma_DEG_summary.csv"

out_tab <- "15_tables/cross_SSRI_integration"
out_fig <- "14_figures/cross_SSRI_integration"

dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

message("Reading GSE166297 fluoxetine module statistics...")
flx <- read.csv(gse166_file, check.names = FALSE)

flx$dataset <- "GSE166297"
flx$drug <- "Fluoxetine"
flx$timepoint <- "D10"
flx$treatment <- ifelse(grepl("FLX_IC5", flx$contrast), "FLX_IC5", "FLX_IC20_100")
flx$comparison_label <- paste0("FLX_", flx$treatment)
flx$mean_difference_vs_control <- flx$mean_difference
flx$adjusted_p <- flx$padj_BH

flx_keep <- flx[, c(
  "dataset", "drug", "timepoint", "treatment", "comparison_label",
  "module", "mean_difference_vs_control", "p_value", "adjusted_p"
)]

message("Reading GSE260888 citalopram module statistics...")
cit <- read.csv(gse260_file, check.names = FALSE)

cit$dataset <- "GSE260888"
cit$drug <- "Citalopram"
cit$timepoint <- cit$day_group
cit$comparison_label <- paste0(cit$day_group, "_", cit$treatment)
cit$mean_difference_vs_control <- cit$mean_difference
cit$adjusted_p <- cit$padj_BH

cit_keep <- cit[, c(
  "dataset", "drug", "timepoint", "treatment", "comparison_label",
  "module", "mean_difference_vs_control", "p_value", "adjusted_p"
)]

combined <- rbind(flx_keep, cit_keep)

combined$significant <- !is.na(combined$adjusted_p) & combined$adjusted_p < 0.05
combined$direction <- ifelse(
  combined$mean_difference_vs_control > 0, "Increased",
  ifelse(combined$mean_difference_vs_control < 0, "Decreased", "No change")
)

combined$sig_label <- ifelse(combined$significant, "*", "")

module_order <- c(
  "Neural_progenitor_cell_cycle",
  "Neuronal_differentiation",
  "Axon_neurite_outgrowth",
  "Synaptic_function",
  "Glutamatergic_signalling",
  "GABAergic_signalling",
  "Serotonergic_signalling",
  "Astrocyte_glial_identity",
  "Oligodendrocyte_myelin",
  "Oxidative_stress_apoptosis",
  "Inflammatory_chemokine_signalling"
)

combined$module <- factor(combined$module, levels = rev(module_order))

write.csv(
  combined,
  file.path(out_tab, "SSRI_combined_module_shift_table.csv"),
  row.names = FALSE
)

# Focused D10/D13 table for manuscript
focused <- combined[
  combined$dataset == "GSE166297" |
    (combined$dataset == "GSE260888" & combined$timepoint %in% c("D10", "D13")),
]

focused_order <- c(
  "FLX_FLX_IC5",
  "FLX_FLX_IC20_100",
  "D10_CIT50",
  "D10_CIT100",
  "D10_CIT200",
  "D10_CIT400",
  "D13_CIT50",
  "D13_CIT100",
  "D13_CIT200",
  "D13_CIT400"
)

focused$comparison_label <- factor(focused$comparison_label, levels = focused_order)

write.csv(
  focused,
  file.path(out_tab, "SSRI_focused_D10_D13_module_shift_table.csv"),
  row.names = FALSE
)

# Heatmap: focused D10/D13 comparison
p_heat <- ggplot(focused, aes(x = comparison_label, y = module, fill = mean_difference_vs_control)) +
  geom_tile() +
  geom_text(aes(label = sig_label), size = 5) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white")
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cross-SSRI neurodevelopmental module shifts")

ggsave(
  file.path(out_fig, "SSRI_focused_D10_D13_module_shift_heatmap.pdf"),
  p_heat,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_fig, "SSRI_focused_D10_D13_module_shift_heatmap.png"),
  p_heat,
  width = 11,
  height = 6,
  dpi = 300
)

# All citalopram time points + fluoxetine
all_order <- c(
  "FLX_FLX_IC5",
  "FLX_FLX_IC20_100",
  "D6_CIT50",
  "D6_CIT100",
  "D6_CIT200",
  "D6_CIT400",
  "D10_CIT50",
  "D10_CIT100",
  "D10_CIT200",
  "D10_CIT400",
  "D13_CIT50",
  "D13_CIT100",
  "D13_CIT200",
  "D13_CIT400"
)

combined$comparison_label <- factor(combined$comparison_label, levels = all_order)

p_all <- ggplot(combined, aes(x = comparison_label, y = module, fill = mean_difference_vs_control)) +
  geom_tile() +
  geom_text(aes(label = sig_label), size = 4) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("Cross-SSRI module shifts across all analysed conditions")

ggsave(
  file.path(out_fig, "SSRI_all_conditions_module_shift_heatmap.pdf"),
  p_all,
  width = 13,
  height = 6
)

ggsave(
  file.path(out_fig, "SSRI_all_conditions_module_shift_heatmap.png"),
  p_all,
  width = 13,
  height = 6,
  dpi = 300
)

# Module direction summary
direction_summary <- aggregate(
  significant ~ dataset + drug + timepoint + treatment + direction,
  data = combined,
  FUN = sum
)

colnames(direction_summary)[colnames(direction_summary) == "significant"] <- "n_significant_modules"

write.csv(
  direction_summary,
  file.path(out_tab, "SSRI_significant_module_direction_summary.csv"),
  row.names = FALSE
)

# Correlation between fluoxetine and each citalopram condition
message("Calculating module-shift correlations between fluoxetine and citalopram conditions...")

wide <- reshape(
  combined[, c("module", "comparison_label", "mean_difference_vs_control")],
  idvar = "module",
  timevar = "comparison_label",
  direction = "wide"
)

colnames(wide) <- sub("^mean_difference_vs_control\\.", "", colnames(wide))

flx_cols <- c("FLX_FLX_IC5", "FLX_FLX_IC20_100")
cit_cols <- grep("^D[0-9]+_CIT", colnames(wide), value = TRUE)

cor_summary <- data.frame()

for (fcol in flx_cols) {
  for (ccol in cit_cols) {
    x <- wide[[fcol]]
    y <- wide[[ccol]]
    ok <- !is.na(x) & !is.na(y)

    if (sum(ok) >= 4) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))

      cor_summary <- rbind(
        cor_summary,
        data.frame(
          fluoxetine_condition = fcol,
          citalopram_condition = ccol,
          n_modules = sum(ok),
          correlation_r = unname(ct$estimate),
          p_value = ct$p.value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

cor_summary$padj_BH <- p.adjust(cor_summary$p_value, method = "BH")

write.csv(
  cor_summary,
  file.path(out_tab, "SSRI_module_shift_correlation_summary.csv"),
  row.names = FALSE
)

# Correlation heatmap
if (nrow(cor_summary) > 0) {
  cor_summary$fluoxetine_condition <- factor(cor_summary$fluoxetine_condition, levels = flx_cols)
  cor_summary$citalopram_condition <- factor(cor_summary$citalopram_condition, levels = cit_cols)

  p_cor <- ggplot(cor_summary, aes(x = citalopram_condition, y = fluoxetine_condition, fill = correlation_r)) +
    geom_tile() +
    geom_text(aes(label = round(correlation_r, 2)), size = 3) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Citalopram condition") +
    ylab("Fluoxetine condition") +
    ggtitle("Correlation of module-shift profiles between fluoxetine and citalopram")

  ggsave(
    file.path(out_fig, "SSRI_module_shift_correlation_heatmap.pdf"),
    p_cor,
    width = 10,
    height = 4
  )

  ggsave(
    file.path(out_fig, "SSRI_module_shift_correlation_heatmap.png"),
    p_cor,
    width = 10,
    height = 4,
    dpi = 300
  )
}

# Combined DEG count summary
message("Combining DEG summaries...")

flx_deg <- read.csv(gse166_deg_file, check.names = FALSE)
flx_deg$dataset <- "GSE166297"
flx_deg$drug <- "Fluoxetine"
flx_deg$timepoint <- "D10"
flx_deg$treatment <- ifelse(grepl("FLX_IC5", flx_deg$contrast), "FLX_IC5", "FLX_IC20_100")
flx_deg$total_sig <- flx_deg$total

flx_deg_keep <- flx_deg[, c("dataset", "drug", "timepoint", "treatment", "contrast", "total_sig", "upregulated", "downregulated")]

cit_deg <- read.csv(gse260_deg_file, check.names = FALSE)
cit_deg$dataset <- "GSE260888"
cit_deg$drug <- "Citalopram"
cit_deg$timepoint <- cit_deg$day_group

cit_deg_keep <- cit_deg[, c("dataset", "drug", "timepoint", "treatment", "contrast", "total_sig", "upregulated", "downregulated")]

deg_combined <- rbind(flx_deg_keep, cit_deg_keep)

write.csv(
  deg_combined,
  file.path(out_tab, "SSRI_combined_DEG_count_summary.csv"),
  row.names = FALSE
)

deg_combined$contrast <- factor(deg_combined$contrast, levels = deg_combined$contrast)

p_deg <- ggplot(deg_combined, aes(x = contrast, y = total_sig)) +
  geom_col() +
  geom_text(aes(label = total_sig), vjust = -0.3, size = 3) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  xlab("") +
  ylab("Number of significant DEGs") +
  ggtitle("DEG burden across fluoxetine and citalopram conditions")

ggsave(
  file.path(out_fig, "SSRI_combined_DEG_count_summary.pdf"),
  p_deg,
  width = 12,
  height = 5
)

ggsave(
  file.path(out_fig, "SSRI_combined_DEG_count_summary.png"),
  p_deg,
  width = 12,
  height = 5,
  dpi = 300
)

message("Cross-SSRI integration completed.")
