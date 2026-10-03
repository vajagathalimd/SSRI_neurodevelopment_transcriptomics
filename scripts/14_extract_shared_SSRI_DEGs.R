#!/usr/bin/env Rscript

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

in_dir <- "15_tables/cross_SSRI_gene_level"
out_dir <- "15_tables/cross_SSRI_gene_level"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- c(
  "FLX_IC5_vs_Control_D10__D13_CIT400_vs_CTRL_merged_gene_signature.csv",
  "FLX_IC20_100_vs_Control_D10__D13_CIT400_vs_CTRL_merged_gene_signature.csv",
  "FLX_IC5_vs_Control_D10__D10_CIT200_vs_CTRL_merged_gene_signature.csv",
  "FLX_IC20_100_vs_Control_D10__D10_CIT200_vs_CTRL_merged_gene_signature.csv",
  "FLX_IC20_100_vs_Control_D10__D13_CIT100_vs_CTRL_merged_gene_signature.csv"
)

all_shared <- data.frame()

for (f in pairs) {
  path <- file.path(in_dir, f)
  if (!file.exists(path)) next

  df <- read.csv(path, check.names = FALSE)

  shared <- df[
    df$significant_FLX == TRUE &
    df$significant_CIT == TRUE,
  ]

  if (nrow(shared) == 0) next

  shared$pair <- sub("_merged_gene_signature.csv$", "", f)
  shared$same_direction <- shared$direction_FLX == shared$direction_CIT

  shared <- shared[order(shared$same_direction, -abs(shared$logFC_FLX), -abs(shared$logFC_CIT)), ]

  all_shared <- rbind(all_shared, shared)
}

write.csv(
  all_shared,
  file.path(out_dir, "SSRI_top_pair_shared_significant_DEGs.csv"),
  row.names = FALSE
)

same_dir <- all_shared[all_shared$same_direction == TRUE, ]

write.csv(
  same_dir,
  file.path(out_dir, "SSRI_top_pair_shared_same_direction_DEGs.csv"),
  row.names = FALSE
)

down_down <- same_dir[
  same_dir$direction_FLX == "Down" &
  same_dir$direction_CIT == "Down",
]

write.csv(
  down_down,
  file.path(out_dir, "SSRI_top_pair_shared_downregulated_DEGs.csv"),
  row.names = FALSE
)

cat("Total shared significant rows:", nrow(all_shared), "\n")
cat("Same-direction rows:", nrow(same_dir), "\n")
cat("Down-down rows:", nrow(down_down), "\n")
cat("Top shared downregulated genes:\n")
print(unique(head(down_down$gene_key, 50)))
