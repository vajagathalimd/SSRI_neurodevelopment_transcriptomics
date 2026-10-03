#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

enrich_dir <- "11_pathway_enrichment/GSE166297_gprofiler"
fig_dir <- "14_figures/GSE166297_enrichment"
tab_dir <- "15_tables/GSE166297_enrichment"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(enrich_dir, pattern = "_gprofiler_enrichment.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No g:Profiler enrichment CSV files found.")
}

all_top <- data.frame()

for (file in files) {
  df <- read.csv(file, check.names = FALSE)

  if (nrow(df) == 0) next

  df$p_value <- as.numeric(df$p_value)
  df <- df[!is.na(df$p_value), ]
  df <- df[order(df$p_value), ]

  name <- basename(file)
  name <- sub("_gprofiler_enrichment.csv", "", name)

  df$gene_set <- name

  top20 <- head(df, 20)
  all_top <- rbind(all_top, top20)

  write.csv(
    top20,
    file.path(tab_dir, paste0(name, "_top20_enrichment_terms.csv")),
    row.names = FALSE
  )

  top20$term_name_wrapped <- sapply(top20$term_name, function(x) {
    paste(strwrap(x, width = 55), collapse = "\n")
  })

  top20$term_name_wrapped <- factor(
    top20$term_name_wrapped,
    levels = rev(top20$term_name_wrapped)
  )

  p <- ggplot(top20, aes(x = term_name_wrapped, y = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    xlab("") +
    ylab("-log10 adjusted p-value") +
    ggtitle(paste0("Top enriched terms: ", name))

  ggsave(
    file.path(fig_dir, paste0(name, "_top20_enrichment.pdf")),
    p,
    width = 9,
    height = 6
  )

  ggsave(
    file.path(fig_dir, paste0(name, "_top20_enrichment.png")),
    p,
    width = 9,
    height = 6,
    dpi = 300
  )

  message("Saved enrichment plot: ", name)
}

write.csv(
  all_top,
  file.path(tab_dir, "GSE166297_all_gene_sets_top20_enrichment_terms.csv"),
  row.names = FALSE
)

message("All enrichment plots completed.")
