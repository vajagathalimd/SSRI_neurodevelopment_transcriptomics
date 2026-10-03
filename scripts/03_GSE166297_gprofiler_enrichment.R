#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gprofiler2)
  library(ggplot2)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

de_dir  <- "09_differential_expression/GSE166297"
tab_dir <- "11_pathway_enrichment/GSE166297_gprofiler"
fig_dir <- "14_figures/GSE166297_enrichment"

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

get_gene_vector <- function(df) {
  if ("gene_name" %in% colnames(df) && any(!is.na(df$gene_name) & df$gene_name != "")) {
    genes <- ifelse(is.na(df$gene_name) | df$gene_name == "", df$gene_id, df$gene_name)
  } else {
    genes <- df$gene_id
  }

  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- sub("\\..*$", "", genes)
  return(genes)
}

run_gprofiler <- function(genes, name) {
  message("Running enrichment for ", name, " | genes: ", length(genes))

  if (length(genes) < 5) {
    message("Skipped: too few genes.")
    return(NULL)
  }

  res <- tryCatch(
    gost(
      query = genes,
      organism = "hsapiens",
      correction_method = "g_SCS",
      significant = TRUE,
      sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP")
    ),
    error = function(e) {
      message("g:Profiler error for ", name, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    message("No significant enrichment found for ", name)
    return(NULL)
  }

  enrich <- res$result
  enrich <- enrich[order(enrich$p_value), ]

  write.csv(
    enrich,
    file.path(tab_dir, paste0(name, "_gprofiler_enrichment.csv")),
    row.names = FALSE
  )

  top <- head(enrich, 20)

  p <- ggplot(top, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    xlab("") +
    ylab("-log10 adjusted p-value") +
    ggtitle(paste0("Top enriched pathways: ", name))

  ggsave(
    file.path(fig_dir, paste0(name, "_top20_enrichment.pdf")),
    p,
    width = 8,
    height = 6
  )

  ggsave(
    file.path(fig_dir, paste0(name, "_top20_enrichment.png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )

  message("Saved enrichment for ", name, " | significant terms: ", nrow(enrich))
  return(enrich)
}

ic5_sig <- read.csv(
  file.path(de_dir, "FLX_IC5_vs_Control_D10_significant_FDR05_log2FC058.csv"),
  check.names = FALSE
)

ic20_sig <- read.csv(
  file.path(de_dir, "FLX_IC20_100_vs_Control_D10_significant_FDR05_log2FC058.csv"),
  check.names = FALSE
)

ic5_genes <- get_gene_vector(ic5_sig)
ic20_genes <- get_gene_vector(ic20_sig)
shared_genes <- intersect(ic5_genes, ic20_genes)

ic5_up <- get_gene_vector(ic5_sig[ic5_sig$log2FoldChange > 0, ])
ic5_down <- get_gene_vector(ic5_sig[ic5_sig$log2FoldChange < 0, ])

ic20_up <- get_gene_vector(ic20_sig[ic20_sig$log2FoldChange > 0, ])
ic20_down <- get_gene_vector(ic20_sig[ic20_sig$log2FoldChange < 0, ])

input_summary <- data.frame(
  gene_set = c(
    "FLX_IC5_all",
    "FLX_IC20_100_all",
    "Shared",
    "FLX_IC5_up",
    "FLX_IC5_down",
    "FLX_IC20_100_up",
    "FLX_IC20_100_down"
  ),
  n_genes = c(
    length(ic5_genes),
    length(ic20_genes),
    length(shared_genes),
    length(ic5_up),
    length(ic5_down),
    length(ic20_up),
    length(ic20_down)
  )
)

write.csv(
  input_summary,
  file.path(tab_dir, "GSE166297_gprofiler_input_gene_counts.csv"),
  row.names = FALSE
)

run_gprofiler(ic5_genes, "FLX_IC5_all_significant_DEGs")
run_gprofiler(ic20_genes, "FLX_IC20_100_all_significant_DEGs")
run_gprofiler(shared_genes, "Shared_FLX_IC5_and_IC20_100_DEGs")

run_gprofiler(ic5_up, "FLX_IC5_upregulated_DEGs")
run_gprofiler(ic5_down, "FLX_IC5_downregulated_DEGs")

run_gprofiler(ic20_up, "FLX_IC20_100_upregulated_DEGs")
run_gprofiler(ic20_down, "FLX_IC20_100_downregulated_DEGs")

message("g:Profiler enrichment analysis completed.")
