#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

flx_dir <- "09_differential_expression/GSE166297"
cit_dir <- "09_differential_expression/GSE260888"

out_tab <- "15_tables/cross_SSRI_gene_level"
out_fig <- "14_figures/cross_SSRI_gene_level"

dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

clean_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*$", "", x)
  x[is.na(x)] <- ""
  x
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x[toupper(x) %in% c("NA", "NAN", "NULL")] <- ""
  toupper(x)
}

map_to_symbol <- function(ids, supplied_symbol = NULL) {
  ids_clean <- clean_id(ids)

  if (is.null(supplied_symbol)) {
    symbol <- rep("", length(ids_clean))
  } else {
    symbol <- clean_symbol(supplied_symbol)
  }

  missing <- symbol == ""

  # Ensembl mapping, only using keys valid in org.Hs.eg.db
  ensg <- ids_clean
  ensg_keys <- unique(ensg[missing & grepl("^ENSG", ensg)])
  valid_ensg <- intersect(ensg_keys, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL"))

  if (length(valid_ensg) > 0) {
    mapped_ensg <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = valid_ensg,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

    hit <- missing & ensg %in% names(mapped_ensg)
    symbol[hit] <- clean_symbol(mapped_ensg[ensg[hit]])
  }

  # Entrez mapping, only using keys valid in org.Hs.eg.db
  missing <- symbol == ""
  entrez <- ids_clean
  entrez_keys <- unique(entrez[missing & grepl("^[0-9]+$", entrez)])
  valid_entrez <- intersect(entrez_keys, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"))

  if (length(valid_entrez) > 0) {
    mapped_entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = valid_entrez,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    hit <- missing & entrez %in% names(mapped_entrez)
    symbol[hit] <- clean_symbol(mapped_entrez[entrez[hit]])
  }

  # If remaining IDs look like symbols, use them directly
  missing <- symbol == ""
  possible_symbol <- ids_clean
  use_direct <- missing & !grepl("^ENSG", possible_symbol) & !grepl("^[0-9]+$", possible_symbol)
  symbol[use_direct] <- clean_symbol(possible_symbol[use_direct])

  symbol
}

collapse_by_gene_key <- function(df) {
  df <- df[!is.na(df$gene_key) & df$gene_key != "", ]
  df <- df[order(df$gene_key, -abs(df$logFC)), ]
  df <- df[!duplicated(df$gene_key), ]
  df
}

read_flx <- function(file, contrast_name) {
  df <- read.csv(file, check.names = FALSE)

  supplied_symbol <- NULL
  if ("gene_name" %in% colnames(df)) supplied_symbol <- df$gene_name

  gene_key <- map_to_symbol(df$gene_id, supplied_symbol)

  out <- data.frame(
    contrast = contrast_name,
    original_gene_id = df$gene_id,
    gene_key = gene_key,
    logFC = df$log2FoldChange,
    pvalue = df$pvalue,
    padj = df$padj,
    significant = !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= 0.58,
    direction = ifelse(df$log2FoldChange > 0, "Up", ifelse(df$log2FoldChange < 0, "Down", "NoChange")),
    stringsAsFactors = FALSE
  )

  collapse_by_gene_key(out)
}

read_cit <- function(file, contrast_name) {
  df <- read.csv(file, check.names = FALSE)

  supplied_symbol <- NULL
  if ("gene_symbol" %in% colnames(df)) supplied_symbol <- df$gene_symbol

  id_col <- if ("ensembl_id" %in% colnames(df)) df$ensembl_id else df$gene_id

  gene_key <- map_to_symbol(id_col, supplied_symbol)

  out <- data.frame(
    contrast = contrast_name,
    original_gene_id = id_col,
    gene_key = gene_key,
    logFC = df$logFC,
    pvalue = df$P.Value,
    padj = df$adj.P.Val,
    significant = !is.na(df$adj.P.Val) & df$adj.P.Val < 0.05 & abs(df$logFC) >= 0.58,
    direction = ifelse(df$logFC > 0, "Up", ifelse(df$logFC < 0, "Down", "NoChange")),
    stringsAsFactors = FALSE
  )

  collapse_by_gene_key(out)
}

flx_files <- c(
  FLX_IC5_vs_Control_D10 = file.path(flx_dir, "FLX_IC5_vs_Control_D10_all_genes.csv"),
  FLX_IC20_100_vs_Control_D10 = file.path(flx_dir, "FLX_IC20_100_vs_Control_D10_all_genes.csv")
)

cit_files <- list.files(
  cit_dir,
  pattern = "_all_genes.csv$",
  full.names = TRUE
)

cit_names <- basename(cit_files)
cit_names <- sub("_all_genes.csv$", "", cit_names)

flx_list <- list()
for (nm in names(flx_files)) {
  flx_list[[nm]] <- read_flx(flx_files[[nm]], nm)
}

cit_list <- list()
for (i in seq_along(cit_files)) {
  cit_list[[cit_names[i]]] <- read_cit(cit_files[i], cit_names[i])
}

diagnostics <- data.frame()

for (nm in names(flx_list)) {
  diagnostics <- rbind(
    diagnostics,
    data.frame(
      dataset = "GSE166297",
      contrast = nm,
      n_gene_keys = nrow(flx_list[[nm]]),
      n_significant = sum(flx_list[[nm]]$significant),
      example_gene_keys = paste(head(flx_list[[nm]]$gene_key, 10), collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

for (nm in names(cit_list)) {
  diagnostics <- rbind(
    diagnostics,
    data.frame(
      dataset = "GSE260888",
      contrast = nm,
      n_gene_keys = nrow(cit_list[[nm]]),
      n_significant = sum(cit_list[[nm]]$significant),
      example_gene_keys = paste(head(cit_list[[nm]]$gene_key, 10), collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

write.csv(
  diagnostics,
  file.path(out_tab, "SSRI_gene_key_mapping_diagnostics.csv"),
  row.names = FALSE
)

cor_summary <- data.frame()
overlap_summary <- data.frame()

for (f_name in names(flx_list)) {
  f_df <- flx_list[[f_name]]

  for (c_name in names(cit_list)) {
    c_df <- cit_list[[c_name]]

    merged <- merge(
      f_df[, c("gene_key", "original_gene_id", "logFC", "padj", "significant", "direction")],
      c_df[, c("gene_key", "original_gene_id", "logFC", "padj", "significant", "direction")],
      by = "gene_key",
      suffixes = c("_FLX", "_CIT")
    )

    merged <- merged[!is.na(merged$logFC_FLX) & !is.na(merged$logFC_CIT), ]

    if (nrow(merged) >= 10) {
      pearson <- suppressWarnings(cor.test(merged$logFC_FLX, merged$logFC_CIT, method = "pearson"))
      spearman <- suppressWarnings(cor.test(merged$logFC_FLX, merged$logFC_CIT, method = "spearman"))

      cor_summary <- rbind(
        cor_summary,
        data.frame(
          fluoxetine_contrast = f_name,
          citalopram_contrast = c_name,
          n_common_genes = nrow(merged),
          pearson_r = unname(pearson$estimate),
          pearson_p = pearson$p.value,
          spearman_rho = unname(spearman$estimate),
          spearman_p = spearman$p.value,
          stringsAsFactors = FALSE
        )
      )
    }

    flx_sig <- merged$significant_FLX
    cit_sig <- merged$significant_CIT

    shared_sig <- merged[flx_sig & cit_sig, ]
    same_dir <- shared_sig[shared_sig$direction_FLX == shared_sig$direction_CIT, ]
    opposite_dir <- shared_sig[shared_sig$direction_FLX != shared_sig$direction_CIT, ]

    down_down <- shared_sig[shared_sig$direction_FLX == "Down" & shared_sig$direction_CIT == "Down", ]
    up_up <- shared_sig[shared_sig$direction_FLX == "Up" & shared_sig$direction_CIT == "Up", ]

    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        fluoxetine_contrast = f_name,
        citalopram_contrast = c_name,
        common_gene_universe = nrow(merged),
        flx_sig = sum(flx_sig, na.rm = TRUE),
        cit_sig = sum(cit_sig, na.rm = TRUE),
        shared_sig = nrow(shared_sig),
        shared_same_direction = nrow(same_dir),
        shared_opposite_direction = nrow(opposite_dir),
        shared_down_down = nrow(down_down),
        shared_up_up = nrow(up_up),
        stringsAsFactors = FALSE
      )
    )

    important_cit <- c(
      "D10_CIT200_vs_CTRL",
      "D10_CIT400_vs_CTRL",
      "D13_CIT100_vs_CTRL",
      "D13_CIT400_vs_CTRL",
      "D6_CIT400_vs_CTRL"
    )

    if (c_name %in% important_cit) {
      write.csv(
        merged,
        file.path(out_tab, paste0(f_name, "__", c_name, "_merged_gene_signature.csv")),
        row.names = FALSE
      )
    }
  }
}

cor_summary$pearson_padj_BH <- p.adjust(cor_summary$pearson_p, method = "BH")
cor_summary$spearman_padj_BH <- p.adjust(cor_summary$spearman_p, method = "BH")

write.csv(
  cor_summary,
  file.path(out_tab, "SSRI_gene_level_logFC_correlation_summary.csv"),
  row.names = FALSE
)

write.csv(
  overlap_summary,
  file.path(out_tab, "SSRI_significant_DEG_overlap_summary.csv"),
  row.names = FALSE
)

if (nrow(cor_summary) > 0) {
  cor_summary$citalopram_contrast <- factor(
    cor_summary$citalopram_contrast,
    levels = unique(cor_summary$citalopram_contrast)
  )

  p1 <- ggplot(cor_summary, aes(x = citalopram_contrast, y = fluoxetine_contrast, fill = spearman_rho)) +
    geom_tile() +
    geom_text(aes(label = round(spearman_rho, 2)), size = 3) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    xlab("Citalopram contrast") +
    ylab("Fluoxetine contrast") +
    ggtitle("Gene-level SSRI signature correlation: Spearman rho")

  ggsave(file.path(out_fig, "SSRI_gene_level_spearman_correlation_heatmap.pdf"), p1, width = 10, height = 4)
  ggsave(file.path(out_fig, "SSRI_gene_level_spearman_correlation_heatmap.png"), p1, width = 10, height = 4, dpi = 300)

  p2 <- ggplot(cor_summary, aes(x = citalopram_contrast, y = fluoxetine_contrast, fill = pearson_r)) +
    geom_tile() +
    geom_text(aes(label = round(pearson_r, 2)), size = 3) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    xlab("Citalopram contrast") +
    ylab("Fluoxetine contrast") +
    ggtitle("Gene-level SSRI signature correlation: Pearson r")

  ggsave(file.path(out_fig, "SSRI_gene_level_pearson_correlation_heatmap.pdf"), p2, width = 10, height = 4)
  ggsave(file.path(out_fig, "SSRI_gene_level_pearson_correlation_heatmap.png"), p2, width = 10, height = 4, dpi = 300)
}

if (nrow(overlap_summary) > 0) {
  overlap_summary$pair <- paste(overlap_summary$fluoxetine_contrast, overlap_summary$citalopram_contrast, sep = " | ")

  p3 <- ggplot(overlap_summary, aes(x = pair, y = shared_sig)) +
    geom_col() +
    geom_text(aes(label = shared_sig), vjust = -0.3, size = 2.7) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
    xlab("") +
    ylab("Shared significant DEGs") +
    ggtitle("Shared significant DEGs between fluoxetine and citalopram contrasts")

  ggsave(file.path(out_fig, "SSRI_shared_significant_DEG_overlap_barplot.pdf"), p3, width = 13, height = 5)
  ggsave(file.path(out_fig, "SSRI_shared_significant_DEG_overlap_barplot.png"), p3, width = 13, height = 5, dpi = 300)
}

message("Cross-SSRI gene-level signature integration completed.")
