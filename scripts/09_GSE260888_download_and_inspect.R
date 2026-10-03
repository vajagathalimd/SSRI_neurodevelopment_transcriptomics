#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GEOquery)
})

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

geo_id <- "GSE260888"

out_raw  <- "03_raw_downloads/GSE260888_processed"
out_meta <- "02_metadata/GSE260888"

dir.create(out_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(out_meta, recursive = TRUE, showWarnings = FALSE)

message("Downloading GEO series matrix for ", geo_id)
gse <- getGEO(geo_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

if (length(gse) > 1) {
  eset <- gse[[1]]
} else {
  eset <- gse[[1]]
}

pheno <- pData(eset)

write.csv(
  pheno,
  file.path(out_meta, "GSE260888_series_matrix_metadata.csv"),
  row.names = TRUE
)

message("Metadata dimensions:")
print(dim(pheno))

message("Metadata columns:")
print(colnames(pheno))

message("Downloading supplementary files for ", geo_id)
getGEOSuppFiles(geo_id, baseDir = out_raw, makeDirectory = TRUE)

message("Downloaded files:")
print(list.files(out_raw, recursive = TRUE, full.names = TRUE))

message("Done.")
