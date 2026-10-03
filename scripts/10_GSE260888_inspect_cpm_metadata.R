#!/usr/bin/env Rscript

project_dir <- "/mnt/d/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

cpm_file <- "03_raw_downloads/GSE260888_processed/GSE260888/GSE260888_cpm.csv.gz"
meta_file <- "02_metadata/GSE260888/GSE260888_series_matrix_metadata.csv"

out_dir <- "02_metadata/GSE260888"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading CPM matrix...")
cpm <- read.csv(gzfile(cpm_file), check.names = FALSE, stringsAsFactors = FALSE)

message("CPM dimensions:")
print(dim(cpm))

message("First 20 CPM columns:")
print(colnames(cpm)[1:min(20, ncol(cpm))])

message("Last 20 CPM columns:")
print(tail(colnames(cpm), 20))

write.csv(
  data.frame(column_number = seq_along(colnames(cpm)), column_name = colnames(cpm)),
  file.path(out_dir, "GSE260888_cpm_column_names.csv"),
  row.names = FALSE
)

message("Preview CPM first 5 rows and first 8 columns:")
print(cpm[1:min(5, nrow(cpm)), 1:min(8, ncol(cpm))])

message("Reading GEO metadata...")
meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)

message("Metadata dimensions:")
print(dim(meta))

message("Metadata columns:")
print(colnames(meta))

write.csv(
  data.frame(column_number = seq_along(colnames(meta)), column_name = colnames(meta)),
  file.path(out_dir, "GSE260888_metadata_column_names.csv"),
  row.names = FALSE
)

message("Preview metadata first 5 rows and selected columns:")
print(meta[1:min(5, nrow(meta)), 1:min(10, ncol(meta))])

# Try to identify GSM/sample columns in expression matrix
sample_cols <- grep("^GSM", colnames(cpm), value = TRUE)

message("Detected GSM sample columns in CPM:")
print(length(sample_cols))
print(head(sample_cols, 20))

write.csv(
  data.frame(sample_id = sample_cols),
  file.path(out_dir, "GSE260888_detected_cpm_sample_columns.csv"),
  row.names = FALSE
)

# Try to identify GSM/sample ids in metadata
gsm_like_cols <- names(meta)[sapply(meta, function(x) any(grepl("^GSM", as.character(x))))]

message("Metadata columns containing GSM-like values:")
print(gsm_like_cols)

if (length(gsm_like_cols) > 0) {
  for (cc in gsm_like_cols) {
    tmp <- meta[, cc]
    tmp <- tmp[grepl("^GSM", as.character(tmp))]
    message("Column: ", cc)
    print(head(tmp, 20))
  }
}

message("Inspection completed.")
