
# ============================================================
# GSE260888 - edgeR quasi-likelihood day-matched reanalysis
# PLOS ONE Major Revision - Reviewer Comment 2
# ============================================================

options(stringsAsFactors = FALSE)

project_dir <- "/mnt/g/D/SSRI_neurodevelopment_transcriptomics"
setwd(project_dir)

rev_dir <- "20_PLOS_revision/comment2_GSE260888_edgeR"
raw_dir <- file.path(rev_dir, "raw_featurecounts")
all_dir <- file.path(rev_dir, "results/all_genes")
sig_dir <- file.path(rev_dir, "results/significant")
rob_dir <- file.path(rev_dir, "results/robustness_vs_old_limma")
fig_dir <- file.path(rev_dir, "figures")
log_dir <- file.path(rev_dir, "logs")

dirs <- c(all_dir, sig_dir, rob_dir, fig_dir, log_dir)
for (d in dirs) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# Package check
# ------------------------------------------------------------

if (!requireNamespace("edgeR", quietly=TRUE)) {
    if (!requireNamespace("BiocManager", quietly=TRUE)) {
        install.packages(
            "BiocManager",
            repos="https://cloud.r-project.org"
        )
    }

    BiocManager::install(
        "edgeR",
        ask=FALSE,
        update=FALSE
    )
}

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

cat("R version:", R.version.string, "\n")
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
cat("limma version:", as.character(packageVersion("limma")), "\n\n")


# ------------------------------------------------------------
# 1. Read all featureCounts files
# ------------------------------------------------------------

fc_files <- sort(
    list.files(
        raw_dir,
        pattern="_hisat2_FC\\.txt\\.gz$",
        full.names=TRUE
    )
)

cat("Number of featureCounts files:", length(fc_files), "\n")

if (length(fc_files) != 89) {
    stop(
        paste(
            "Expected 89 featureCounts files but found",
            length(fc_files)
        )
    )
}


read_featurecounts <- function(f) {

    x <- read.delim(
        gzfile(f),
        header=TRUE,
        comment.char="#",
        check.names=FALSE,
        stringsAsFactors=FALSE
    )

    if (!"Geneid" %in% colnames(x)) {
        stop(
            paste(
                "Geneid column not found in",
                basename(f),
                "Columns:",
                paste(colnames(x), collapse=", ")
            )
        )
    }

    # Individual GEO featureCounts files contain one count column
    # after the feature annotation columns.
    count_col <- ncol(x)

    count_values <- suppressWarnings(
        as.numeric(x[[count_col]])
    )

    if (anyNA(count_values)) {
        stop(
            paste(
                "Non-numeric count values detected in",
                basename(f)
            )
        )
    }

    data.frame(
        gene_id=x$Geneid,
        count=count_values,
        stringsAsFactors=FALSE
    )
}


first <- read_featurecounts(fc_files[1])
gene_ids <- first$gene_id

count_matrix <- matrix(
    0,
    nrow=length(gene_ids),
    ncol=length(fc_files),
    dimnames=list(gene_ids, NULL)
)


sample_ids <- character(length(fc_files))

for (i in seq_along(fc_files)) {

    f <- fc_files[i]
    dat <- read_featurecounts(f)

    if (!identical(dat$gene_id, gene_ids)) {
        stop(
            paste(
                "Gene ordering differs in",
                basename(f)
            )
        )
    }

    count_matrix[, i] <- dat$count

    sample_ids[i] <- sub(
        "_hisat2_FC\\.txt\\.gz$",
        "",
        sub(
            "^GSM[0-9]+_",
            "",
            basename(f)
        )
    )
}

colnames(count_matrix) <- sample_ids

cat(
    "Count matrix dimensions:",
    nrow(count_matrix),
    "genes x",
    ncol(count_matrix),
    "samples\n"
)

cat(
    "Counts integer-valued:",
    all(count_matrix == round(count_matrix)),
    "\n\n"
)


# ------------------------------------------------------------
# 2. Metadata
# ------------------------------------------------------------

meta_file <- paste0(
    "08_processed_expression/GSE260888/",
    "GSE260888_parsed_sample_metadata.csv"
)

meta <- read.csv(
    meta_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

# Match raw featureCounts files to CPM-derived metadata using the
# unique numeric sample identifier rather than the complete sample name.
# This avoids a known CPM-header naming typo for sample 291
# (D13-CIT40-D in the CPM header versus D13-CIT50-D in the raw file).

meta$sample_number <- as.character(meta$sample_number)

raw_sample_ids <- colnames(count_matrix)

raw_sample_numbers <- sub(
    "-.*$",
    "",
    raw_sample_ids
)

if (anyDuplicated(meta$sample_number)) {
    stop("Duplicate sample_number values detected in metadata.")
}

if (anyDuplicated(raw_sample_numbers)) {
    stop("Duplicate sample numbers detected in raw count filenames.")
}

if (!all(raw_sample_numbers %in% meta$sample_number)) {

    missing_numbers <- setdiff(
        raw_sample_numbers,
        meta$sample_number
    )

    stop(
        paste(
            "Raw-count sample numbers missing from metadata:",
            paste(missing_numbers, collapse=", ")
        )
    )
}

meta <- meta[
    match(
        raw_sample_numbers,
        meta$sample_number
    ),
    ,
    drop=FALSE
]

# Preserve both identifiers for auditability
meta$cpm_sample_id <- meta$sample_id
meta$raw_count_sample_id <- raw_sample_ids

# Verify biological metadata against the raw-count filenames
raw_parts <- do.call(
    rbind,
    strsplit(
        raw_sample_ids,
        "-",
        fixed=TRUE
    )
)

raw_day <- raw_parts[,2]
raw_treatment <- raw_parts[,3]

if (!all(raw_day == meta$day_group)) {

    bad <- which(
        raw_day != meta$day_group
    )

    stop(
        paste(
            "Day mismatch between raw files and metadata at:",
            paste(raw_sample_ids[bad], collapse=", ")
        )
    )
}

if (!all(raw_treatment == as.character(meta$treatment))) {

    bad <- which(
        raw_treatment != as.character(meta$treatment)
    )

    stop(
        paste(
            "Treatment mismatch between raw files and metadata at:",
            paste(raw_sample_ids[bad], collapse=", ")
        )
    )
}

cat(
    "Metadata matched to raw featureCounts by unique sample number: PASS\n"
)

cat(
    "Raw filename versus CPM-derived metadata identifier differences:\n"
)

id_diff <- data.frame(
    sample_number=meta$sample_number,
    raw_count_sample_id=meta$raw_count_sample_id,
    cpm_sample_id=meta$cpm_sample_id,
    stringsAsFactors=FALSE
)

id_diff <- id_diff[
    id_diff$raw_count_sample_id != id_diff$cpm_sample_id,
    ,
    drop=FALSE
]

print(id_diff)

write.csv(
    id_diff,
    file.path(
        rev_dir,
        "GSE260888_sample_identifier_discrepancies.csv"
    ),
    row.names=FALSE
)

cat("\n")

cat("Exact sample numbers:\n")
print(
    with(
        meta,
        table(
            day_group,
            treatment
        )
    )
)


# ------------------------------------------------------------
# 3. Restrict treatment analysis to D6, D10 and D13
# ------------------------------------------------------------

keep_samples <- meta$day_group %in% c(
    "D6",
    "D10",
    "D13"
)

counts_de <- count_matrix[, keep_samples, drop=FALSE]
meta_de <- meta[keep_samples, , drop=FALSE]

meta_de$group <- paste(
    meta_de$day_group,
    meta_de$treatment,
    sep="_"
)

desired_groups <- c(
    "D6_CTRL",
    "D6_CIT50",
    "D6_CIT100",
    "D6_CIT200",
    "D6_CIT400",

    "D10_CTRL",
    "D10_CIT50",
    "D10_CIT100",
    "D10_CIT200",
    "D10_CIT400",

    "D13_CTRL",
    "D13_CIT50",
    "D13_CIT100",
    "D13_CIT200",
    "D13_CIT400"
)

meta_de$group <- factor(
    meta_de$group,
    levels=desired_groups
)

if (anyNA(meta_de$group)) {
    stop("Unexpected day/treatment group detected.")
}


# ------------------------------------------------------------
# 4. DGEList + low-expression filtering + TMM
# ------------------------------------------------------------

design <- model.matrix(
    ~ 0 + group,
    data=meta_de
)

colnames(design) <- levels(meta_de$group)

y <- DGEList(
    counts=counts_de,
    samples=meta_de
)

genes_before <- nrow(y)

keep_genes <- filterByExpr(
    y,
    design=design
)

y <- y[
    keep_genes,
    ,
    keep.lib.sizes=FALSE
]

genes_after <- nrow(y)

cat(
    "\nGenes before filterByExpr:",
    genes_before,
    "\n"
)

cat(
    "Genes retained after filterByExpr:",
    genes_after,
    "\n"
)

cat(
    "Genes removed:",
    genes_before - genes_after,
    "\n\n"
)

y <- calcNormFactors(
    y,
    method="TMM"
)

write.csv(
    data.frame(
        sample_id=colnames(y),
        day_group=meta_de$day_group,
        treatment=meta_de$treatment,
        dose_nM=meta_de$dose_nM,
        replicate=meta_de$replicate,
        library_size=y$samples$lib.size,
        norm_factor=y$samples$norm.factors
    ),
    file.path(
        rev_dir,
        "GSE260888_edgeR_sample_and_normalization_metadata.csv"
    ),
    row.names=FALSE
)


# ------------------------------------------------------------
# 5. edgeR quasi-likelihood fit
# ------------------------------------------------------------

y <- estimateDisp(
    y,
    design
)

fit <- glmQLFit(
    y,
    design
)

cat(
    "edgeR dispersion estimation: completed\n"
)

cat(
    "edgeR quasi-likelihood fitting: completed\n\n"
)


# ------------------------------------------------------------
# 6. Annotation
# ------------------------------------------------------------

annotation_file <- paste0(
    "08_processed_expression/GSE260888/",
    "GSE260888_gene_annotation_Ensembl_to_symbol.csv"
)

annotation <- read.csv(
    annotation_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

annotation$ensembl_clean <- sub(
    "\\..*$",
    "",
    annotation$ensembl_id
)

annotation <- annotation[
    !duplicated(annotation$ensembl_clean),
]

symbol_lookup <- setNames(
    annotation$gene_symbol,
    annotation$ensembl_clean
)


# ------------------------------------------------------------
# 7. Day-matched contrasts
# ------------------------------------------------------------

days <- c(
    "D6",
    "D10",
    "D13"
)

doses <- c(
    "CIT50",
    "CIT100",
    "CIT200",
    "CIT400"
)

summary_results <- data.frame()

result_cache <- list()


for (day_i in days) {

    for (dose_i in doses) {

        contrast_name <- paste0(
            day_i,
            "_",
            dose_i,
            "_vs_CTRL"
        )

        treated_group <- paste0(
            day_i,
            "_",
            dose_i
        )

        control_group <- paste0(
            day_i,
            "_CTRL"
        )

        contrast_vector <- makeContrasts(
            contrasts=paste0(
                treated_group,
                "-",
                control_group
            ),
            levels=design
        )

        qlf <- glmQLFTest(
            fit,
            contrast=contrast_vector
        )

        res <- topTags(
            qlf,
            n=Inf,
            sort.by="none"
        )$table

        res$gene_id <- rownames(res)

        res$ensembl_id <- sub(
            "\\..*$",
            "",
            res$gene_id
        )

        res$gene_symbol <- unname(
            symbol_lookup[
                res$ensembl_id
            ]
        )

        res$significant_FDR05 <- (
            !is.na(res$FDR) &
            res$FDR < 0.05
        )

        res$significant_FDR05_log2FC058 <- (
            !is.na(res$FDR) &
            res$FDR < 0.05 &
            abs(res$logFC) >= 0.58
        )

        res <- res[
            ,
            c(
                "gene_id",
                "ensembl_id",
                "gene_symbol",
                "logFC",
                "logCPM",
                "F",
                "PValue",
                "FDR",
                "significant_FDR05",
                "significant_FDR05_log2FC058"
            )
        ]

        # Order by statistical evidence for output
        res_sorted <- res[
            order(
                res$FDR,
                res$PValue
            ),
        ]

        write.csv(
            res_sorted,
            file.path(
                all_dir,
                paste0(
                    contrast_name,
                    "_edgeR_QL_all_genes.csv"
                )
            ),
            row.names=FALSE
        )

        sig <- res_sorted[
            res_sorted$significant_FDR05_log2FC058,
            ,
            drop=FALSE
        ]

        write.csv(
            sig,
            file.path(
                sig_dir,
                paste0(
                    contrast_name,
                    "_edgeR_QL_significant_FDR05_log2FC058.csv"
                )
            ),
            row.names=FALSE
        )

        n_up <- sum(
            sig$logFC > 0,
            na.rm=TRUE
        )

        n_down <- sum(
            sig$logFC < 0,
            na.rm=TRUE
        )

        n_fdr_only <- sum(
            res$significant_FDR05,
            na.rm=TRUE
        )

        n_ctrl <- sum(
            meta_de$day_group == day_i &
            meta_de$treatment == "CTRL"
        )

        n_treat <- sum(
            meta_de$day_group == day_i &
            meta_de$treatment == dose_i
        )

        summary_results <- rbind(
            summary_results,
            data.frame(
                contrast=contrast_name,
                day=day_i,
                treatment=dose_i,
                dose_nM=as.numeric(
                    sub("CIT", "", dose_i)
                ),
                n_control=n_ctrl,
                n_treated=n_treat,
                tested_genes=nrow(res),
                significant_FDR05=n_fdr_only,
                significant_FDR05_log2FC058=nrow(sig),
                upregulated=n_up,
                downregulated=n_down,
                stringsAsFactors=FALSE
            )
        )

        result_cache[[contrast_name]] <- res

        cat(
            contrast_name,
            " | CTRL n=",
            n_ctrl,
            " | treated n=",
            n_treat,
            " | FDR<0.05=",
            n_fdr_only,
            " | FDR<0.05 & |log2FC|>=0.58=",
            nrow(sig),
            " | UP=",
            n_up,
            " | DOWN=",
            n_down,
            "\n",
            sep=""
        )
    }
}


write.csv(
    summary_results,
    file.path(
        rev_dir,
        "GSE260888_edgeR_QL_DEG_summary.csv"
    ),
    row.names=FALSE
)


# ------------------------------------------------------------
# 8. Exact sample-size table
# ------------------------------------------------------------

sample_count_table <- as.data.frame(
    with(
        meta_de,
        table(
            day_group,
            treatment
        )
    )
)

colnames(sample_count_table) <- c(
    "day",
    "treatment",
    "n"
)

sample_count_table <- sample_count_table[
    sample_count_table$n > 0,
]

write.csv(
    sample_count_table,
    file.path(
        rev_dir,
        "GSE260888_exact_group_sample_sizes.csv"
    ),
    row.names=FALSE
)


# ------------------------------------------------------------
# 9. Compare corrected edgeR results with previous limma output
#    This is an internal robustness audit.
# ------------------------------------------------------------

old_de_dir <- "09_differential_expression/GSE260888"

robustness_summary <- data.frame()

candidate_genes <- c(
    "BDNF",
    "GABRA2",
    "THRB",
    "SEMA3D",
    "LGI1",
    "CDH13",
    "SH2D5",
    "SCN4A",
    "ATF5",
    "ASNS",
    "SESN2",
    "TRIB3",
    "SLC1A5",
    "ASS1",
    "PHGDH",
    "GPT2",
    "CHAC1",
    "PCK2",
    "CEBPB",
    "JDP2",
    "STC2",
    "TNFRSF12A",
    "IL11",
    "CRYAB",
    "SERPINB9",
    "EMP1"
)

candidate_comparison <- data.frame()


for (contrast_name in names(result_cache)) {

    old_file <- file.path(
        old_de_dir,
        paste0(
            contrast_name,
            "_all_genes.csv"
        )
    )

    if (!file.exists(old_file)) {
        warning(
            paste(
                "Old limma result not found:",
                old_file
            )
        )
        next
    }

    old <- read.csv(
        old_file,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )

    old$ensembl_id_clean <- sub(
        "\\..*$",
        "",
        old$ensembl_id
    )

    new <- result_cache[[contrast_name]]

    merged <- merge(
        new,
        old,
        by.x="ensembl_id",
        by.y="ensembl_id_clean",
        suffixes=c("_edgeR", "_limma")
    )

    # Identify old significance using manuscript rule
    old_sig <- (
        !is.na(merged$adj.P.Val) &
        merged$adj.P.Val < 0.05 &
        abs(merged$logFC_limma) >= 0.58
    )

    new_sig <- (
        !is.na(merged$FDR) &
        merged$FDR < 0.05 &
        abs(merged$logFC_edgeR) >= 0.58
    )

    both_sig <- old_sig & new_sig

    pearson_r <- suppressWarnings(
        cor(
            merged$logFC_edgeR,
            merged$logFC_limma,
            method="pearson",
            use="complete.obs"
        )
    )

    spearman_r <- suppressWarnings(
        cor(
            merged$logFC_edgeR,
            merged$logFC_limma,
            method="spearman",
            use="complete.obs"
        )
    )

    direction_agreement_all <- mean(
        sign(merged$logFC_edgeR) ==
        sign(merged$logFC_limma),
        na.rm=TRUE
    )

    direction_agreement_both_sig <- if (
        sum(both_sig) > 0
    ) {
        mean(
            sign(
                merged$logFC_edgeR[both_sig]
            ) ==
            sign(
                merged$logFC_limma[both_sig]
            ),
            na.rm=TRUE
        )
    } else {
        NA_real_
    }

    robustness_summary <- rbind(
        robustness_summary,
        data.frame(
            contrast=contrast_name,
            common_genes=nrow(merged),
            pearson_logFC=pearson_r,
            spearman_logFC=spearman_r,
            old_limma_significant=sum(old_sig),
            new_edgeR_significant=sum(new_sig),
            significant_in_both=sum(both_sig),
            direction_agreement_all=direction_agreement_all,
            direction_agreement_both_significant=direction_agreement_both_sig,
            stringsAsFactors=FALSE
        )
    )

    gene_col <- if (
        "gene_symbol_edgeR" %in% colnames(merged)
    ) {
        "gene_symbol_edgeR"
    } else {
        "gene_symbol"
    }

    candidate_rows <- merged[
        merged[[gene_col]] %in% candidate_genes,
        ,
        drop=FALSE
    ]

    if (nrow(candidate_rows) > 0) {

        candidate_rows$contrast <- contrast_name

        candidate_comparison <- rbind(
            candidate_comparison,
            candidate_rows
        )
    }
}


write.csv(
    robustness_summary,
    file.path(
        rob_dir,
        "edgeR_vs_previous_limma_summary.csv"
    ),
    row.names=FALSE
)

write.csv(
    candidate_comparison,
    file.path(
        rob_dir,
        "candidate_cross_SSRI_genes_edgeR_vs_previous_limma.csv"
    ),
    row.names=FALSE
)


# ------------------------------------------------------------
# 10. Save normalized logCPM matrix
# ------------------------------------------------------------

logcpm <- cpm(
    y,
    log=TRUE,
    prior.count=2
)

write.csv(
    data.frame(
        gene_id=rownames(logcpm),
        logcpm,
        check.names=FALSE
    ),
    file.path(
        rev_dir,
        "GSE260888_edgeR_TMM_logCPM_filtered.csv"
    ),
    row.names=FALSE
)


# ------------------------------------------------------------
# 11. MDS diagnostic plot
# ------------------------------------------------------------

pdf(
    file.path(
        fig_dir,
        "GSE260888_edgeR_MDS.pdf"
    ),
    width=9,
    height=7
)

plotMDS(
    y,
    labels=paste(
        meta_de$day_group,
        meta_de$treatment,
        sep="_"
    )
)

dev.off()


png(
    file.path(
        fig_dir,
        "GSE260888_edgeR_MDS.png"
    ),
    width=2400,
    height=1900,
    res=300
)

plotMDS(
    y,
    labels=paste(
        meta_de$day_group,
        meta_de$treatment,
        sep="_"
    )
)

dev.off()


# ------------------------------------------------------------
# 12. Session information
# ------------------------------------------------------------

sink(
    file.path(
        log_dir,
        "GSE260888_edgeR_sessionInfo.txt"
    )
)

sessionInfo()

sink()


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("GSE260888 edgeR QL REANALYSIS COMPLETED SUCCESSFULLY\n")
cat("====================================================\n\n")

cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filterByExpr:", genes_after, "\n\n")

cat("FINAL DEG SUMMARY\n")
print(
    summary_results,
    row.names=FALSE
)

cat("\nROBUSTNESS VS PREVIOUS LIMMA\n")
print(
    robustness_summary,
    row.names=FALSE
)

cat("\nOutputs saved to:\n")
cat(
    normalizePath(
        rev_dir,
        mustWork=FALSE
    ),
    "\n"
)

