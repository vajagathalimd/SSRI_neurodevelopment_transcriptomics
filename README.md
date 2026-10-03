# SSRI neurodevelopment transcriptomics

This repository contains the analysis scripts, processed outputs, figure source material, revision analyses, and submission-support files for the manuscript:

**Integrative transcriptomic analysis of SSRI exposure datasets identifies neurodevelopmental and stress-response signatures in human neural differentiation systems**

## Public datasets

- Fluoxetine dataset: GEO accession **GSE166297**
- Citalopram dataset: GEO accession **GSE260888**

The repository does not duplicate raw FASTQ files or large workflow intermediates. These files are available from GEO or can be regenerated using the included scripts and inputs.

## Main contents

- `scripts/` - original analysis scripts for processing, differential expression, enrichment, module scoring, figures, and submission organization.
- `metadata/` and `nfcore_inputs/` - sample metadata and nf-core/rnaseq input files.
- `processed_expression/` - processed expression matrices used for downstream analyses.
- `differential_expression/` - differential-expression output tables.
- `pathway_enrichment/`, `gene_modules/`, `cross_dataset_integration/`, and `validation/` - downstream analysis outputs.
- `tables/` - manuscript and supplementary table source outputs.
- `figures/` - manuscript figure files.
- `plos_revision/` - revision-specific analyses for the PLOS ONE response, including the edgeR citalopram reanalysis, formal cross-SSRI overlap tests, and module random-set benchmark.
- `submission_documents/` - revised manuscript, tracked manuscript, response to reviewers, and supplementary files.
- `software_logs/` - available software/session and workflow logs.

## Key software versions reported in the manuscript

- nf-core/rnaseq v3.26.0
- Nextflow v26.04.1
- Salmon v1.10.3
- DESeq2 v1.42.1
- R v4.3.3
- edgeR v4.0.16
- limma v3.58.1
- g:Profiler API using organism `hsapiens`, annotated-domain background, and g:SCS multiple-testing correction
- GRCh38 NCBI iGenomes reference bundle used by nf-core/rnaseq

## Notes for PLOS ONE resubmission

This folder is prepared as the GitHub repository content for review and is available at:

https://github.com/vajagathalimd/SSRI_neurodevelopment_transcriptomics

The manuscript Data Availability statement and response letter cite this public GitHub repository directly.
