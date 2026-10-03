#!/usr/bin/env bash

set -u

BASE="/mnt/d/SSRI_neurodevelopment_transcriptomics"
cd "$BASE" || exit 1

OUT="19_submission_package"
MF="$OUT/main_figures"
SF="$OUT/supplementary_figures"
MT="$OUT/main_tables"
ST="$OUT/supplementary_tables"

mkdir -p "$MF" "$SF" "$MT" "$ST"

copy_one () {
  src="$1"
  dest="$2"
  if [ -f "$src" ]; then
    cp "$src" "$dest"
    echo "Copied: $src -> $dest"
  else
    echo "MISSING: $src"
  fi
}

copy_glob () {
  pattern="$1"
  destdir="$2"
  mkdir -p "$destdir"
  shopt -s nullglob
  files=( $pattern )
  if [ ${#files[@]} -eq 0 ]; then
    echo "MISSING GLOB: $pattern"
  else
    for f in "${files[@]}"; do
      cp "$f" "$destdir/"
      echo "Copied: $f -> $destdir/"
    done
  fi
  shopt -u nullglob
}

echo "=============================="
echo "Creating main figure folder"
echo "=============================="

copy_one "14_figures/manuscript_combined_figures/Figure_1_GSE166297_workflow.pdf" "$MF/Figure_1_workflow.pdf"
copy_one "14_figures/manuscript_combined_figures/Figure_1_GSE166297_workflow.png" "$MF/Figure_1_workflow.png"

copy_one "14_figures/manuscript_combined_figures/Figure_2_GSE166297_DEG_overview.pdf" "$MF/Figure_2_fluoxetine_DEG_overview.pdf"
copy_one "14_figures/manuscript_combined_figures/Figure_2_GSE166297_DEG_overview.png" "$MF/Figure_2_fluoxetine_DEG_overview.png"

copy_one "14_figures/manuscript_combined_figures/Figure_3_GSE166297_enrichment_overview.pdf" "$MF/Figure_3_fluoxetine_enrichment_overview.pdf"
copy_one "14_figures/manuscript_combined_figures/Figure_3_GSE166297_enrichment_overview.png" "$MF/Figure_3_fluoxetine_enrichment_overview.png"

copy_one "14_figures/manuscript_combined_figures/Figure_4_GSE166297_neurodevelopmental_module_shifts.pdf" "$MF/Figure_4_fluoxetine_module_shifts.pdf"
copy_one "14_figures/manuscript_combined_figures/Figure_4_GSE166297_neurodevelopmental_module_shifts.png" "$MF/Figure_4_fluoxetine_module_shifts.png"

copy_one "14_figures/manuscript_combined_figures/Figure_5_cross_SSRI_comparative_integration.pdf" "$MF/Figure_5_cross_SSRI_comparative_integration.pdf"
copy_one "14_figures/manuscript_combined_figures/Figure_5_cross_SSRI_comparative_integration.png" "$MF/Figure_5_cross_SSRI_comparative_integration.png"


echo "=============================="
echo "Creating supplementary figure folder"
echo "=============================="

copy_one "14_figures/manuscript_combined_figures/Supplementary_Figure_GSE166297_sample_level_module_heatmap.pdf" "$SF/Supplementary_Figure_1_GSE166297_sample_level_module_heatmap.pdf"
copy_one "14_figures/manuscript_combined_figures/Supplementary_Figure_GSE166297_sample_level_module_heatmap.png" "$SF/Supplementary_Figure_1_GSE166297_sample_level_module_heatmap.png"

copy_one "14_figures/GSE260888/GSE260888_PCA_all_samples.pdf" "$SF/Supplementary_Figure_2_GSE260888_PCA_all_samples.pdf"
copy_one "14_figures/GSE260888/GSE260888_PCA_all_samples.png" "$SF/Supplementary_Figure_2_GSE260888_PCA_all_samples.png"

mkdir -p "$SF/Supplementary_Figure_3_GSE260888_citalopram_volcano_plots"
copy_glob "14_figures/GSE260888/*_volcano.pdf" "$SF/Supplementary_Figure_3_GSE260888_citalopram_volcano_plots"
copy_glob "14_figures/GSE260888/*_volcano.png" "$SF/Supplementary_Figure_3_GSE260888_citalopram_volcano_plots"

copy_one "14_figures/GSE260888/GSE260888_condition_mean_module_heatmap.pdf" "$SF/Supplementary_Figure_4_GSE260888_condition_mean_module_heatmap.pdf"
copy_one "14_figures/GSE260888/GSE260888_condition_mean_module_heatmap.png" "$SF/Supplementary_Figure_4_GSE260888_condition_mean_module_heatmap.png"

copy_one "14_figures/cross_SSRI_gene_level/SSRI_gene_level_pearson_correlation_heatmap.pdf" "$SF/Supplementary_Figure_5_cross_SSRI_pearson_correlation_heatmap.pdf"
copy_one "14_figures/cross_SSRI_gene_level/SSRI_gene_level_pearson_correlation_heatmap.png" "$SF/Supplementary_Figure_5_cross_SSRI_pearson_correlation_heatmap.png"


echo "=============================="
echo "Creating main tables"
echo "=============================="

cat > "$MT/Table_1_dataset_characteristics_and_analysis_design.csv" <<'CSV'
dataset,compound,model_system,groups_or_timepoints,sample_size,input_data,processing_method,differential_expression_method,role_in_study
GSE166297,Fluoxetine,Human neural progenitor-derived neuron-astrocyte differentiation cultures,"Control_D10, FLX_IC5, FLX_IC20_100","n=8 per group; total n=24",Raw RNA-seq FASTQ,"nf-core/rnaseq; fastp; FastQC; Salmon; GRCh38",DESeq2,Primary fluoxetine dataset
GSE260888,Citalopram,Human neural differentiation cultures,"D6, D10, D13 day-matched controls and CIT50/CIT100/CIT200/CIT400; D0 controls included","D6: 6 per group; D10 and D13: mostly 5-6 per group",Processed CPM matrix,"Parsed CPM matrix; Ensembl-to-symbol mapping; logCPM-based analysis",limma,Independent citalopram comparison dataset
CSV

copy_one "15_tables/cross_SSRI_integration/SSRI_combined_DEG_count_summary.csv" "$MT/Table_2_differential_expression_summary_across_SSRI_contrasts.csv"

copy_one "15_tables/GSE166297_final_summary/GSE166297_neurodevelopmental_module_statistics.csv" "$MT/Table_3_fluoxetine_neurodevelopmental_module_statistics.csv"

copy_one "15_tables/cross_SSRI_gene_level/SSRI_shared_downregulated_gene_recurrence_summary.csv" "$MT/Table_4_recurrent_shared_downregulated_genes_across_SSRI_contrasts.csv"


echo "=============================="
echo "Creating supplementary tables"
echo "=============================="

copy_one "02_metadata/GSE166297_pilot_metadata.csv" "$ST/Supplementary_Table_1_GSE166297_curated_sample_metadata.csv"

mkdir -p "$ST/Supplementary_Table_2_GSE166297_full_DESeq2_results"
copy_one "09_differential_expression/GSE166297/FLX_IC5_vs_Control_D10_all_genes.csv" "$ST/Supplementary_Table_2_GSE166297_full_DESeq2_results/FLX_IC5_vs_Control_D10_all_genes.csv"
copy_one "09_differential_expression/GSE166297/FLX_IC20_100_vs_Control_D10_all_genes.csv" "$ST/Supplementary_Table_2_GSE166297_full_DESeq2_results/FLX_IC20_100_vs_Control_D10_all_genes.csv"

mkdir -p "$ST/Supplementary_Table_3_GSE166297_significant_DEGs"
copy_one "09_differential_expression/GSE166297/FLX_IC5_vs_Control_D10_significant_FDR05_log2FC058.csv" "$ST/Supplementary_Table_3_GSE166297_significant_DEGs/FLX_IC5_vs_Control_D10_significant_DEGs.csv"
copy_one "09_differential_expression/GSE166297/FLX_IC20_100_vs_Control_D10_significant_FDR05_log2FC058.csv" "$ST/Supplementary_Table_3_GSE166297_significant_DEGs/FLX_IC20_100_vs_Control_D10_significant_DEGs.csv"

mkdir -p "$ST/Supplementary_Table_4_GSE166297_DEG_overlap_direction"
copy_one "15_tables/GSE166297_final_summary/GSE166297_DEG_overlap_summary.csv" "$ST/Supplementary_Table_4_GSE166297_DEG_overlap_direction/GSE166297_DEG_overlap_summary.csv"
copy_one "15_tables/GSE166297_final_summary/GSE166297_DEG_direction_summary.csv" "$ST/Supplementary_Table_4_GSE166297_DEG_overlap_direction/GSE166297_DEG_direction_summary.csv"

mkdir -p "$ST/Supplementary_Table_5_GSE166297_gProfiler_enrichment"
copy_glob "11_pathway_enrichment/GSE166297_gprofiler/*_gprofiler_enrichment.csv" "$ST/Supplementary_Table_5_GSE166297_gProfiler_enrichment"
copy_glob "11_pathway_enrichment/GSE166297*/*_gprofiler_enrichment.csv" "$ST/Supplementary_Table_5_GSE166297_gProfiler_enrichment"

mkdir -p "$ST/Supplementary_Table_6_curated_module_gene_matching"
copy_one "15_tables/GSE166297_module_scores/GSE166297_module_gene_matching.csv" "$ST/Supplementary_Table_6_curated_module_gene_matching/GSE166297_module_gene_matching.csv"
copy_one "15_tables/GSE260888/GSE260888_module_gene_matching.csv" "$ST/Supplementary_Table_6_curated_module_gene_matching/GSE260888_module_gene_matching.csv"

copy_one "15_tables/GSE166297_final_summary/GSE166297_neurodevelopmental_module_statistics.csv" "$ST/Supplementary_Table_7_GSE166297_full_module_statistics.csv"

copy_one "08_processed_expression/GSE260888/GSE260888_parsed_sample_metadata.csv" "$ST/Supplementary_Table_8_GSE260888_parsed_sample_metadata.csv"

mkdir -p "$ST/Supplementary_Table_9_GSE260888_full_limma_results"
copy_glob "09_differential_expression/GSE260888/*_all_genes.csv" "$ST/Supplementary_Table_9_GSE260888_full_limma_results"

mkdir -p "$ST/Supplementary_Table_10_GSE260888_significant_DEGs_and_summary"
copy_one "15_tables/GSE260888/GSE260888_limma_DEG_summary.csv" "$ST/Supplementary_Table_10_GSE260888_significant_DEGs_and_summary/GSE260888_limma_DEG_summary.csv"
copy_glob "09_differential_expression/GSE260888/*_significant_FDR05_log2FC058.csv" "$ST/Supplementary_Table_10_GSE260888_significant_DEGs_and_summary"

copy_one "15_tables/GSE260888/GSE260888_module_statistics_vs_day_matched_control.csv" "$ST/Supplementary_Table_11_GSE260888_module_statistics_vs_day_matched_control.csv"

mkdir -p "$ST/Supplementary_Table_12_cross_SSRI_module_shift_comparison"
copy_one "15_tables/cross_SSRI_integration/SSRI_combined_module_shift_table.csv" "$ST/Supplementary_Table_12_cross_SSRI_module_shift_comparison/SSRI_combined_module_shift_table.csv"
copy_one "15_tables/cross_SSRI_integration/SSRI_focused_D10_D13_module_shift_table.csv" "$ST/Supplementary_Table_12_cross_SSRI_module_shift_comparison/SSRI_focused_D10_D13_module_shift_table.csv"

copy_one "15_tables/cross_SSRI_gene_level/SSRI_gene_level_logFC_correlation_summary.csv" "$ST/Supplementary_Table_13_cross_SSRI_gene_level_logFC_correlation_summary.csv"

copy_one "15_tables/cross_SSRI_gene_level/SSRI_significant_DEG_overlap_summary.csv" "$ST/Supplementary_Table_14_shared_significant_DEG_overlap_across_SSRI_contrasts.csv"

copy_one "15_tables/cross_SSRI_gene_level/SSRI_shared_downregulated_gene_recurrence_summary.csv" "$ST/Supplementary_Table_15_recurrent_shared_downregulated_SSRI_genes.csv"


echo "=============================="
echo "Creating manifest"
echo "=============================="

cat > "$OUT/README_submission_file_manifest.txt" <<'MANIFEST'
Submission package folders

main_figures:
Figure 1. Transcriptomic reanalysis workflow.
Figure 2. Fluoxetine DEG overview.
Figure 3. Fluoxetine enrichment overview.
Figure 4. Fluoxetine neurodevelopmental module shifts.
Figure 5. Cross-SSRI comparative integration.

supplementary_figures:
Supplementary Figure 1. GSE166297 sample-level module heatmap.
Supplementary Figure 2. GSE260888 PCA.
Supplementary Figure 3. GSE260888 citalopram volcano plots.
Supplementary Figure 4. GSE260888 condition-mean module heatmap.
Supplementary Figure 5. Cross-SSRI Pearson correlation heatmap.

main_tables:
Table 1. Dataset characteristics and analysis design.
Table 2. Differential expression summary across SSRI contrasts.
Table 3. Fluoxetine neurodevelopmental module statistics.
Table 4. Recurrent shared downregulated genes across SSRI contrasts.

supplementary_tables:
Supplementary Table 1. GSE166297 curated sample metadata.
Supplementary Table 2. Full DESeq2 results for GSE166297.
Supplementary Table 3. Significant DEGs from GSE166297.
Supplementary Table 4. GSE166297 DEG overlap and direction summary.
Supplementary Table 5. g:Profiler enrichment results for fluoxetine DEG sets.
Supplementary Table 6. Curated module gene matching.
Supplementary Table 7. GSE166297 module-score statistics.
Supplementary Table 8. GSE260888 parsed sample metadata.
Supplementary Table 9. Full limma results for GSE260888.
Supplementary Table 10. GSE260888 significant DEG summary.
Supplementary Table 11. GSE260888 module-score statistics.
Supplementary Table 12. Cross-SSRI module-shift comparison.
Supplementary Table 13. Cross-SSRI gene-level logFC correlation summary.
Supplementary Table 14. Shared significant DEG overlap across SSRI contrasts.
Supplementary Table 15. Recurrent shared downregulated SSRI-responsive genes.
MANIFEST

echo "=============================="
echo "Submission package created:"
echo "$OUT"
echo "=============================="

find "$OUT" -maxdepth 3 -type f | sort
