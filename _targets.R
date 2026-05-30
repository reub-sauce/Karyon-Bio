library(targets)
source("R/00_packages.R")
source("R/01_utils.R")
source("R/02_download_and_metadata.R")
source("R/03_primary_seurat_qc.R")
source("R/04_annotation.R")
source("R/05_pseudobulk_de.R")
source("R/06_pathway_analysis.R")
source("R/07_validation_gse207310.R")
source("R/08_biomarker_prioritization.R")
source("R/09_figures_and_reports.R")

tar_option_set(
  packages = c(
    "Seurat", "tidyverse", "data.table", "Matrix", "yaml", "GEOquery",
    "Biobase", "SingleCellExperiment", "scater", "scran", "scDblFinder", "edgeR",
    "DESeq2", "limma", "fgsea", "msigdbr", "clusterProfiler", "org.Hs.eg.db",
    "AnnotationDbi", "BiocParallel", "patchwork", "ggrepel", "pheatmap",
    "janitor", "harmony", "scales"
  ),
  error = "continue"
)

list(
  tar_target(config, read_config(), deployment = "main"),
  tar_target(initialized_project_dirs, init_project_dirs(config)),

  tar_target(primary_geo_metadata, get_geo_metadata(config$project$primary_geo, config$paths$table_dir)),
  tar_target(primary_downloaded_files, download_geo_supplementary(config$project$primary_geo, config$paths$raw_dir)),
  tar_target(primary_sample_table, build_gse136103_sample_table(primary_downloaded_files, primary_geo_metadata, config)),
  tar_target(primary_dataset_summary, summarize_primary_dataset(primary_sample_table, config)),
  tar_target(primary_seurat_list, read_gse136103_samples(primary_sample_table, config)),
  tar_target(primary_qc, run_primary_qc_pipeline(primary_seurat_list, config)),
  tar_target(primary_annotated, annotate_primary(primary_qc, config)),
  tar_target(primary_de, run_all_pseudobulk_de(primary_annotated, config)),
  tar_target(primary_pathway, run_pathway_analysis(primary_de, config)),

  tar_target(validation_geo_metadata, get_geo_metadata(config$project$validation_geo, config$paths$table_dir)),
  tar_target(validation_downloaded_files, download_geo_supplementary(config$project$validation_geo, config$paths$raw_dir)),
  tar_target(validation_result, run_gse207310_validation(validation_downloaded_files, validation_geo_metadata, primary_de, config)),

  tar_target(biomarker_result, prioritize_biomarkers(primary_annotated, primary_de, primary_pathway, validation_result, config)),
  tar_target(report_paths, render_all_reports(biomarker_result, config))
)
