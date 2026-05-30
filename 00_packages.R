required_packages <- c(
  "Seurat", "SeuratObject", "tidyverse", "data.table", "Matrix", "yaml",
  "GEOquery", "Biobase", "SingleCellExperiment", "scater", "scran", "scDblFinder",
  "edgeR", "DESeq2", "limma", "fgsea", "msigdbr", "clusterProfiler", "org.Hs.eg.db",
  "AnnotationDbi", "BiocParallel", "patchwork", "ggrepel", "pheatmap",
  "targets", "tarchetypes", "janitor", "harmony", "scales"
)

load_required_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      "\nRun source('scripts/install_packages.R') first.",
      call. = FALSE
    )
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

load_required_packages()
