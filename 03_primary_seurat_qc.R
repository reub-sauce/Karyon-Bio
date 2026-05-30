# Primary GSE136103 scRNA-seq loading, QC, doublet removal, normalization, integration.

read_one_mtx_sample <- function(row) {
  if (is.na(row$barcode_path) || is.na(row$feature_path)) {
    stop("Missing barcode or feature file for sample ", row$sample_id, call. = FALSE)
  }
  message("Reading ", row$sample_id)
  counts <- Seurat::ReadMtx(
    mtx = row$matrix_path,
    cells = row$barcode_path,
    features = row$feature_path,
    feature.column = 2,
    unique.features = TRUE
  )

  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "GSE136103",
    min.cells = 3,
    min.features = 100
  )

  obj <- Seurat::RenameCells(obj, add.cell.id = row$sample_id)
  obj$geo_accession <- row$geo_accession
  obj$sample_id <- row$sample_id
  obj$sample_label <- row$sample_label
  obj$donor_id <- row$donor_id
  obj$disease_status <- as.character(row$disease_status)
  obj$fibrosis_stage_harmonized <- row$fibrosis_stage_harmonized
  obj$tissue <- row$tissue
  obj$species <- row$species
  obj$fraction <- row$fraction
  obj
}

read_gse136103_samples <- function(sample_table, config = read_config()) {
  primary_samples <- sample_table %>% filter(include_primary)
  if (nrow(primary_samples) == 0) stop("No primary human liver samples selected.", call. = FALSE)
  objs <- purrr::map(seq_len(nrow(primary_samples)), ~ read_one_mtx_sample(primary_samples[.x, ]))
  names(objs) <- primary_samples$sample_id
  objs
}

add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^HB[ABDEGQMZ]")
  obj
}

merge_seurat_list <- function(seurat_list) {
  if (length(seurat_list) == 1) return(ensure_joined_layers(seurat_list[[1]]))

  # Seurat v5 no longer exports merge() as Seurat::merge.
  # Use the S3 merge method instead, which dispatches correctly for Seurat objects
  # in both Seurat v4 and v5.
  merged <- seurat_list[[1]]
  if (length(seurat_list) > 1) {
    for (i in seq_along(seurat_list)[-1]) {
      merged <- merge(x = merged, y = seurat_list[[i]], project = "GSE136103")
    }
  }

  # Seurat v5 keeps one counts/data layer per sample after merge. Many downstream
  # functions used here, including Seurat::as.SingleCellExperiment() and
  # GetAssayData(), expect a single joined layer. JoinLayers() is a no-op for
  # Seurat v4 / non-layered assays, so this keeps the pipeline cross-version safe.
  merged <- ensure_joined_layers(merged)
  merged
}

ensure_joined_layers <- function(obj, assay = "RNA") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) return(obj)
  if (!assay %in% names(obj@assays)) return(obj)

  # Only Seurat v5 Assay5 objects have Layers(). If Layers() or JoinLayers() is not
  # available, return unchanged for Seurat v4 compatibility.
  layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  if (length(layers) <= 1) return(obj)

  message("Joining Seurat v5 assay layers for assay ", assay)
  obj <- tryCatch(
    SeuratObject::JoinLayers(obj, assay = assay),
    error = function(e) {
      if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
        get("JoinLayers", envir = asNamespace("Seurat"))(obj, assay = assay)
      } else {
        stop(e)
      }
    }
  )
  obj
}

calculate_qc_thresholds <- function(meta, config = read_config()) {
  qc <- config$qc
  meta %>%
    group_by(sample_id) %>%
    summarise(
      min_features = qc$min_features_hard,
      min_counts = qc$min_counts_hard,
      max_features = mad_upper(nFeature_RNA, n_mad = qc$max_features_mad, hard_max = Inf),
      max_counts = mad_upper(nCount_RNA, n_mad = qc$max_counts_mad, hard_max = Inf),
      max_percent_mt = mad_upper(percent.mt, n_mad = qc$max_percent_mt_mad, hard_max = qc$max_percent_mt_hard),
      max_percent_hb = qc$max_percent_hb_hard,
      .groups = "drop"
    )
}

plot_qc_violin <- function(obj, title = NULL) {
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb")
  suppressWarnings(
    Seurat::VlnPlot(obj, features = features, group.by = "sample_id", ncol = 2, pt.size = 0) +
      patchwork::plot_annotation(title = title %||% "QC metrics")
  )
}

plot_qc_scatter <- function(obj) {
  df <- obj@meta.data %>% tibble::rownames_to_column("cell")
  ggplot(df, aes(nCount_RNA, nFeature_RNA, color = disease_status)) +
    geom_point(alpha = 0.25, size = 0.2) +
    facet_wrap(~ sample_id, scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw(base_size = 9) +
    labs(x = "UMIs / counts", y = "Detected genes", color = "Disease", title = "GSE136103 QC: genes vs counts")
}

run_doublet_detection <- function(obj, config = read_config()) {
  if (!isTRUE(config$qc$remove_doublets)) {
    obj$doublet_class <- "not_run"
    obj$doublet_score <- NA_real_
    return(obj)
  }
  message_header("Running scDblFinder doublet detection")
  obj <- ensure_joined_layers(obj)
  sce <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")
  sce$sample_id <- obj$sample_id
  set.seed(123)
  bp <- BiocParallel::SerialParam()
  sce <- scDblFinder::scDblFinder(sce, samples = "sample_id", BPPARAM = bp)
  obj$doublet_score <- sce$scDblFinder.score
  obj$doublet_class <- as.character(sce$scDblFinder.class)
  obj
}

filter_by_qc <- function(obj, thresholds) {
  meta <- obj@meta.data %>% tibble::rownames_to_column("cell")
  meta2 <- meta %>% left_join(thresholds, by = "sample_id") %>%
    mutate(
      pass_feature_count = nFeature_RNA >= min_features & nFeature_RNA <= max_features,
      pass_umi_count = nCount_RNA >= min_counts & nCount_RNA <= max_counts,
      pass_mito = percent.mt <= max_percent_mt,
      pass_hb = percent.hb <= max_percent_hb,
      pass_doublet = is.na(doublet_class) | doublet_class != "doublet",
      pass_qc = pass_feature_count & pass_umi_count & pass_mito & pass_hb & pass_doublet
    )
  obj$pass_feature_count <- meta2$pass_feature_count[match(colnames(obj), meta2$cell)]
  obj$pass_umi_count <- meta2$pass_umi_count[match(colnames(obj), meta2$cell)]
  obj$pass_mito <- meta2$pass_mito[match(colnames(obj), meta2$cell)]
  obj$pass_hb <- meta2$pass_hb[match(colnames(obj), meta2$cell)]
  obj$pass_doublet <- meta2$pass_doublet[match(colnames(obj), meta2$cell)]
  obj$pass_qc <- meta2$pass_qc[match(colnames(obj), meta2$cell)]

  filtered <- subset(obj, subset = pass_qc)
  list(obj_with_flags = obj, filtered = filtered, cell_qc = meta2)
}

summarize_qc <- function(cell_qc, config = read_config()) {
  summary <- cell_qc %>%
    group_by(sample_id, donor_id, disease_status, fraction) %>%
    summarise(
      n_cells_prefilter = n(),
      n_cells_postfilter = sum(pass_qc, na.rm = TRUE),
      pct_retained = 100 * n_cells_postfilter / n_cells_prefilter,
      median_features_prefilter = median(nFeature_RNA, na.rm = TRUE),
      median_counts_prefilter = median(nCount_RNA, na.rm = TRUE),
      median_percent_mt_prefilter = median(percent.mt, na.rm = TRUE),
      n_doublets = sum(doublet_class == "doublet", na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safely(summary, file.path(config$paths$table_dir, "qc_summary_by_sample.csv"))
  summary
}

normalize_and_integrate <- function(obj, config = read_config()) {
  message_header("Normalizing and integrating primary data")
  obj <- ensure_joined_layers(obj)
  n_dims <- config$integration$dims
  variable_features <- config$normalization$variable_features

  obj$disease_status <- factor(obj$disease_status, levels = c("healthy", "cirrhotic"))
  obj <- Seurat::NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = variable_features, verbose = FALSE)
  obj <- Seurat::ScaleData(object = obj, vars.to.regress = unlist(config$normalization$regress_vars), verbose = FALSE)
  obj <- Seurat::RunPCA(object = obj, features = Seurat::VariableFeatures(obj), npcs = n_dims, verbose = FALSE)

  reduction_use <- "pca"
  if (identical(config$integration$method, "harmony") && requireNamespace("harmony", quietly = TRUE)) {
    batch_var <- config$integration$batch_var
    if (length(batch_var) == 1 && batch_var %in% colnames(obj@meta.data)) {
      message("Running Harmony batch correction using ", batch_var)
      obj <- tryCatch(
        {
          # Harmony argument names changed across package/Seurat versions.
          # Use fully named, current-compatible arguments and avoid ambiguous
          # partial matches such as reduction= or assay.use=.
          harmony::RunHarmony(
            object = obj,
            group.by.vars = batch_var,
            reduction.use = "pca",
            dims.use = seq_len(n_dims),
            reduction.save = "harmony",
            project.dim = FALSE,
            verbose = FALSE
          )
        },
        error = function(e) {
          message("Harmony failed; continuing with PCA-only workflow. Error: ", conditionMessage(e))
          obj
        }
      )
      if ("harmony" %in% names(obj@reductions)) {
        reduction_use <- "harmony"
      }
    } else {
      message("Harmony batch variable not found; continuing with PCA-only workflow.")
    }
  }

  obj <- Seurat::RunUMAP(object = obj, reduction = reduction_use, dims = seq_len(n_dims), verbose = FALSE)
  obj <- Seurat::FindNeighbors(object = obj, reduction = reduction_use, dims = seq_len(n_dims), verbose = FALSE)
  obj <- Seurat::FindClusters(object = obj, resolution = config$integration$cluster_resolution, verbose = FALSE)
  obj$primary_cluster <- as.character(Seurat::Idents(obj))
  obj@misc$integration_reduction <- reduction_use
  obj
}

plot_embedding_overview <- function(obj, config = read_config()) {
  p1 <- Seurat::DimPlot(obj, group.by = "disease_status", reduction = "umap") + ggtitle("Disease status")
  p2 <- Seurat::DimPlot(obj, group.by = "sample_id", reduction = "umap") + ggtitle("Sample") + theme(legend.position = "none")
  p3 <- Seurat::DimPlot(obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE) + ggtitle("Clusters")
  plot <- p1 + p2 + p3
  save_plot_safely(plot, file.path(config$paths$figure_dir, "umap_disease_sample_cluster.pdf"), width = 15, height = 5)
}

run_primary_qc_pipeline <- function(seurat_list, config = read_config()) {
  message_header("Running primary QC pipeline")
  seurat_list <- purrr::map(seurat_list, add_qc_metrics)
  merged <- merge_seurat_list(seurat_list)

  save_plot_safely(plot_qc_violin(merged, "Prefilter QC metrics"), file.path(config$paths$figure_dir, "qc_violin_prefilter.pdf"), width = 14, height = 8)
  save_plot_safely(plot_qc_scatter(merged), file.path(config$paths$figure_dir, "qc_scatter_nFeature_nCount.pdf"), width = 14, height = 10)

  thresholds <- calculate_qc_thresholds(merged@meta.data, config)
  write_csv_safely(thresholds, file.path(config$paths$table_dir, "qc_thresholds_by_sample.csv"))

  merged <- run_doublet_detection(merged, config)
  qc_res <- filter_by_qc(merged, thresholds)
  qc_summary <- summarize_qc(qc_res$cell_qc, config)
  filtered <- qc_res$filtered

  save_plot_safely(plot_qc_violin(filtered, "Postfilter QC metrics"), file.path(config$paths$figure_dir, "qc_violin_postfilter.pdf"), width = 14, height = 8)

  integrated <- normalize_and_integrate(filtered, config)
  plot_embedding_overview(integrated, config)
  save_rds_safely(integrated, file.path(config$paths$rds_dir, "gse136103_qc_integrated.rds"))

  list(
    seurat = integrated,
    qc_thresholds = thresholds,
    qc_summary = qc_summary,
    cell_qc = qc_res$cell_qc
  )
}
