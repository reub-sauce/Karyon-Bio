# Marker-based cell type annotation and required compartment validation.


# Use the project-wide sanitizer from R/01_utils.R. Redefine here as a
# source-order guard so later scripts also get the robust version.
sanitize_table <- function(x) {
  if (is.null(x)) return(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE)
  for (nm in names(x)) {
    col <- x[[nm]]
    if (is.data.frame(col)) {
      x[[nm]] <- apply(col, 1, function(v) paste(as.character(v), collapse = ";"))
    } else if (is.matrix(col) || is.array(col)) {
      if (length(dim(col)) > 1 && nrow(col) == nrow(x)) {
        x[[nm]] <- apply(col, 1, function(v) paste(as.character(v), collapse = ";"))
      } else {
        x[[nm]] <- as.vector(col)
      }
    } else if (is.list(col)) {
      x[[nm]] <- vapply(col, function(v) {
        if (length(v) == 0) return(NA_character_)
        if (all(is.na(v))) return(NA_character_)
        paste(as.character(unlist(v, use.names = FALSE)), collapse = ";")
      }, character(1))
    }
  }
  tibble::as_tibble(x)
}

sanitize_seurat_metadata <- function(obj) {
  md <- sanitize_table(obj@meta.data)
  md <- as.data.frame(md, stringsAsFactors = FALSE)
  rownames(md) <- colnames(obj)
  obj@meta.data <- md
  obj
}

load_marker_table <- function(path = "resources/cell_type_markers.csv") {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(gene = clean_gene_symbols(gene))
}

marker_list_by_cell_type <- function(marker_table = load_marker_table(), role = c("canonical", "state", "all")) {
  role <- match.arg(role)
  mt <- marker_table
  if (role != "all") mt <- mt %>% filter(marker_role == role)
  mt %>%
    group_by(cell_type) %>%
    summarise(markers = list(unique(gene)), .groups = "drop") %>%
    tibble::deframe()
}

add_marker_module_scores <- function(obj, marker_sets, prefix = "MS_") {
  genes_present <- rownames(obj)
  for (nm in names(marker_sets)) {
    genes <- intersect(marker_sets[[nm]], genes_present)
    if (length(genes) >= 2) {
      obj <- Seurat::AddModuleScore(obj, features = list(genes), name = paste0(prefix, nm, "_"), search = FALSE)
      new_col <- grep(paste0("^", prefix, nm, "_"), colnames(obj@meta.data), value = TRUE)
      colnames(obj@meta.data)[match(new_col, colnames(obj@meta.data))] <- paste0(prefix, nm)
    }
  }
  obj
}

cluster_annotation_scores <- function(obj, marker_sets, cluster_var = "seurat_clusters", prefix = "MS_") {
  obj <- sanitize_seurat_metadata(obj)
  score_cols <- paste0(prefix, names(marker_sets))
  score_cols <- intersect(score_cols, colnames(obj@meta.data))
  if (length(score_cols) == 0) {
    stop("No marker module-score columns were found for annotation. Check marker genes and object rownames.")
  }

  md <- sanitize_table(obj@meta.data)
  if (!cluster_var %in% names(md)) {
    stop("Cluster variable not found in Seurat metadata: ", cluster_var)
  }
  md[[cluster_var]] <- as.character(md[[cluster_var]])
  for (cc in score_cols) md[[cc]] <- suppressWarnings(as.numeric(md[[cc]]))

  md %>%
    dplyr::group_by(.data[[cluster_var]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "cell_type", values_to = "mean_score") %>%
    dplyr::mutate(
      cell_type = stringr::str_remove(cell_type, paste0("^", prefix)),
      mean_score = suppressWarnings(as.numeric(mean_score))
    ) %>%
    dplyr::group_by(.data[[cluster_var]]) %>%
    dplyr::mutate(rank_within_cluster = rank(-mean_score, ties.method = "first")) %>%
    dplyr::ungroup()
}


initial_annotation_from_scores <- function(scores, cluster_var = "seurat_clusters") {
  scores %>%
    filter(rank_within_cluster == 1) %>%
    transmute(
      cluster = as.character(.data[[cluster_var]]),
      predicted_cell_type = cell_type,
      best_marker_score = mean_score,
      n_cells = n_cells
    )
}

load_manual_annotations <- function(path = "config/manual_cluster_annotations.csv") {
  if (!file.exists(path)) return(tibble::tibble())
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (nrow(x) == 0) return(tibble::tibble())
  x %>% filter(!is.na(cluster), cluster != "", !is.na(cell_type), cell_type != "") %>% mutate(cluster = as.character(cluster))
}

apply_annotations <- function(obj, auto_annotations, marker_table, config = read_config()) {
  obj <- sanitize_seurat_metadata(obj)
  auto_annotations <- sanitize_table(auto_annotations)
  marker_table <- sanitize_table(marker_table)
  manual <- load_manual_annotations(config$annotation$manual_annotation_file)
  ann <- auto_annotations %>%
    rename(cell_type = predicted_cell_type) %>%
    mutate(cluster = as.character(cluster))

  if (nrow(manual) > 0) {
    ann <- ann %>%
      dplyr::select(cluster, auto_cell_type = cell_type, best_marker_score, n_cells) %>%
      left_join(manual %>% dplyr::select(cluster, manual_cell_type = cell_type, manual_compartment = compartment, notes), by = "cluster") %>%
      mutate(
        cell_type = ifelse(!is.na(manual_cell_type) & manual_cell_type != "", manual_cell_type, auto_cell_type),
        annotation_source = ifelse(!is.na(manual_cell_type) & manual_cell_type != "", "manual", "marker_score")
      )
  } else {
    ann <- ann %>% mutate(annotation_source = "marker_score", notes = NA_character_)
  }

  comp_lookup <- marker_table %>% dplyr::distinct(cell_type, compartment)
  ann <- ann %>% dplyr::left_join(comp_lookup, by = "cell_type") %>% sanitize_table()
  if (!"manual_compartment" %in% names(ann)) ann$manual_compartment <- NA_character_
  ann$compartment <- ifelse(!is.na(ann$manual_compartment) & ann$manual_compartment != "", ann$manual_compartment, ann$compartment)
  ann$compartment[is.na(ann$compartment) | ann$compartment == ""] <- ann$cell_type[is.na(ann$compartment) | ann$compartment == ""]

  ann$cluster <- as.character(unlist(ann$cluster))
  ann$cell_type <- as.character(unlist(ann$cell_type))
  ann$compartment <- as.character(unlist(ann$compartment))
  cluster_to_type <- stats::setNames(ann$cell_type, ann$cluster)
  cluster_to_comp <- stats::setNames(ann$compartment, ann$cluster)
  obj$cell_type <- unname(cluster_to_type[as.character(obj$seurat_clusters)])
  obj$compartment <- unname(cluster_to_comp[as.character(obj$seurat_clusters)])
  obj$cell_type[is.na(obj$cell_type)] <- paste0("Cluster_", obj$seurat_clusters[is.na(obj$cell_type)])
  obj$compartment[is.na(obj$compartment)] <- obj$cell_type[is.na(obj$compartment)]
  obj@misc$cell_type_annotation <- ann
  obj
}

score_required_cell_states <- function(obj, config = read_config()) {
  obj <- sanitize_seurat_metadata(obj)
  gene_sets <- yaml::read_yaml("resources/compartment_gene_sets.yml")
  state_sets <- gene_sets$fibrosis_mechanism_sets
  state_sets <- purrr::map(state_sets, clean_gene_symbols)
  obj <- add_marker_module_scores(obj, state_sets, prefix = "STATE_")

  obj <- sanitize_seurat_metadata(obj)
  meta <- obj@meta.data
  state_cols <- grep("^STATE_", colnames(meta), value = TRUE)
  if (length(state_cols) > 0) {
    state_mat <- as.data.frame(meta[, state_cols, drop = FALSE], stringsAsFactors = FALSE)
    for (cc in state_cols) state_mat[[cc]] <- suppressWarnings(as.numeric(state_mat[[cc]]))
    max_state <- apply(as.matrix(state_mat), 1, function(x) {
      if (all(is.na(x))) return(NA_character_)
      colnames(state_mat)[which.max(x)]
    })
    obj$fibrosis_program_top <- stringr::str_remove(max_state, "^STATE_")
  }

  if (!"fibrosis_program_top" %in% colnames(obj@meta.data)) obj$fibrosis_program_top <- NA_character_
  obj$cell_state <- as.character(obj$cell_type)
  obj$cell_state[obj$compartment == "Macrophage_monocyte" & obj$fibrosis_program_top == "scar_macrophage"] <- "Scar_associated_Macrophage"
  obj$cell_state[obj$compartment == "HSC_mesenchymal" & obj$fibrosis_program_top == "activated_stellate"] <- "Activated_HSC_Myofibroblast"
  obj$cell_state[obj$compartment == "Endothelial" & obj$fibrosis_program_top == "endothelial_remodeling"] <- "ACKR1_PLVAP_Endothelial"
  obj
}

plot_annotation_figures <- function(obj, marker_table = load_marker_table(), config = read_config()) {
  obj <- sanitize_seurat_metadata(obj)
  marker_table <- sanitize_table(marker_table)
  fig_dir <- config$paths$figure_dir

  p_cell <- Seurat::DimPlot(obj, group.by = "cell_type", reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Major liver cell types")
  save_plot_safely(p_cell, file.path(fig_dir, "umap_cell_type.pdf"), width = 10, height = 7)

  p_comp <- Seurat::DimPlot(obj, group.by = "compartment", reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Compartments")
  save_plot_safely(p_comp, file.path(fig_dir, "umap_compartment.pdf"), width = 10, height = 7)

  major_markers <- marker_table %>%
    filter(marker_role == "canonical") %>%
    group_by(cell_type) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    pull(gene) %>% unique() %>% intersect(rownames(obj))

  if (length(major_markers) > 0) {
    p_dot <- Seurat::DotPlot(obj, features = major_markers, group.by = "cell_type") +
      Seurat::RotatedAxis() + ggtitle("Canonical cell-type markers")
    save_plot_safely(p_dot, file.path(fig_dir, "dotplot_major_celltype_markers.pdf"), width = 14, height = 8)
  }

  req_sets <- yaml::read_yaml("resources/compartment_gene_sets.yml")$required_compartments
  req_markers <- unlist(purrr::map(req_sets, "markers"), use.names = FALSE) %>% unique() %>% clean_gene_symbols() %>% intersect(rownames(obj))
  if (length(req_markers) > 0) {
    p_req <- Seurat::DotPlot(obj, features = req_markers, group.by = "compartment") +
      Seurat::RotatedAxis() + ggtitle("Required disease-relevant compartment marker validation")
    save_plot_safely(p_req, file.path(fig_dir, "dotplot_required_compartments.pdf"), width = 14, height = 7)
  }

  comp <- sanitize_table(obj@meta.data) %>%
    dplyr::count(donor_id, disease_status, cell_type, name = "n_cells") %>%
    group_by(donor_id) %>%
    mutate(frac_cells = n_cells / sum(n_cells)) %>%
    ungroup()
  write_csv_safely(comp, file.path(config$paths$table_dir, "cell_type_counts_by_sample.csv"))

  p_bar <- comp %>%
    ggplot(aes(x = donor_id, y = frac_cells, fill = cell_type)) +
    geom_col(width = 0.9) +
    facet_grid(~ disease_status, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 10) +
    labs(x = "Donor", y = "Fraction of retained cells", fill = "Cell type", title = "Cell-type composition by donor")
  save_plot_safely(p_bar, file.path(fig_dir, "cell_type_composition_by_donor.pdf"), width = 12, height = 7)

  state_cols <- grep("^STATE_", colnames(obj@meta.data), value = TRUE)
  if (length(state_cols) > 0) {
    state_df <- sanitize_table(obj@meta.data) %>%
      tibble::rownames_to_column("cell") %>%
      filter(compartment %in% c("HSC_mesenchymal", "Macrophage_monocyte", "Endothelial")) %>%
      dplyr::select(cell, disease_status, donor_id, compartment, all_of(state_cols)) %>%
      pivot_longer(cols = all_of(state_cols), names_to = "program", values_to = "score") %>%
      mutate(program = stringr::str_remove(program, "^STATE_"))

    state_df <- state_df %>%
      dplyr::mutate(
        disease_status = factor(as.character(disease_status), levels = c("healthy", "cirrhotic")),
        program = stringr::str_replace_all(program, "_", " ")
      )

    state_method <- tibble::tibble(
      item = c("Score method", "Expression layer", "Program gene sets", "Interpretation"),
      description = c(
        "Seurat AddModuleScore was applied to each fibrosis mechanism gene set. For each cell, the displayed score is the average expression of the program genes minus the average expression of matched control genes binned by expression level.",
        "Scores are computed from the normalized RNA assay after QC and normalization; they are used for cell-state visualization and heuristic labeling, not as donor-aware differential-expression statistics.",
        "Gene sets are defined in resources/compartment_gene_sets.yml under fibrosis_mechanism_sets.",
        "Higher values indicate stronger relative activity of the corresponding fibrosis mechanism program within that cell."
      )
    )
    write_csv_safely(state_method, file.path(config$paths$table_dir, "cell_state_score_method.csv"))

    p_state <- state_df %>%
      ggplot(aes(x = disease_status, y = score, fill = disease_status)) +
      geom_boxplot(outlier.size = 0.1, alpha = 0.85, width = 0.72) +
      facet_grid(program ~ compartment, scales = "free_y") +
      theme_bw(base_size = 9) +
      theme(legend.position = "top") +
      labs(
        x = NULL,
        y = "Module score",
        fill = "Disease status",
        title = "Fibrosis program scores in required compartments",
        subtitle = "Scores use Seurat AddModuleScore: program-gene expression minus matched control-gene expression",
        caption = "Gene sets are defined in resources/compartment_gene_sets.yml. Scores support cell-state interpretation; donor-aware DE is performed separately using pseudobulk counts."
      )
    save_plot_safely(p_state, file.path(fig_dir, "cell_state_scores_required_compartments.pdf"), width = 13, height = 10)
  }

  invisible(TRUE)
}

annotate_primary <- function(qc_result, config = read_config()) {
  message_header("Annotating primary dataset")
  obj <- sanitize_seurat_metadata(qc_result$seurat)
  marker_table <- sanitize_table(load_marker_table())
  canonical_sets <- marker_list_by_cell_type(marker_table, role = "canonical")
  obj <- add_marker_module_scores(obj, canonical_sets, prefix = "MS_")

  scores <- cluster_annotation_scores(obj, canonical_sets, cluster_var = "seurat_clusters", prefix = "MS_")
  write_csv_safely(scores, file.path(config$paths$table_dir, "cluster_annotation_scores_long.csv"))

  auto_ann <- initial_annotation_from_scores(scores, cluster_var = "seurat_clusters")
  write_csv_safely(auto_ann, file.path(config$paths$table_dir, "cluster_annotation_scores.csv"))

  obj <- apply_annotations(obj, auto_ann, marker_table, config)
  obj <- score_required_cell_states(obj, config)
  plot_annotation_figures(obj, marker_table, config)

  save_rds_safely(obj, file.path(config$paths$rds_dir, "gse136103_annotated.rds"))
  obj
}
