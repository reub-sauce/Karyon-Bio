# Donor-aware pseudobulk differential expression and cell proportion testing.
# v11: deliberately avoids fragile dplyr/tibble operations during pseudobulk
# assembly because Seurat/GEO metadata can contain list-like columns on some
# systems. The DE construction below first coerces all relevant metadata to
# plain base R vectors, then uses base split/aggregate operations.

atomic_chr <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.data.frame(x) || is.matrix(x) || is.array(x)) x <- as.vector(x)
  if (is.list(x)) {
    return(vapply(x, function(v) {
      if (length(v) == 0) return(NA_character_)
      paste(as.character(unlist(v, use.names = FALSE)), collapse = ";")
    }, character(1)))
  }
  as.character(x)
}

atomic_num <- function(x) {
  suppressWarnings(as.numeric(atomic_chr(x)))
}

metadata_for_de <- function(obj) {
  md <- sanitize_table(obj@meta.data)
  md <- as.data.frame(md, stringsAsFactors = FALSE, optional = TRUE)
  md$cell <- rownames(obj@meta.data)
  if (is.null(md$cell) || length(md$cell) != ncol(obj)) md$cell <- colnames(obj)

  required <- c(
    "cell", "cell_type", "cell_state", "donor_id", "sample_id", "disease_status",
    "fibrosis_stage_harmonized", "fraction"
  )
  for (nm in required) {
    if (!nm %in% names(md)) md[[nm]] <- NA_character_
    md[[nm]] <- atomic_chr(md[[nm]])
  }
  md$nCount_RNA <- if ("nCount_RNA" %in% names(md)) atomic_num(md$nCount_RNA) else NA_real_
  md$nFeature_RNA <- if ("nFeature_RNA" %in% names(md)) atomic_num(md$nFeature_RNA) else NA_real_

  md <- md[md$cell %in% colnames(obj), , drop = FALSE]
  md <- md[match(intersect(colnames(obj), md$cell), md$cell), , drop = FALSE]
  rownames(md) <- md$cell
  md
}

write_pseudobulk_debug <- function(md, config = read_config(), label = "pseudobulk_metadata_debug") {
  out <- md[, intersect(c("cell", "cell_type", "cell_state", "donor_id", "sample_id", "disease_status", "fraction"), names(md)), drop = FALSE]
  write_csv_safely(head(out, 2000), file.path(config$paths$table_dir, paste0(label, ".csv")))
  invisible(out)
}

make_pseudobulk <- function(obj, group_var = "cell_type", sample_var = "donor_id", min_cells = 20, config = read_config()) {
  meta <- metadata_for_de(obj)
  if (!group_var %in% names(meta)) stop("Missing group_var in metadata: ", group_var)
  if (!sample_var %in% names(meta)) stop("Missing sample_var in metadata: ", sample_var)

  meta$group_value <- atomic_chr(meta[[group_var]])
  meta$sample_value <- atomic_chr(meta[[sample_var]])
  meta <- meta[!is.na(meta$group_value) & meta$group_value != "" & !is.na(meta$sample_value) & meta$sample_value != "", , drop = FALSE]
  meta$pb_id <- paste(meta$group_value, meta$sample_value, sep = "__")

  if (nrow(meta) == 0) stop("No cells remain after filtering metadata for pseudobulk.")
  tab <- table(meta$pb_id)
  keep_pb <- names(tab)[tab >= min_cells]
  meta <- meta[meta$pb_id %in% keep_pb, , drop = FALSE]
  if (nrow(meta) == 0) stop("No pseudobulk groups have at least ", min_cells, " cells.")

  write_pseudobulk_debug(meta, config, paste0("pseudobulk_metadata_debug_", group_var))

  counts <- get_assay_data_compat(obj, assay = "RNA", layer = "counts")
  common_cells <- intersect(colnames(counts), meta$cell)
  meta <- meta[match(common_cells, meta$cell), , drop = FALSE]
  counts <- counts[, common_cells, drop = FALSE]

  pb_factor <- factor(meta$pb_id, levels = unique(meta$pb_id))
  design <- Matrix::sparse.model.matrix(~ 0 + pb_factor)
  colnames(design) <- levels(pb_factor)
  rownames(design) <- meta$cell
  pb_counts <- counts %*% design

  split_idx <- split(seq_len(nrow(meta)), meta$pb_id)
  pb_meta <- do.call(rbind, lapply(names(split_idx), function(id) {
    ii <- split_idx[[id]]
    data.frame(
      pb_id = id,
      group = meta$group_value[ii][1],
      donor_id = meta$donor_id[ii][1],
      sample_ids = paste(sort(unique(meta$sample_id[ii])), collapse = ";"),
      disease_status = meta$disease_status[ii][1],
      fibrosis_stage_harmonized = meta$fibrosis_stage_harmonized[ii][1],
      fraction = paste(sort(unique(meta$fraction[ii])), collapse = ";"),
      n_cells = length(ii),
      stringsAsFactors = FALSE
    )
  }))
  pb_meta <- sanitize_table(pb_meta)

  pb_counts <- pb_counts[, pb_meta$pb_id, drop = FALSE]
  list(counts = pb_counts, metadata = pb_meta, group_var = group_var)
}

run_edgeR_for_group <- function(pb, group_name, config = read_config()) {
  meta <- as.data.frame(sanitize_table(pb$metadata), stringsAsFactors = FALSE)
  meta <- meta[meta$group == group_name, , drop = FALSE]
  if (nrow(meta) == 0) return(tibble::tibble())

  meta$pb_id <- atomic_chr(meta$pb_id)
  meta$disease_status <- atomic_chr(meta$disease_status)
  meta$donor_id <- atomic_chr(meta$donor_id)
  meta$fraction <- atomic_chr(meta$fraction)
  meta$n_cells <- atomic_num(meta$n_cells)

  keep <- meta$disease_status %in% c(config$pseudobulk$reference_level, config$pseudobulk$disease_contrast_level)
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) == 0) return(tibble::tibble())

  meta$disease_status <- factor(meta$disease_status, levels = c(config$pseudobulk$reference_level, config$pseudobulk$disease_contrast_level))
  donor_tab <- unique(meta[, c("donor_id", "disease_status"), drop = FALSE])
  donors_per_group <- table(donor_tab$disease_status)
  if (length(donors_per_group) < 2 || any(donors_per_group < config$pseudobulk$min_donors_per_group)) {
    message("Skipping ", group_name, ": not enough donors in both disease groups.")
    return(tibble::tibble())
  }

  counts <- as.matrix(pb$counts[, meta$pb_id, drop = FALSE])
  storage.mode(counts) <- "integer"
  samples_df <- as.data.frame(meta, stringsAsFactors = FALSE)
  rownames(samples_df) <- samples_df$pb_id

  y <- edgeR::DGEList(counts = counts, samples = samples_df)
  keep_genes <- edgeR::filterByExpr(y, group = meta$disease_status)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(y) < 10) return(tibble::tibble())
  y <- edgeR::calcNormFactors(y)

  design <- stats::model.matrix(~ disease_status, data = meta)
  if (length(unique(meta$fraction)) > 1) {
    design_try <- stats::model.matrix(~ fraction + disease_status, data = meta)
    if (qr(design_try)$rank == ncol(design_try)) design <- design_try
  }

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  coef_name <- grep(paste0("disease_status", config$pseudobulk$disease_contrast_level), colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    message("Skipping ", group_name, ": could not identify disease coefficient.")
    return(tibble::tibble())
  }

  qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab <- tibble::rownames_to_column(as.data.frame(tab), "gene")
  tab$analysis_group <- group_name
  tab$group_var <- pb$group_var
  tab$contrast <- paste(config$pseudobulk$disease_contrast_level, "vs", config$pseudobulk$reference_level)
  tab$n_pseudobulk_samples <- nrow(meta)
  tab$n_donors_healthy <- sum(meta$disease_status == config$pseudobulk$reference_level, na.rm = TRUE)
  tab$n_donors_cirrhotic <- sum(meta$disease_status == config$pseudobulk$disease_contrast_level, na.rm = TRUE)
  tab$mean_cells_per_pseudobulk <- mean(meta$n_cells, na.rm = TRUE)
  tab$direction <- ifelse(tab$logFC > 0, "up_in_cirrhotic", ifelse(tab$logFC < 0, "down_in_cirrhotic", "flat"))
  sanitize_table(tab)
}

run_pseudobulk_de <- function(obj, group_var = "cell_type", config = read_config()) {
  message_header("Running pseudobulk DE by ", group_var)
  pb <- make_pseudobulk(
    obj,
    group_var = group_var,
    sample_var = config$pseudobulk$sample_var,
    min_cells = config$pseudobulk$min_cells_per_pseudobulk,
    config = config
  )
  groups <- sort(unique(atomic_chr(pb$metadata$group)))
  pieces <- lapply(groups, function(g) {
    tryCatch(
      run_edgeR_for_group(pb, g, config),
      error = function(e) {
        message("Skipping ", g, " because edgeR failed: ", conditionMessage(e))
        tibble::tibble()
      }
    )
  })
  sanitize_table(dplyr::bind_rows(pieces))
}

run_cell_proportion_tests <- function(obj, config = read_config()) {
  message_header("Running donor-level cell proportion tests")
  meta <- metadata_for_de(obj)
  meta <- meta[!is.na(meta$donor_id) & !is.na(meta$disease_status) & !is.na(meta$cell_type), , drop = FALSE]
  if (nrow(meta) == 0) return(tibble::tibble())

  totals <- as.data.frame(table(meta$donor_id, meta$disease_status), stringsAsFactors = FALSE)
  names(totals) <- c("donor_id", "disease_status", "n_total")
  totals <- totals[totals$n_total > 0, , drop = FALSE]
  cell_counts <- as.data.frame(table(meta$donor_id, meta$disease_status, meta$cell_type), stringsAsFactors = FALSE)
  names(cell_counts) <- c("donor_id", "disease_status", "cell_type", "n_cells")
  cell_counts <- merge(cell_counts, totals, by = c("donor_id", "disease_status"), all.x = TRUE)
  cell_counts$n_other <- cell_counts$n_total - cell_counts$n_cells

  results <- lapply(split(cell_counts, cell_counts$cell_type), function(df) {
    df <- df[df$disease_status %in% c(config$pseudobulk$reference_level, config$pseudobulk$disease_contrast_level), , drop = FALSE]
    if (length(unique(df$disease_status)) < 2 || nrow(df) < 4) return(NULL)
    df$disease_status <- factor(df$disease_status, levels = c(config$pseudobulk$reference_level, config$pseudobulk$disease_contrast_level))
    fit <- try(stats::glm(cbind(n_cells, n_other) ~ disease_status, data = df, family = quasibinomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    sm <- summary(fit)$coefficients
    coef_name <- grep(paste0("disease_status", config$pseudobulk$disease_contrast_level), rownames(sm), value = TRUE)
    if (length(coef_name) != 1) return(NULL)
    data.frame(
      cell_type = unique(df$cell_type)[1],
      estimate_log_odds = sm[coef_name, "Estimate"],
      std_error = sm[coef_name, "Std. Error"],
      statistic = sm[coef_name, "t value"],
      p_value = sm[coef_name, "Pr(>|t|)"],
      mean_fraction_healthy = mean(df$n_cells[df$disease_status == config$pseudobulk$reference_level] / df$n_total[df$disease_status == config$pseudobulk$reference_level], na.rm = TRUE),
      mean_fraction_cirrhotic = mean(df$n_cells[df$disease_status == config$pseudobulk$disease_contrast_level] / df$n_total[df$disease_status == config$pseudobulk$disease_contrast_level], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  res <- sanitize_table(dplyr::bind_rows(results))
  if (nrow(res) > 0) res$fdr <- p.adjust(res$p_value, method = "BH")
  write_csv_safely(res, file.path(config$paths$table_dir, "cell_proportion_tests.csv"))
  res
}

required_compartment_filter <- function(de_tbl) {
  if (is.null(de_tbl) || nrow(de_tbl) == 0) return(tibble::tibble())
  de_tbl <- sanitize_table(de_tbl)
  de_tbl %>% dplyr::filter(
    analysis_group %in% c("HSC_mesenchymal", "Macrophage_monocyte", "Endothelial") |
      stringr::str_detect(analysis_group, "HSC|Mesenchymal|Macrophage|Monocyte|Endothelial|Scar|Myofibroblast|ACKR1|PLVAP")
  )
}

run_all_pseudobulk_de <- function(obj, config = read_config()) {
  obj <- sanitize_seurat_metadata(obj)
  de_cell_type <- run_pseudobulk_de(obj, group_var = "cell_type", config = config)
  de_cell_state <- run_pseudobulk_de(obj, group_var = "cell_state", config = config)
  all_de <- sanitize_table(dplyr::bind_rows(de_cell_type, de_cell_state))

  write_csv_safely(all_de, file.path(config$paths$table_dir, "pseudobulk_de_all_cell_types.csv"))
  req <- required_compartment_filter(all_de) %>% sanitize_table()
  write_csv_safely(req, file.path(config$paths$table_dir, "pseudobulk_de_required_compartments.csv"))
  prop <- run_cell_proportion_tests(obj, config) %>% sanitize_table()

  list(all_de = all_de, required_de = req, proportion_tests = prop)
}
