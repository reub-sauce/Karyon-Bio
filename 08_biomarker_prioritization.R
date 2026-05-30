# Biomarker / target prioritization score.
# v14: robust to missing validation counts/results, missing genes, Seurat v5 layer
# behavior, duplicated validation/tag rows, cell-name mismatches, and empty expression
# specificity tables. This prevents rowMeans()/rename() errors when optional GSE207310
# validation is unavailable or when one grouping variable is absent after filtering.

load_translational_tags <- function(path = "resources/translational_gene_tags.csv") {
  tags <- sanitize_table(readr::read_csv(path, show_col_types = FALSE)) %>%
    dplyr::mutate(gene = clean_gene_symbols(gene))
  # Enforce one row per gene to avoid many-to-many joins during ranking.
  tags %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      translational_class = dplyr::first(stats::na.omit(translational_class)) %||% NA_character_,
      secreted_or_surface = paste(sort(unique(stats::na.omit(secreted_or_surface))), collapse = ";"),
      interpretation_note = paste(sort(unique(stats::na.omit(interpretation_note))), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      secreted_or_surface = dplyr::na_if(secreted_or_surface, ""),
      interpretation_note = dplyr::na_if(interpretation_note, "")
    )
}

get_expression_matrix_for_scoring <- function(obj, genes, assay = "RNA") {
  genes <- intersect(clean_gene_symbols(genes), rownames(obj))
  if (length(genes) == 0) return(NULL)

  # Prefer normalized data when available; fall back to counts if the data layer is
  # empty or absent. The layer warnings in Seurat v5 are harmless, but an empty data
  # layer would create downstream zero-column/zero-row issues.
  mat <- tryCatch(get_assay_data_compat(obj, assay = assay, layer = "data"), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    mat <- get_assay_data_compat(obj, assay = assay, layer = "counts")
  }
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0 || ncol(mat) == 0) return(NULL)
  mat[genes, , drop = FALSE]
}

compute_expression_specificity <- function(obj, genes, group_var = "cell_type") {
  data <- get_expression_matrix_for_scoring(obj, genes, assay = "RNA")
  if (is.null(data)) return(tibble::tibble())

  genes <- rownames(data)
  meta <- sanitize_table(obj@meta.data) %>%
    tibble::rownames_to_column("cell")

  if (!group_var %in% names(meta)) return(tibble::tibble())
  if (!"disease_status" %in% names(meta)) meta$disease_status <- NA_character_

  meta <- meta %>%
    dplyr::transmute(
      cell = as.character(cell),
      group = as.character(.data[[group_var]]),
      disease_status = as.character(disease_status)
    ) %>%
    dplyr::filter(!is.na(group), group != "", cell %in% colnames(data))

  groups <- sort(unique(meta$group))
  if (length(groups) == 0) return(tibble::tibble())

  safe_row_means <- function(m) {
    if (is.null(m) || ncol(m) == 0) return(rep(NA_real_, length(genes)))
    as.numeric(Matrix::rowMeans(m))
  }

  avg <- purrr::map_dfr(groups, function(g) {
    cells <- intersect(meta$cell[meta$group == g], colnames(data))
    if (length(cells) == 0) return(tibble::tibble())
    sub <- data[, cells, drop = FALSE]
    tibble::tibble(
      gene = genes,
      group = g,
      avg_expr = safe_row_means(sub),
      pct_expr = safe_row_means(sub > 0)
    )
  })

  if (nrow(avg) == 0) return(tibble::tibble())

  tau <- avg %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      tau_specificity = {
        x <- as.numeric(avg_expr)
        mx <- suppressWarnings(max(x, na.rm = TRUE))
        if (!is.finite(mx) || mx <= 0 || length(x) <= 1) 0 else sum(1 - x / mx, na.rm = TRUE) / (length(x) - 1)
      },
      .groups = "drop"
    )

  disease_pct <- purrr::map_dfr(groups, function(g) {
    cells <- intersect(meta$cell[meta$group == g & meta$disease_status == "cirrhotic"], colnames(data))
    if (length(cells) == 0) {
      return(tibble::tibble(gene = genes, group = g, pct_expr_cirrhotic = NA_real_))
    }
    tibble::tibble(
      gene = genes,
      group = g,
      pct_expr_cirrhotic = safe_row_means(data[, cells, drop = FALSE] > 0)
    )
  })

  avg %>%
    dplyr::left_join(tau, by = "gene") %>%
    dplyr::left_join(disease_pct, by = c("gene", "group")) %>%
    sanitize_table()
}


empty_expression_specificity <- function() {
  tibble::tibble(
    gene = character(),
    analysis_group = character(),
    avg_expr = numeric(),
    pct_expr = numeric(),
    tau_specificity = numeric(),
    pct_expr_cirrhotic = numeric()
  )
}

expression_specificity_for_groupvar <- function(obj, genes, group_var) {
  out <- tryCatch(
    compute_expression_specificity(obj, genes, group_var = group_var),
    error = function(e) {
      warning("Expression specificity failed for group_var=", group_var, ": ", conditionMessage(e))
      tibble::tibble()
    }
  )
  out <- sanitize_table(out)
  if (is.null(out) || nrow(out) == 0 || !"group" %in% names(out)) {
    return(empty_expression_specificity())
  }
  if (!"gene" %in% names(out)) out$gene <- character()
  if (!"avg_expr" %in% names(out)) out$avg_expr <- NA_real_
  if (!"pct_expr" %in% names(out)) out$pct_expr <- NA_real_
  if (!"tau_specificity" %in% names(out)) out$tau_specificity <- 0
  if (!"pct_expr_cirrhotic" %in% names(out)) out$pct_expr_cirrhotic <- NA_real_
  out %>%
    dplyr::mutate(
      gene = clean_gene_symbols(gene),
      analysis_group = as.character(group)
    ) %>%
    dplyr::select(gene, analysis_group, avg_expr, pct_expr, tau_specificity, pct_expr_cirrhotic) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(analysis_group), analysis_group != "") %>%
    sanitize_table()
}

prepare_candidate_table <- function(de_result, config = read_config()) {
  de <- sanitize_table(de_result$required_de)
  if (nrow(de) == 0) return(tibble::tibble())
  de %>%
    dplyr::mutate(gene = clean_gene_symbols(gene)) %>%
    dplyr::filter(!grepl("^MT-|^RPL|^RPS|^HB[ABDEGQMZ]", gene)) %>%
    dplyr::filter(!is.na(FDR), !is.na(logFC)) %>%
    dplyr::filter(FDR <= config$pseudobulk$fdr_cutoff, abs(logFC) >= config$pseudobulk$abs_logfc_cutoff) %>%
    dplyr::mutate(candidate_id = paste(analysis_group, gene, sep = "__")) %>%
    dplyr::distinct(candidate_id, .keep_all = TRUE)
}

mechanism_score_table <- function(pathway_result) {
  mech <- sanitize_table(pathway_result$curated_mechanism_hits)
  if (is.null(mech) || nrow(mech) == 0) {
    return(tibble::tibble(gene = character(), analysis_group = character(), mechanism_score_raw = numeric(), mechanism_labels = character()))
  }
  mech %>%
    dplyr::mutate(gene = clean_gene_symbols(gene)) %>%
    dplyr::group_by(gene, analysis_group) %>%
    dplyr::summarise(
      mechanism_score_raw = dplyr::n_distinct(mechanism),
      mechanism_labels = paste(sort(unique(mechanism)), collapse = ";"),
      .groups = "drop"
    )
}

validation_score_table <- function(validation_result) {
  val <- tryCatch(sanitize_table(validation_result$validation_results), error = function(e) tibble::tibble())
  if (is.null(val) || nrow(val) == 0 || !"gene" %in% names(val)) {
    return(tibble::tibble(
      gene = character(),
      validation_score_raw = numeric(),
      validation_log2fc = numeric(),
      validation_padj = numeric(),
      validation_model = character(),
      validation_contrast = character(),
      validation_direction = character()
    ))
  }

  if (!"padj" %in% names(val)) val$padj <- NA_real_
  if (!"log2FoldChange" %in% names(val)) val$log2FoldChange <- NA_real_
  if (!"validation_model" %in% names(val)) val$validation_model <- NA_character_
  if (!"validation_contrast" %in% names(val)) val$validation_contrast <- NA_character_
  if (!"validation_direction" %in% names(val)) val$validation_direction <- NA_character_

  val %>%
    dplyr::mutate(
      gene = clean_gene_symbols(gene),
      validation_padj = suppressWarnings(as.numeric(padj)),
      validation_log2fc = suppressWarnings(as.numeric(log2FoldChange)),
      validation_score_raw = rescale_01(abs(validation_log2fc)) * 0.5 + rescale_01(safe_neg_log10(validation_padj)) * 0.5
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::arrange(validation_padj, dplyr::desc(abs(validation_log2fc)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, validation_score_raw, validation_log2fc, validation_padj, validation_model, validation_contrast, validation_direction)
}

assign_translational_relevance <- function(tbl) {
  tbl %>%
    dplyr::mutate(
      likely_use_case = dplyr::case_when(
        translational_class %in% c("diagnostic_marker") & grepl("secreted|surface|plasma|ECM|ecm", secreted_or_surface, ignore.case = TRUE) ~ "diagnostic / tissue or blood biomarker",
        translational_class %in% c("therapeutic_target") ~ "therapeutic target candidate",
        translational_class %in% c("state_marker") ~ "cell-state validation marker",
        !is.na(validation_padj) & validation_padj <= 0.1 ~ "external validation candidate",
        TRUE ~ "future validation candidate"
      ),
      translational_note = dplyr::case_when(
        !is.na(interpretation_note) & interpretation_note != "" ~ interpretation_note,
        direction == "up_in_cirrhotic" ~ "Upregulated in cirrhosis in a disease-relevant compartment",
        TRUE ~ "Disease-associated candidate requiring orthogonal validation"
      )
    )
}

prioritize_biomarkers <- function(obj, de_result, pathway_result, validation_result, config = read_config()) {
  message_header("Prioritizing biomarkers and therapeutic targets")
  candidates <- sanitize_table(prepare_candidate_table(de_result, config))
  if (nrow(candidates) == 0) {
    warning("No candidates passed DE thresholds; relaxing to top FDR-ranked required-compartment genes.")
    candidates <- sanitize_table(de_result$required_de) %>%
      dplyr::mutate(gene = clean_gene_symbols(gene)) %>%
      dplyr::filter(!grepl("^MT-|^RPL|^RPS|^HB[ABDEGQMZ]", gene)) %>%
      dplyr::arrange(FDR) %>%
      dplyr::group_by(analysis_group) %>%
      dplyr::slice_head(n = 50) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(candidate_id = paste(analysis_group, gene, sep = "__")) %>%
      dplyr::distinct(candidate_id, .keep_all = TRUE)
  }

  if (nrow(candidates) == 0) {
    empty <- tibble::tibble()
    write_csv_safely(empty, file.path(config$paths$table_dir, "ranked_biomarkers_full.csv"))
    write_csv_safely(empty, file.path(config$paths$table_dir, "ranked_biomarkers_top20.csv"))
    return(list(full = empty, top20 = empty))
  }

  genes <- unique(candidates$gene)
  expr <- dplyr::bind_rows(
    expression_specificity_for_groupvar(obj, genes, "cell_type"),
    expression_specificity_for_groupvar(obj, genes, "cell_state")
  ) %>%
    dplyr::distinct(gene, analysis_group, .keep_all = TRUE) %>%
    sanitize_table()

  if (is.null(expr) || nrow(expr) == 0) {
    expr <- candidates %>%
      dplyr::distinct(gene, analysis_group) %>%
      dplyr::mutate(
        avg_expr = NA_real_,
        pct_expr = NA_real_,
        tau_specificity = 0,
        pct_expr_cirrhotic = NA_real_
      ) %>%
      sanitize_table()
  }

  mech <- mechanism_score_table(pathway_result)
  val <- validation_score_table(validation_result)
  tags <- load_translational_tags()
  weights <- config$biomarker_score$weights

  expr <- sanitize_table(expr) %>% dplyr::distinct(gene, analysis_group, .keep_all = TRUE)
  mech <- sanitize_table(mech) %>% dplyr::distinct(gene, analysis_group, .keep_all = TRUE)
  val <- sanitize_table(val) %>% dplyr::distinct(gene, .keep_all = TRUE)
  tags <- sanitize_table(tags) %>% dplyr::distinct(gene, .keep_all = TRUE)

  scored <- sanitize_table(candidates) %>%
    dplyr::distinct(candidate_id, .keep_all = TRUE) %>%
    dplyr::left_join(expr, by = c("gene", "analysis_group")) %>%
    dplyr::left_join(mech, by = c("gene", "analysis_group")) %>%
    dplyr::left_join(val, by = "gene") %>%
    dplyr::left_join(tags, by = "gene") %>%
    dplyr::mutate(
      logFC = suppressWarnings(as.numeric(logFC)),
      FDR = suppressWarnings(as.numeric(FDR)),
      effect_score = rescale_01(abs(logFC)),
      significance_score = rescale_01(safe_neg_log10(FDR)),
      specificity_score = dplyr::coalesce(suppressWarnings(as.numeric(tau_specificity)), 0),
      prevalence_score = dplyr::coalesce(suppressWarnings(as.numeric(pct_expr_cirrhotic)), suppressWarnings(as.numeric(pct_expr)), 0),
      mechanism_score = rescale_01(dplyr::coalesce(suppressWarnings(as.numeric(mechanism_score_raw)), 0)),
      translational_score = dplyr::case_when(
        !is.na(translational_class) & translational_class == "therapeutic_target" ~ 1.0,
        !is.na(translational_class) & translational_class == "diagnostic_marker" ~ 0.9,
        !is.na(translational_class) & translational_class == "state_marker" ~ 0.5,
        TRUE ~ 0.2
      ),
      validation_direction_consistent = dplyr::case_when(
        is.na(validation_log2fc) ~ NA,
        direction == "up_in_cirrhotic" & validation_log2fc > 0 ~ TRUE,
        direction == "down_in_cirrhotic" & validation_log2fc < 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      validation_score = dplyr::case_when(
        validation_direction_consistent %in% TRUE ~ dplyr::coalesce(validation_score_raw, 0),
        validation_direction_consistent %in% FALSE ~ 0.1 * dplyr::coalesce(validation_score_raw, 0),
        TRUE ~ 0
      ),
      total_score = weights$effect * effect_score +
        weights$significance * significance_score +
        weights$specificity * specificity_score +
        weights$prevalence * prevalence_score +
        weights$mechanism * mechanism_score +
        weights$translational * translational_score +
        weights$validation * validation_score,
      rank = dplyr::dense_rank(dplyr::desc(total_score))
    ) %>%
    assign_translational_relevance() %>%
    dplyr::arrange(rank, FDR) %>%
    sanitize_table()

  full_path <- file.path(config$paths$table_dir, "ranked_biomarkers_full.csv")
  top_path <- file.path(config$paths$table_dir, "ranked_biomarkers_top20.csv")
  write_csv_safely(scored, full_path)
  top <- sanitize_table(scored) %>% dplyr::slice_head(n = config$biomarker_score$top_n)
  write_csv_safely(top, top_path)

  if (nrow(top) > 0) {
    p <- top %>%
      dplyr::mutate(label = paste(gene, analysis_group, sep = " / ")) %>%
      ggplot2::ggplot(ggplot2::aes(x = total_score, y = reorder(label, total_score))) +
      ggplot2::geom_col() +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(x = "Prioritization score", y = NULL, title = "Top ranked fibrosis biomarkers / targets")
    save_plot_safely(p, file.path(config$paths$figure_dir, "ranked_biomarkers_top20.pdf"), width = 10, height = 7)
  }

  list(full = scored, top20 = top)
}
