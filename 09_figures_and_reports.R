# Report rendering helpers.

render_quarto_if_available <- function(input, output_dir = "results/reports", output_file = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- output_file %||% basename(sub("\\.qmd$", ".html", input))
  output_path <- file.path(output_dir, output_file)

  render_fallback <- function(reason = NULL) {
    if (!is.null(reason)) {
      message("Quarto rendering was not available or failed: ", conditionMessage(reason))
    }
    message("Writing a lightweight fallback HTML report for: ", input)
    write_fallback_html_report(output_path = output_path)
    normalizePath(output_path, mustWork = FALSE)
  }

  if (!(requireNamespace("quarto", quietly = TRUE) && isTRUE(quarto::quarto_available()))) {
    message("Quarto CLI was not detected.")
    return(render_fallback())
  }

  message("Rendering polished HTML report with Quarto: ", input)

  # Some R quarto package versions do not support `output_dir`. Build the
  # render call from available formal arguments and copy/move the result into
  # results/reports afterward when needed.
  qformals <- names(formals(quarto::quarto_render))

  rendered <- tryCatch({
    if ("output_dir" %in% qformals) {
      args <- list(
        input = input,
        output_dir = output_dir,
        output_file = output_file,
        quiet = TRUE
      )
      args <- args[names(args) %in% qformals]
      do.call(quarto::quarto_render, args)
    } else {
      # Older quarto R package: render beside the template, then copy.
      args <- list(
        input = input,
        output_file = output_file,
        quiet = TRUE
      )
      args <- args[names(args) %in% qformals]
      do.call(quarto::quarto_render, args)

      candidate_paths <- unique(c(
        file.path(dirname(input), output_file),
        file.path(getwd(), output_file),
        output_file
      ))
      existing <- candidate_paths[file.exists(candidate_paths)]
      if (length(existing) > 0) {
        file.copy(existing[[1]], output_path, overwrite = TRUE)
      }
    }
    TRUE
  }, error = function(e) e)

  if (inherits(rendered, "error")) {
    return(render_fallback(rendered))
  }

  if (!file.exists(output_path)) {
    # Last-resort search in case quarto changed the output name or location.
    possible <- list.files(
      path = c(dirname(input), getwd(), output_dir),
      pattern = paste0("^", tools::file_path_sans_ext(basename(output_file)), ".*\\.html$"),
      full.names = TRUE,
      recursive = FALSE
    )
    possible <- possible[file.exists(possible)]
    if (length(possible) > 0) {
      file.copy(possible[[1]], output_path, overwrite = TRUE)
    }
  }

  if (!file.exists(output_path)) {
    message("Quarto did not create the expected output file; using fallback HTML.")
    write_fallback_html_report(output_path = output_path)
  }

  normalizePath(output_path, mustWork = FALSE)
}

html_escape_simple <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

read_csv_report <- function(path) {
  if (file.exists(path)) {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    tibble::tibble()
  }
}

first_existing_col <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

sum_first_existing_col <- function(x, candidates) {
  if (is.null(x) || nrow(x) == 0) return(0)
  col <- first_existing_col(x, candidates)
  if (is.na(col)) return(0)
  suppressWarnings(sum(as.numeric(x[[col]]), na.rm = TRUE))
}

html_table_report <- function(x, max_rows = 20) {
  if (is.null(x) || nrow(x) == 0) return("<p><em>No table available.</em></p>")
  x <- head(as.data.frame(x), max_rows)
  x[] <- lapply(x, function(col) html_escape_simple(as.character(col)))
  header <- paste0("<tr>", paste0("<th>", names(x), "</th>", collapse = ""), "</tr>")
  rows <- apply(x, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", header, paste(rows, collapse = "\n"), "</table>")
}

fallback_metric_card <- function(value, label) {
  value <- ifelse(is.na(value) || length(value) == 0, 0, value)
  paste0("<div class='metric-card'><div class='metric-value'>", html_escape_simple(format(value, big.mark = ",")), "</div><div class='metric-label'>", html_escape_simple(label), "</div></div>")
}

fallback_report_asset_dir <- function() {
  dir.create("results/reports/report_assets", recursive = TRUE, showWarnings = FALSE)
  "results/reports/report_assets"
}

fallback_pdf_preview <- function(path) {
  if (!file.exists(path)) return(NULL)
  asset_dir <- fallback_report_asset_dir()
  out <- file.path(asset_dir, paste0(tools::file_path_sans_ext(basename(path)), "_page1.png"))
  if (file.exists(out)) return(out)
  if (requireNamespace("pdftools", quietly = TRUE)) {
    converted <- tryCatch({
      pdftools::pdf_convert(path, format = "png", pages = 1, dpi = 160, filenames = out, verbose = FALSE)
    }, error = function(e) NULL)
    if (!is.null(converted) && file.exists(out)) return(out)
  }
  NULL
}

fallback_figure_embed <- function(path, title = basename(path), height = 520) {
  if (!file.exists(path)) return("")
  ext <- tolower(tools::file_ext(path))
  title <- html_escape_simple(title)

  # Most browsers do not reliably render local PDFs inside <object>/<embed>.
  # When pdftools is installed, convert page 1 to a PNG preview and link to the
  # source PDF. This makes the HTML report reviewer-friendly and self-contained
  # enough for local inspection.
  if (ext == "pdf") {
    preview <- fallback_pdf_preview(path)
    pdf_src <- file.path("..", "figures", basename(path))
    if (!is.null(preview) && file.exists(preview)) {
      preview_src <- file.path("report_assets", basename(preview))
      return(paste0("<figure class='report-figure'><a href='", pdf_src, "'><img src='", preview_src, "' alt='", title, "'></a><figcaption><a href='", pdf_src, "'>", title, "</a></figcaption></figure>"))
    }
    return(paste0("<figure class='report-figure pdf-inline'><iframe src='", pdf_src, "' title='", title, "' style='width:100%;height:", height, "px;border:0;border-radius:10px;background:#f8fafc;'></iframe><figcaption><a href='", pdf_src, "'>", title, "</a></figcaption></figure>"))
  }

  src <- file.path("..", "figures", basename(path))
  if (ext %in% c("png", "jpg", "jpeg", "gif", "svg")) {
    paste0("<figure class='report-figure'><img src='", src, "' alt='", title, "'><figcaption>", title, "</figcaption></figure>")
  } else {
    paste0("<p><a href='", src, "'>", title, "</a></p>")
  }
}

report_numeric_column <- function(x, candidates, default = NA_real_) {
  if (is.null(x) || nrow(x) == 0) return(rep(default, max(0, nrow(x))))
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) return(rep(default, nrow(x)))
  suppressWarnings(as.numeric(x[[hit[[1]]]]))
}

make_validation_figures_from_tables <- function() {
  # Build reviewer-friendly validation figures from CSV outputs. These figures
  # are generated during reporting so they are available even when the Quarto
  # render falls back to lightweight HTML.
  fig_dir <- "results/figures"
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  top <- read_csv_report("results/tables/ranked_biomarkers_top20.csv") %>% sanitize_table()
  de <- read_csv_report("results/tables/pseudobulk_de_required_compartments.csv") %>% sanitize_table()
  val <- read_csv_report("results/tables/gse207310_validation_results.csv") %>% sanitize_table()

  if (nrow(val) == 0 || !("gene" %in% names(val))) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  val$validation_log2fc <- report_numeric_column(val, c("log2FoldChange", "logFC", "validation_log2fc"))
  val$validation_padj <- report_numeric_column(val, c("padj", "adj.P.Val", "pvalue", "P.Value", "validation_padj"))

  val2 <- val %>%
    dplyr::mutate(
      gene = clean_gene_symbols(.data$gene),
      validation_log2fc = .data$validation_log2fc,
      validation_padj = .data$validation_padj
    ) %>%
    dplyr::filter(!is.na(.data$gene), .data$gene != "") %>%
    dplyr::arrange(.data$validation_padj) %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE)

  de2 <- de
  if (nrow(de2) > 0 && all(c("gene", "logFC") %in% names(de2))) {
    de2 <- de2 %>%
      dplyr::mutate(
        gene = clean_gene_symbols(.data$gene),
        discovery_logFC = suppressWarnings(as.numeric(.data$logFC)),
        discovery_FDR = suppressWarnings(as.numeric(.data$FDR))
      ) %>%
      dplyr::filter(!is.na(.data$gene), .data$gene != "") %>%
      dplyr::arrange(.data$discovery_FDR) %>%
      dplyr::distinct(.data$gene, .keep_all = TRUE)
  } else {
    de2 <- tibble::tibble(gene = character(), discovery_logFC = numeric(), discovery_FDR = numeric())
  }

  top_genes <- if (nrow(top) > 0 && "gene" %in% names(top)) unique(clean_gene_symbols(top$gene)) else character()
  discovery_genes <- unique(de2$gene)
  validation_genes <- unique(val2$gene)
  overlap_genes <- intersect(discovery_genes, validation_genes)
  top_overlap_genes <- intersect(top_genes, validation_genes)

  overlap_tbl <- tibble::tibble(
    category = factor(
      c("Discovery genes", "GSE207310 genes", "Discovery-validation overlap", "Top-20 biomarker overlap"),
      levels = c("Discovery genes", "GSE207310 genes", "Discovery-validation overlap", "Top-20 biomarker overlap")
    ),
    n_genes = c(length(discovery_genes), length(validation_genes), length(overlap_genes), length(top_overlap_genes))
  )
  write_csv_safely(overlap_tbl, file.path("results/tables", "gse207310_validation_overlap_summary.csv"))

  p1 <- ggplot2::ggplot(overlap_tbl, ggplot2::aes(x = category, y = n_genes)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Discovery and GSE207310 validation gene overlap",
      x = NULL,
      y = "Number of genes",
      caption = "Discovery genes are from required-compartment pseudobulk DE; validation genes are genes parsed from GSE207310 bulk RNA-seq."
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "validation_discovery_overlap.pdf"), p1, width = 9, height = 4.8)

  scatter <- de2 %>% dplyr::inner_join(val2, by = "gene")
  if (nrow(scatter) > 2) {
    scatter <- scatter %>%
      dplyr::mutate(
        in_top20 = .data$gene %in% top_genes,
        validation_support = dplyr::case_when(
          !is.na(.data$validation_padj) & .data$validation_padj < 0.05 & sign(.data$discovery_logFC) == sign(.data$validation_log2fc) ~ "FDR<0.05, same direction",
          sign(.data$discovery_logFC) == sign(.data$validation_log2fc) ~ "same direction",
          TRUE ~ "opposite/unclear"
        )
      )
    p2 <- ggplot2::ggplot(scatter, ggplot2::aes(x = discovery_logFC, y = validation_log2fc)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(shape = in_top20), alpha = 0.75, size = 2) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "grey35") +
      ggplot2::labs(
        title = "Discovery effect size versus GSE207310 validation effect size",
        subtitle = "Positive validation log2FC indicates higher expression in the more severe validation group/model.",
        x = "GSE136103 discovery logFC (cirrhotic vs healthy)",
        y = "GSE207310 validation log2FC",
        shape = "Top-20 candidate"
      ) +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(file.path(fig_dir, "validation_discovery_vs_gse207310_logfc.pdf"), p2, width = 8, height = 6)
  }

  if (length(top_genes) > 0) {
    top_val <- val2 %>%
      dplyr::filter(.data$gene %in% top_genes) %>%
      dplyr::mutate(
        gene = factor(.data$gene, levels = rev(top_genes[top_genes %in% .data$gene])),
        neg_log10_padj = -log10(pmax(.data$validation_padj, 1e-300))
      )
    if (nrow(top_val) > 0) {
      p3 <- ggplot2::ggplot(top_val, ggplot2::aes(x = gene, y = validation_log2fc)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
        ggplot2::geom_col(width = 0.75) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "GSE207310 validation signal for top discovery biomarkers",
          x = NULL,
          y = "Validation log2FC",
          caption = "Genes shown are top ranked discovery biomarkers that were detected in GSE207310."
        ) +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(file.path(fig_dir, "validation_top_biomarkers_logfc.pdf"), p3, width = 8, height = max(4, 0.28 * nrow(top_val) + 2))
    }
  }

  invisible(NULL)
}

fallback_figure_gallery <- function() {
  fig_dir <- "results/figures"
  if (!dir.exists(fig_dir)) return("<p><em>No figures available.</em></p>")
  figs <- list.files(fig_dir, pattern = "\\.(pdf|png|jpg|jpeg|svg)$", full.names = TRUE, ignore.case = TRUE)
  preferred <- file.path(fig_dir, c(
    "qc_violin_prefilter.pdf", "qc_violin_postfilter.pdf", "qc_scatter_nFeature_nCount.pdf",
    "umap_disease_sample_cluster.pdf", "umap_cell_type.pdf", "umap_compartment.pdf",
    "dotplot_major_celltype_markers.pdf", "dotplot_required_compartments.pdf",
    "cell_type_composition_by_donor.pdf", "cell_state_scores_required_compartments.pdf",
    "pathway_top_terms_by_compartment.pdf", "ranked_biomarkers_top20.pdf",
    "validation_discovery_overlap.pdf", "validation_discovery_vs_gse207310_logfc.pdf",
    "validation_top_biomarkers_logfc.pdf"
  ))
  figs <- unique(c(preferred[file.exists(preferred)], figs))
  if (length(figs) == 0) return("<p><em>No figures available.</em></p>")
  titles <- tools::file_path_sans_ext(basename(figs))
  titles <- tools::toTitleCase(gsub("_", " ", titles))
  paste(vapply(seq_along(figs), function(i) fallback_figure_embed(figs[[i]], titles[[i]]), character(1)), collapse = "\n")
}

write_fallback_html_report <- function(output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  top <- read_csv_report("results/tables/ranked_biomarkers_top20.csv")
  de <- read_csv_report("results/tables/pseudobulk_de_required_compartments.csv")
  fgsea <- read_csv_report("results/tables/fgsea_all_results.csv")
  qc <- read_csv_report("results/tables/qc_summary_by_sample.csv")
  libs <- read_csv_report("results/tables/gse136103_library_summary.csv")
  val <- read_csv_report("results/tables/gse207310_validation_results.csv")
  val_meta <- read_csv_report("results/tables/gse207310_curated_metadata.csv")
  val_map <- read_csv_report("results/tables/gse207310_expression_file_mapping.csv")
  val_matrix_summary <- read_csv_report("results/tables/gse207310_expression_matrix_summary.csv")
  val_id_map <- read_csv_report("results/tables/gse207310_gene_id_mapping.csv")
  val_overlap <- read_csv_report("results/tables/gse207310_validation_overlap_summary.csv")

  make_validation_figures_from_tables()
  val_overlap <- read_csv_report("results/tables/gse207310_validation_overlap_summary.csv")

  n_libraries <- if (nrow(libs) > 0) nrow(libs) else 0
  n_samples <- if (nrow(qc) > 0 && "sample_id" %in% names(qc)) length(unique(qc$sample_id)) else 0
  n_cells_pre <- sum_first_existing_col(qc, c("n_cells_prefilter", "n_cells_pre_qc", "cells_before_qc", "n_cells_before_qc"))
  n_cells_post <- sum_first_existing_col(qc, c("n_cells_postfilter", "n_cells_post_qc", "cells_retained", "n_cells_after_qc"))
  pct_retained <- if (isTRUE(n_cells_pre > 0)) round(100 * n_cells_post / n_cells_pre, 1) else 0
  n_validation_results <- if (nrow(val) > 0) nrow(val) else 0
  validation_status <- if (nrow(val) > 0) {
    "GSE207310 validation results were incorporated into biomarker prioritization."
  } else if (nrow(val_meta) > 0 || nrow(val_map) > 0 || nrow(val_matrix_summary) > 0) {
    "GSE207310 expression or metadata files were parsed, but no usable validation contrast was available. Use config/gse207310_validation_metadata_override.csv if phenotype labels need manual harmonization."
  } else {
    "GSE207310 validation was skipped or unavailable in this run."
  }

  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>Karyon Bio Liver Fibrosis Report</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;max-width:1200px;margin:40px auto;line-height:1.55;color:#0f172a;padding:0 20px}",
    "h1{font-size:2.1rem}h2{border-top:1px solid #e5e7eb;padding-top:24px;margin-top:32px}",
    ".hero{background:linear-gradient(135deg,#0f172a,#1e3a8a,#0f766e);color:white;border-radius:18px;padding:24px;margin-bottom:24px}",
    ".metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin:18px 0 22px}",
    ".metric-card{border:1px solid #e5e7eb;border-radius:16px;padding:16px;background:white;box-shadow:0 6px 18px rgba(15,23,42,.06)}.metric-value{font-size:1.7rem;font-weight:800;color:#1e3a8a}.metric-label{color:#475569}",
    "table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #e5e7eb;padding:6px 8px}th{background:#f1f5f9;text-align:left}",
    ".figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px}.report-figure{border:1px solid #e5e7eb;border-radius:16px;padding:12px;background:white;box-shadow:0 4px 16px rgba(15,23,42,.05)}.report-figure img{width:100%;height:auto;border-radius:10px}.report-figure figcaption{font-size:.88rem;color:#475569;margin-top:8px}",
    ".card{border-left:5px solid #2563eb;background:#eff6ff;border-radius:12px;padding:16px;margin:16px 0}.toc{border:1px solid #e5e7eb;background:#f8fafc;border-radius:16px;padding:16px 22px;margin:20px 0}.toc ol{margin-bottom:0;padding-left:1.4rem}.toc li{margin:.25rem 0}.toc a{text-decoration:none;color:#1e3a8a}.pdf-link,.pdf-inline{text-align:center;color:#475569}.report-figure img{width:100%;height:auto;border-radius:10px}</style></head><body>",
    "<div class='hero'><h1>Cell-Type-Specific Biomarker Discovery in Human Liver Fibrosis</h1><p>Automated HTML report generated by the Karyon Bio pipeline.</p></div>",
    "<div class='metric-grid'>", fallback_metric_card(n_libraries, "curated primary libraries"), fallback_metric_card(n_samples, "samples with QC summaries"), fallback_metric_card(n_cells_pre, "cells before QC"), fallback_metric_card(n_cells_post, "cells retained post-QC"), fallback_metric_card(paste0(pct_retained, "%"), "overall QC retention"), fallback_metric_card(nrow(de), "required-compartment DE rows"), fallback_metric_card(nrow(fgsea), "pathway results"), fallback_metric_card(nrow(top), "ranked top candidates"), fallback_metric_card(n_validation_results, "GSE207310 validation rows"), "</div>",
    "<div class='toc'><strong>Table of contents</strong><ol><li><a href='#executive-interpretation'>Executive interpretation</a></li><li><a href='#dataset-curation-and-metadata'>Dataset curation and metadata</a></li><li><a href='#qc-and-preprocessing'>QC and preprocessing</a></li><li><a href='#cell-type-annotation-and-required-compartments'>Cell type annotation and required compartments</a></li><li><a href='#fibrosiscirrhosis-associated-genes-and-cell-states'>Fibrosis/cirrhosis-associated genes and cell states</a></li><li><a href='#pathway-and-mechanism-analysis'>Pathway and mechanism analysis</a></li><li><a href='#optional-gse207310-validation'>Optional GSE207310 validation</a></li><li><a href='#biomarker-and-therapeutic-target-prioritization'>Biomarker and therapeutic target prioritization</a></li><li><a href='#translational-interpretation'>Translational interpretation</a></li><li><a href='#limitations-and-next-validation-steps'>Limitations and next validation steps</a></li><li><a href='#reproducibility-and-deliverables'>Reproducibility and deliverables</a></li></ol></div>",
    "<div class='card'><strong>Workflow:</strong> GSE136103 discovery with optional GSE207310 validation; QC, annotation, pseudobulk DE, pathway analysis, and biomarker prioritization.</div>",
    "<h2 id='executive-interpretation'>Executive interpretation</h2><p>The core discovery analysis completed on GSE136103. The ranked candidates combine donor-aware pseudobulk DE, pathway mechanism membership, cell-type specificity, translational annotation, and optional validation when available.</p>",
    "<h2 id='dataset-curation-and-metadata'>Dataset curation and metadata</h2>", html_table_report(libs, 30),
    "<h2 id='qc-and-preprocessing'>QC and preprocessing</h2>", html_table_report(qc, 30),
    "<h2>Key figures</h2><div class='figure-grid'>", fallback_figure_gallery(), "</div>",
    "<h2 id='cell-type-annotation-and-required-compartments'>Cell type annotation and required compartments</h2><p>Cell-type annotations are supported by UMAP and marker dot plot figures above.</p>",
    "<div class='card'><strong>Fibrosis program score definition:</strong> The Cell State Scores Required Compartments figure uses Seurat AddModuleScore on curated fibrosis mechanism gene sets from resources/compartment_gene_sets.yml. For each cell, the score is the average normalized expression of program genes minus matched control genes binned by expression level. Higher scores indicate stronger relative activity of that program. The updated box plots are filled by disease status: healthy versus cirrhotic.</div>",
    "<h2 id='fibrosiscirrhosis-associated-genes-and-cell-states'>Fibrosis/cirrhosis-associated genes and cell states</h2>", html_table_report(de, 30),
    "<h2 id='pathway-and-mechanism-analysis'>Pathway and mechanism analysis</h2>", html_table_report(fgsea, 30),
    "<h2 id='optional-gse207310-validation'>Optional GSE207310 validation</h2><p>GSE207310 is treated as an optional bulk liver RNA-seq validation dataset. The parser supports combined expression matrices and per-sample gzipped TXT/TSV/CSV files from GSE207310_RAW.tar. It infers gene and expression columns, maps numeric Entrez IDs or Ensembl IDs to HGNC gene symbols using org.Hs.eg.db when possible, uses DESeq2 for integer-like counts, and falls back to limma for normalized expression values. Because this dataset is bulk tissue, it provides whole-liver disease-association support but not cell-type-specific localization.</p><div class='card'><strong>Validation status:</strong> ", html_escape_simple(validation_status), "</div><h3>Validation parsing audit</h3>", html_table_report(val_matrix_summary, 20), "<h4>Gene identifier mapping audit</h4>", html_table_report(val_id_map, 30), html_table_report(val_map, 30), html_table_report(val_meta, 30), "<h3>Discovery-validation overlap statistics</h3>", html_table_report(val_overlap, 20), "<h3>Validation results</h3>", html_table_report(val, 30),
    "<h2 id='biomarker-and-therapeutic-target-prioritization'>Biomarker and therapeutic target prioritization</h2>", html_table_report(top, 20),
    "<h2 id='translational-interpretation'>Translational interpretation</h2><p>Prioritize disease-upregulated, secreted/surface/ECM candidates in HSC/mesenchymal, macrophage/monocyte, and endothelial compartments.</p>",
    "<h2 id='limitations-and-next-validation-steps'>Limitations and next validation steps</h2><p>Validate top candidates with orthogonal localization assays, spatial transcriptomics, immunostaining, and blood/tissue biomarker assays.</p>",
    "<h2 id='reproducibility-and-deliverables'>Reproducibility and deliverables</h2><p>Primary deliverables are written to results/reports, results/tables, and results/figures.</p>",
    "</body></html>"
  )
  writeLines(html, output_path)
  invisible(output_path)
}

render_all_reports <- function(biomarker_result, config = read_config()) {
  message_header("Rendering final reports")
  dir.create(config$paths$report_dir, recursive = TRUE, showWarnings = FALSE)

  final <- render_quarto_if_available(
    "reports/final_report_template.qmd",
    output_dir = config$paths$report_dir,
    output_file = "final_report.html"
  )

  exec <- render_quarto_if_available(
    "reports/executive_summary_template.qmd",
    output_dir = config$paths$report_dir,
    output_file = "executive_summary.html"
  )

  manifest <- tibble::tibble(
    report = c("final_report", "executive_summary"),
    path = c(final, exec),
    generated_at = as.character(Sys.time())
  )
  readr::write_csv(manifest, file.path(config$paths$report_dir, "report_manifest.csv"))
  manifest
}
