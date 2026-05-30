`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

read_config <- function(path = "config/config.yml") {
  yaml::read_yaml(path)
}

project_dirs <- function(config = read_config()) {
  unlist(config$paths, use.names = TRUE)
}

init_project_dirs <- function(config = read_config()) {
  dirs <- project_dirs(config)
  purrr::walk(dirs, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))
  invisible(dirs)
}

write_csv_safely <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(x) || tibble::is_tibble(x)) x <- sanitize_table(x)
  readr::write_csv(x, path)
  path
}

save_rds_safely <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
  path
}

save_plot_safely <- function(plot, path, width = 9, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, device = cairo_pdf)
  path
}

save_plot_png_safely <- function(plot, path, width = 9, height = 6, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi)
  path
}

clean_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x
}


# Normalize validation gene identifiers to HGNC-style gene symbols when possible.
# GSE207310 per-sample files may use numeric Entrez IDs or Ensembl IDs, while
# GSE136103 discovery results use gene symbols. Without this mapping, validation
# can appear to have no overlap with discovery even when the biology overlaps.
normalize_gene_identifiers <- function(x, source = "validation", write_mapping = FALSE, mapping_path = NULL) {
  raw <- as.character(x)
  raw <- trimws(raw)
  raw <- sub("\\.\\d+$", "", raw) # remove Ensembl version suffixes
  raw[raw == ""] <- NA_character_

  cleaned <- clean_gene_symbols(raw)
  mapped <- rep(NA_character_, length(raw))
  id_type <- rep("symbol_or_unmapped", length(raw))

  if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    nonmissing <- unique(raw[!is.na(raw)])

    numeric_keys <- nonmissing[grepl("^[0-9]+$", nonmissing)]
    if (length(numeric_keys) > 0) {
      numeric_map <- tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = numeric_keys,
          keytype = "ENTREZID",
          column = "SYMBOL",
          multiVals = "first"
        ),
        error = function(e) NULL
      )
      if (!is.null(numeric_map)) {
        idx <- raw %in% names(numeric_map)
        mapped[idx] <- unname(numeric_map[raw[idx]])
        id_type[idx] <- "ENTREZID"
      }
    }

    ensembl_keys <- nonmissing[grepl("^ENSG[0-9]+", nonmissing, ignore.case = TRUE)]
    if (length(ensembl_keys) > 0) {
      ensembl_map <- tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = ensembl_keys,
          keytype = "ENSEMBL",
          column = "SYMBOL",
          multiVals = "first"
        ),
        error = function(e) NULL
      )
      if (!is.null(ensembl_map)) {
        idx <- raw %in% names(ensembl_map)
        mapped[idx] <- unname(ensembl_map[raw[idx]])
        id_type[idx] <- "ENSEMBL"
      }
    }

    # For already-symbol-like keys, keep the cleaned symbol. This preserves
    # discovery data and validation files that already provide HGNC symbols.
    symbol_like <- is.na(mapped) & !is.na(cleaned) & !grepl("^[0-9]+$", raw) & !grepl("^ENSG[0-9]+", raw, ignore.case = TRUE)
    mapped[symbol_like] <- cleaned[symbol_like]
    id_type[symbol_like] <- "SYMBOL"
  } else {
    mapped <- cleaned
  }

  mapped <- clean_gene_symbols(mapped)
  out <- ifelse(!is.na(mapped) & mapped != "", mapped, cleaned)

  if (isTRUE(write_mapping) && !is.null(mapping_path)) {
    map_tbl <- tibble::tibble(
      original_gene_id = raw,
      normalized_gene_symbol = out,
      inferred_id_type = id_type
    ) %>%
      dplyr::distinct() %>%
      dplyr::arrange(.data$inferred_id_type, .data$original_gene_id)
    write_csv_safely(map_tbl, mapping_path)
  }
  out
}


rescale_01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

mad_upper <- function(x, n_mad = 3, hard_max = Inf) {
  val <- stats::median(x, na.rm = TRUE) + n_mad * stats::mad(x, na.rm = TRUE)
  min(val, hard_max, na.rm = TRUE)
}

mad_lower <- function(x, n_mad = 3, hard_min = -Inf) {
  val <- stats::median(x, na.rm = TRUE) - n_mad * stats::mad(x, na.rm = TRUE)
  max(val, hard_min, na.rm = TRUE)
}

make_unique_clean_names <- function(x) {
  x <- janitor::make_clean_names(x)
  make.unique(x, sep = "_")
}

std_error <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

safe_factor <- function(x, levels = NULL, ref = NULL) {
  out <- factor(x, levels = levels)
  if (!is.null(ref) && ref %in% levels(out)) out <- stats::relevel(out, ref = ref)
  out
}

looks_like_gene_column <- function(x) {
  x <- as.character(x)
  mean(grepl("^[A-Za-z0-9_.-]+$", x) & nchar(x) <= 40, na.rm = TRUE) > 0.8
}

safe_neg_log10 <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- -log10(pmax(p[ok], .Machine$double.xmin))
  out
}

message_header <- function(...) {
  msg <- paste0(..., collapse = "")
  message("\n", paste(rep("=", nchar(msg)), collapse = ""))
  message(msg)
  message(paste(rep("=", nchar(msg)), collapse = ""))
}

get_assay_data_compat <- function(obj, assay = "RNA", layer = "counts") {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = layer)
  )
}

# Convert list-like columns in plain data frames/tibbles to atomic columns.
# GEOquery/Bioconductor/Seurat can create list columns, which break dplyr,
# readr::write_csv(), edgeR design matrices, and ggplot operations with errors
# like: "Argument 'x' is not a vector: list".
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

safe_write_ready <- function(x) {
  if (is.data.frame(x) || tibble::is_tibble(x)) sanitize_table(x) else x
}
