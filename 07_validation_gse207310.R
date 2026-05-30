# Orthogonal validation using GSE207310 bulk liver RNA-seq.
# This script is intentionally flexible because GEO supplementary TXT formats vary.
# v20 adds robust support for GSE207310_RAW.tar archives containing per-sample .txt.gz files.

read_txt_maybe_gz <- function(path) {
  # GEO supplementary TXT files vary: some are tab-delimited matrices, some are
  # per-sample RSEM/featureCounts/HTSeq-style tables, and some include comment
  # lines. Try several safe readers before giving up.
  out <- tryCatch(
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, comment = "#", name_repair = "unique_quiet"),
    error = function(e) NULL
  )
  if (!is.null(out) && ncol(out) >= 2) return(sanitize_table(out))

  out <- tryCatch(
    readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE,
                      comment = "#", col_names = FALSE, name_repair = "unique_quiet"),
    error = function(e) NULL
  )
  if (!is.null(out) && ncol(out) >= 2) return(sanitize_table(out))

  out <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE, comment = "#", name_repair = "unique_quiet"),
    error = function(e) NULL
  )
  if (!is.null(out) && ncol(out) >= 2) return(sanitize_table(out))

  stop("Could not read validation TXT/TSV/CSV file: ", path, call. = FALSE)
}

looks_like_ensembl <- function(x) {
  x <- as.character(x)
  mean(grepl("^ENSG[0-9]+", x), na.rm = TRUE) > 0.3
}

infer_gene_column <- function(tbl) {
  tbl <- sanitize_table(tbl)
  names(tbl) <- make_unique_clean_names(names(tbl))

  # Prefer gene-symbol/name columns over Ensembl ID columns so validation can
  # match candidate gene symbols from GSE136103.
  priority_patterns <- c(
    "^gene_name$", "^gene_symbol$", "^symbol$", "^external_gene_name$",
    "^hgnc_symbol$", "^gene$", "^geneid$", "^gene_id$", "ensembl"
  )
  for (pat in priority_patterns) {
    hits <- names(tbl)[grepl(pat, names(tbl), ignore.case = TRUE)]
    if (length(hits) > 0) return(hits[1])
  }

  non_numeric <- names(tbl)[!vapply(tbl, is.numeric, logical(1))]
  gene_like <- non_numeric[vapply(tbl[non_numeric], looks_like_gene_column, logical(1))]
  if (length(gene_like) > 0) return(gene_like[1])

  # HTSeq-style two-column files often arrive without headers. In that case,
  # the first column is usually the gene identifier and the second is the count.
  if (ncol(tbl) >= 2) return(names(tbl)[1])

  stop("Could not infer gene column in validation count/expression table.", call. = FALSE)
}

infer_expression_column <- function(tbl, gene_col = NULL) {
  tbl <- sanitize_table(tbl)
  names(tbl) <- make_unique_clean_names(names(tbl))
  if (!is.null(gene_col) && gene_col %in% names(tbl)) exclude <- gene_col else exclude <- character()

  # Parse numeric-looking character columns. This helps files where counts are
  # read as character because of occasional non-numeric rows.
  candidate_cols <- setdiff(names(tbl), exclude)
  for (cc in candidate_cols) {
    if (!is.numeric(tbl[[cc]]) && !is.factor(tbl[[cc]])) {
      parsed <- suppressWarnings(readr::parse_number(as.character(tbl[[cc]])))
      if (sum(!is.na(parsed)) > 0.7 * length(parsed)) tbl[[cc]] <- parsed
    }
  }

  numeric_cols <- setdiff(names(tbl)[vapply(tbl, is.numeric, logical(1))], exclude)
  if (length(numeric_cols) == 0) return(NULL)

  # Prefer count-like columns. Avoid genomic coordinate / length columns.
  avoid <- grepl("start|end|width|length|chr|chrom|strand|gc|eff_length|effective_length", numeric_cols, ignore.case = TRUE)
  usable <- numeric_cols[!avoid]
  if (length(usable) == 0) usable <- numeric_cols

  priorities <- c(
    "^expected_count$", "^expected_counts$", "^raw_count$", "^raw_counts$",
    "^count$", "^counts$", "read_count", "numreads", "unstranded",
    "tpm$", "fpkm$", "rpkm$"
  )
  for (pat in priorities) {
    hits <- usable[grepl(pat, usable, ignore.case = TRUE)]
    if (length(hits) > 0) return(hits[1])
  }

  # FeatureCounts often uses the last sample-specific numeric column after
  # Geneid/Chr/Start/End/Strand/Length. Otherwise, choose the column with the
  # largest informative sum.
  sums <- vapply(tbl[usable], function(x) sum(as.numeric(x), na.rm = TRUE), numeric(1))
  usable[which.max(sums)]
}

standardize_count_matrix <- function(tbl) {
  tbl <- sanitize_table(tbl)
  names(tbl) <- make_unique_clean_names(names(tbl))
  gene_col <- infer_gene_column(tbl)

  candidate_cols <- setdiff(names(tbl), gene_col)
  for (cc in candidate_cols) {
    if (!is.numeric(tbl[[cc]]) && !is.factor(tbl[[cc]])) {
      parsed <- suppressWarnings(readr::parse_number(as.character(tbl[[cc]])))
      if (sum(!is.na(parsed)) > 0.7 * length(parsed)) tbl[[cc]] <- parsed
    }
  }

  numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, gene_col)
  numeric_cols <- numeric_cols[!grepl("start|end|width|length|chr|chrom|strand|gc|eff_length|effective_length", numeric_cols, ignore.case = TRUE)]
  if (length(numeric_cols) == 0) stop("Could not infer sample expression columns in count table.", call. = FALSE)

  mat_tbl <- tbl %>%
    mutate(gene = normalize_gene_identifiers(.data[[gene_col]], write_mapping = TRUE, mapping_path = file.path(read_config()$paths$table_dir, "gse207310_gene_id_mapping.csv"))) %>%
    filter(!is.na(gene), gene != "", !grepl("^__", gene)) %>%
    group_by(gene) %>%
    summarise(across(all_of(numeric_cols), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

  out <- as.matrix(mat_tbl[, numeric_cols, drop = FALSE])
  rownames(out) <- mat_tbl$gene
  storage.mode(out) <- "numeric"
  out
}

sample_id_from_validation_file <- function(path) {
  base <- basename(path)
  gsm <- stringr::str_extract(base, "GSM[0-9]+")
  if (!is.na(gsm) && nzchar(gsm)) return(gsm)
  # Remove common compression/table suffixes.
  base <- sub("\\.gz$", "", base, ignore.case = TRUE)
  base <- sub("\\.(txt|tsv|csv)$", "", base, ignore.case = TRUE)
  make_unique_clean_names(base)
}

combine_per_sample_txt_files <- function(txt_files) {
  parsed <- purrr::map(txt_files, function(f) {
    tbl <- tryCatch(read_txt_maybe_gz(f), error = function(e) {
      warning("Skipping unreadable validation file ", basename(f), ": ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)

    tbl <- sanitize_table(tbl)
    names(tbl) <- make_unique_clean_names(names(tbl))
    gene_col <- tryCatch(infer_gene_column(tbl), error = function(e) NULL)
    if (is.null(gene_col)) return(NULL)

    count_col <- infer_expression_column(tbl, gene_col)
    if (is.null(count_col)) return(NULL)

    tibble::tibble(
      gene = normalize_gene_identifiers(tbl[[gene_col]], write_mapping = FALSE),
      original_gene_id = as.character(tbl[[gene_col]]),
      sample_id = sample_id_from_validation_file(f),
      count = suppressWarnings(as.numeric(tbl[[count_col]])),
      source_file = basename(f),
      expression_column = count_col
    ) %>%
      filter(!is.na(gene), gene != "", !grepl("^__", gene), !is.na(count)) %>%
      group_by(gene, original_gene_id, sample_id, source_file, expression_column) %>%
      summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  }) %>% purrr::compact()

  if (length(parsed) == 0) return(NULL)
  long <- dplyr::bind_rows(parsed) %>% sanitize_table()

  # Record exactly which file/column contributed to each validation sample.
  mapping <- long %>%
    dplyr::distinct(sample_id, source_file, expression_column) %>%
    dplyr::arrange(sample_id)
  write_csv_safely(mapping, file.path(read_config()$paths$table_dir, "gse207310_expression_file_mapping.csv"))

  gene_mapping <- long %>%
    dplyr::distinct(original_gene_id, gene) %>%
    dplyr::rename(normalized_gene_symbol = gene) %>%
    dplyr::mutate(
      inferred_id_type = dplyr::case_when(
        grepl("^[0-9]+$", .data$original_gene_id) ~ "ENTREZID_or_numeric",
        grepl("^ENSG[0-9]+", .data$original_gene_id, ignore.case = TRUE) ~ "ENSEMBL",
        TRUE ~ "SYMBOL_or_other"
      )
    ) %>%
    dplyr::arrange(.data$inferred_id_type, .data$original_gene_id)
  write_csv_safely(gene_mapping, file.path(read_config()$paths$table_dir, "gse207310_gene_id_mapping.csv"))

  wide <- long %>%
    dplyr::select(gene, sample_id, count) %>%
    dplyr::group_by(gene, sample_id) %>%
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = count, values_fill = 0)

  mat <- as.matrix(wide[, setdiff(names(wide), "gene"), drop = FALSE])
  rownames(mat) <- wide$gene
  storage.mode(mat) <- "numeric"
  mat
}

read_gse207310_counts <- function(downloaded_files, config = read_config()) {
  message_header("Reading GSE207310 count/expression tables")
  txt_files <- downloaded_files$path[grepl("\\.txt(\\.gz)?$|\\.tsv(\\.gz)?$|\\.csv(\\.gz)?$", downloaded_files$path, ignore.case = TRUE)]
  txt_files <- txt_files[!grepl("series_matrix", basename(txt_files), ignore.case = TRUE)]
  txt_files <- unique(txt_files[file.exists(txt_files)])

  if (length(txt_files) == 0) {
    warning("No GSE207310 TXT/TSV/CSV supplementary count files found.")
    return(NULL)
  }

  message("Found ", length(txt_files), " candidate validation expression file(s).")

  # Case 1: a single count/expression matrix with genes x samples.
  if (length(txt_files) == 1) {
    tbl <- read_txt_maybe_gz(txt_files[1])
    mat <- tryCatch(standardize_count_matrix(tbl), error = function(e) {
      warning("Single validation matrix parsing failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (!is.null(mat) && ncol(mat) >= 2) {
      write_csv_safely(
        tibble::tibble(sample_id = colnames(mat), source_file = basename(txt_files[1]), expression_column = colnames(mat)),
        file.path(config$paths$table_dir, "gse207310_expression_file_mapping.csv")
      )
      return(mat)
    }
  }

  # Case 2: one gzipped TXT/TSV/CSV file per sample. This is the common layout
  # for GSE207310_RAW.tar and similar GEO bulk RNA-seq submissions.
  mat <- combine_per_sample_txt_files(txt_files)
  if (is.null(mat) || ncol(mat) == 0) {
    warning("Could not combine GSE207310 per-sample TXT files into an expression matrix.")
    return(NULL)
  }

  id_map_path <- file.path(config$paths$table_dir, "gse207310_gene_id_mapping.csv")
  mapped_summary <- tibble::tibble(metric = c("genes", "samples"), value = c(nrow(mat), ncol(mat)))
  if (file.exists(id_map_path)) {
    id_map <- tryCatch(readr::read_csv(id_map_path, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(id_map) && nrow(id_map) > 0 && "inferred_id_type" %in% names(id_map)) {
      extra <- id_map %>%
        dplyr::count(.data$inferred_id_type, name = "value") %>%
        dplyr::mutate(metric = paste0("gene_id_type__", .data$inferred_id_type)) %>%
        dplyr::select(metric, value)
      mapped_summary <- dplyr::bind_rows(mapped_summary, extra)
    }
  }
  write_csv_safely(mapped_summary, file.path(config$paths$table_dir, "gse207310_expression_matrix_summary.csv"))
  mat
}

extract_gse207310_validation_metadata <- function(geo_metadata, config = read_config()) {
  message_header("Curating GSE207310 validation metadata")
  pheno <- sanitize_table(geo_metadata) %>% mutate(sample_id = geo_accession)

  text_cols <- names(pheno)[vapply(pheno, function(x) is.character(x) || is.factor(x), logical(1))]
  pheno$all_metadata_text <- apply(pheno[, text_cols, drop = FALSE], 1, function(x) paste(x, collapse = "; "))
  lower_text <- tolower(pheno$all_metadata_text)

  pheno$validation_group <- dplyr::case_when(
    grepl("no[-_ ]?nafld|non[-_ ]?nafld|control|healthy|no nash", lower_text) ~ "no_NAFLD_or_no_NASH",
    grepl("nash|mash|steatohepatitis", lower_text) ~ "NASH_or_MASH",
    grepl("nafld", lower_text) ~ "NAFLD_unspecified",
    TRUE ~ NA_character_
  )

  severity_cols <- names(pheno)[grepl("fibrosis|nas|saf|steatosis|activity|balloon|inflam", names(pheno), ignore.case = TRUE)]
  pheno$severity_numeric <- NA_real_
  if (length(severity_cols) > 0) {
    # Use the first numeric-looking severity column as default. Users can override in config/gse207310_validation_metadata_override.csv.
    for (cc in severity_cols) {
      val <- readr::parse_number(as.character(pheno[[cc]]))
      if (sum(!is.na(val)) >= 5) {
        pheno$severity_numeric <- val
        pheno$severity_source_column <- cc
        break
      }
    }
  }

  override_file <- config$validation$metadata_override_file
  if (file.exists(override_file)) {
    override <- readr::read_csv(override_file, show_col_types = FALSE) %>% filter(!is.na(sample_id), sample_id != "")
    if (nrow(override) > 0) {
      pheno <- pheno %>%
        left_join(override, by = "sample_id", suffix = c("", "_override")) %>%
        mutate(
          validation_group = ifelse(!is.na(validation_group_override) & validation_group_override != "", validation_group_override, validation_group),
          severity_numeric = ifelse(!is.na(severity_numeric_override), severity_numeric_override, severity_numeric)
        )
    }
  }

  write_csv_safely(pheno, file.path(config$paths$table_dir, "gse207310_curated_metadata.csv"))
  pheno
}

run_limma_validation <- function(expr, pheno, model_label, contrast_label) {
  # Fallback for normalized TPM/FPKM-like matrices or non-integer values.
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  expr <- log2(expr + 1)

  if (model_label == "limma_group") {
    pheno$validation_group <- factor(pheno$validation_group)
    design <- stats::model.matrix(~ validation_group, data = pheno)
    fit <- limma::lmFit(expr, design)
    fit <- limma::eBayes(fit)
    coef_id <- grep("validation_group", colnames(design))[1]
    if (is.na(coef_id)) return(tibble::tibble())
  } else {
    design <- stats::model.matrix(~ severity_numeric, data = pheno)
    fit <- limma::lmFit(expr, design)
    fit <- limma::eBayes(fit)
    coef_id <- "severity_numeric"
  }

  limma::topTable(fit, coef = coef_id, number = Inf, sort.by = "P") %>%
    sanitize_table() %>%
    tibble::rownames_to_column("gene") %>%
    as_tibble() %>%
    transmute(
      gene = normalize_gene_identifiers(gene),
      baseMean = NA_real_,
      log2FoldChange = logFC,
      lfcSE = NA_real_,
      stat = t,
      pvalue = P.Value,
      padj = adj.P.Val,
      validation_model = model_label,
      validation_contrast = contrast_label,
      validation_direction = dplyr::case_when(log2FoldChange > 0 ~ "up_in_more_severe", log2FoldChange < 0 ~ "down_in_more_severe", TRUE ~ "flat")
    )
}

run_deseq2_validation <- function(counts, pheno, candidate_genes = NULL, config = read_config()) {
  if (is.null(counts) || ncol(counts) < 4) {
    warning("GSE207310 counts unavailable or too few samples; validation will be skipped.")
    return(tibble::tibble())
  }

  # Match sample names either by GSM accession or exact column names.
  pheno <- sanitize_table(pheno) %>% mutate(sample_id = as.character(sample_id))
  col_match <- intersect(colnames(counts), pheno$sample_id)
  if (length(col_match) < 4) {
    # Try matching GSM inside column names.
    col_gsm <- stringr::str_extract(colnames(counts), "GSM[0-9]+")
    names(col_gsm) <- colnames(counts)
    matched_cols <- names(col_gsm)[col_gsm %in% pheno$sample_id]
    if (length(matched_cols) >= 4) {
      colnames(counts)[match(matched_cols, colnames(counts))] <- col_gsm[matched_cols]
      col_match <- intersect(colnames(counts), pheno$sample_id)
    }
  }

  if (length(col_match) < 4) {
    warning("Could not match GSE207310 expression columns to GEO sample metadata; validation skipped.")
    write_csv_safely(
      tibble::tibble(expression_column = colnames(counts), parsed_gsm = stringr::str_extract(colnames(counts), "GSM[0-9]+")),
      file.path(config$paths$table_dir, "gse207310_unmatched_expression_columns.csv")
    )
    return(tibble::tibble())
  }

  counts <- counts[, col_match, drop = FALSE]
  pheno <- pheno %>% filter(sample_id %in% col_match) %>% arrange(match(sample_id, colnames(counts)))
  stopifnot(identical(pheno$sample_id, colnames(counts)))

  keep <- rowSums(counts >= 1, na.rm = TRUE) >= max(3, floor(ncol(counts) * 0.2))
  counts <- counts[keep, , drop = FALSE]

  if (!is.null(candidate_genes)) {
    candidate_genes <- clean_gene_symbols(candidate_genes)
  }

  # Determine whether these values look like raw integer counts. If not, use a
  # limma log-expression fallback rather than forcing DESeq2 on TPM/FPKM values.
  values <- as.numeric(counts[seq_len(min(length(counts), 10000))])
  is_integer_like <- mean(abs(values - round(values)) < 1e-6, na.rm = TRUE) > 0.95
  use_deseq <- is_integer_like && requireNamespace("DESeq2", quietly = TRUE)

  if (sum(!is.na(pheno$validation_group)) >= 6 && length(unique(na.omit(pheno$validation_group))) >= 2) {
    pheno2 <- pheno %>% filter(!is.na(validation_group))
    counts2 <- counts[, pheno2$sample_id, drop = FALSE]
    pheno2$validation_group <- factor(pheno2$validation_group)
    ref <- grep("no_NAFLD|no_NASH|control|healthy", levels(pheno2$validation_group), value = TRUE)[1]
    if (!is.na(ref)) pheno2$validation_group <- stats::relevel(pheno2$validation_group, ref = ref)

    if (!use_deseq) {
      return(run_limma_validation(counts2, pheno2, "limma_group", "validation_group"))
    }

    pheno2 <- as.data.frame(pheno2)
    rownames(pheno2) <- pheno2$sample_id
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts2), colData = pheno2, design = ~ validation_group)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    coef_name <- DESeq2::resultsNames(dds)[grepl("validation_group", DESeq2::resultsNames(dds))][1]
    res <- DESeq2::results(dds, name = coef_name)
    tbl <- sanitize_table(as.data.frame(res)) %>%
      tibble::rownames_to_column("gene") %>%
      as_tibble() %>%
      mutate(
        gene = normalize_gene_identifiers(gene),
        validation_model = "DESeq2_group",
        validation_contrast = coef_name,
        validation_direction = dplyr::case_when(log2FoldChange > 0 ~ "up_in_more_severe", log2FoldChange < 0 ~ "down_in_more_severe", TRUE ~ "flat")
      )
    return(tbl)
  }

  if (sum(!is.na(pheno$severity_numeric)) >= 6) {
    pheno2 <- pheno %>% filter(!is.na(severity_numeric))
    counts2 <- counts[, pheno2$sample_id, drop = FALSE]

    if (!use_deseq) {
      return(run_limma_validation(counts2, pheno2, "limma_continuous_severity", "severity_numeric"))
    }

    pheno2 <- as.data.frame(pheno2)
    rownames(pheno2) <- pheno2$sample_id
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts2), colData = pheno2, design = ~ severity_numeric)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, name = "severity_numeric")
    tbl <- sanitize_table(as.data.frame(res)) %>%
      tibble::rownames_to_column("gene") %>%
      as_tibble() %>%
      mutate(
        gene = normalize_gene_identifiers(gene),
        validation_model = "DESeq2_continuous_severity",
        validation_contrast = "severity_numeric",
        validation_direction = dplyr::case_when(log2FoldChange > 0 ~ "positive_with_severity", log2FoldChange < 0 ~ "negative_with_severity", TRUE ~ "flat")
      )
    return(tbl)
  }

  warning("GSE207310 metadata lacks an inferable validation group or severity score. Fill config/gse207310_validation_metadata_override.csv and rerun.")
  tibble::tibble()
}

run_gse207310_validation <- function(downloaded_files, geo_metadata, de_result, config = read_config()) {
  message_header("Running optional validation with GSE207310")
  counts <- read_gse207310_counts(downloaded_files, config)
  pheno <- extract_gse207310_validation_metadata(geo_metadata, config)

  candidate_genes <- sanitize_table(de_result$required_de) %>%
    arrange(FDR) %>%
    pull(gene) %>%
    unique() %>%
    head(config$validation$candidate_gene_limit)

  val <- run_deseq2_validation(counts, pheno, candidate_genes, config)
  if (nrow(val) > 0) {
    val <- sanitize_table(val) %>% arrange(padj, pvalue)
  }
  write_csv_safely(val, file.path(config$paths$table_dir, "gse207310_validation_results.csv"))

  template <- sanitize_table(pheno) %>% dplyr::select(sample_id, validation_group, severity_numeric, any_of("title"), all_metadata_text)
  write_csv_safely(template, file.path(config$paths$table_dir, "gse207310_validation_metadata_review_template.csv"))

  list(validation_results = sanitize_table(val), validation_metadata = sanitize_table(pheno))
}
