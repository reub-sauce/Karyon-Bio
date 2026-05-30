# Dataset download and metadata curation helpers.

get_geo_metadata <- function(geo_accession, out_dir = "results/tables", max_tries = 3) {
  message_header("Fetching GEO metadata: ", geo_accession)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cache_csv <- file.path(out_dir, paste0(tolower(geo_accession), "_geo_metadata.csv"))
  geo_cache_dir <- file.path("data", "raw", "geoquery_cache")
  dir.create(geo_cache_dir, recursive = TRUE, showWarnings = FALSE)

  read_cached <- function(reason = NULL) {
    if (file.exists(cache_csv)) {
      if (!is.null(reason)) message(reason, " Using cached metadata CSV: ", cache_csv)
      return(readr::read_csv(cache_csv, show_col_types = FALSE))
    }
    NULL
  }

  cached <- read_cached()
  # If a previous successful run wrote the metadata, prefer it and avoid NCBI calls.
  if (!is.null(cached) && nrow(cached) > 0) {
    message("Using cached GEO metadata for ", geo_accession, ": ", cache_csv)
    return(cached)
  }

  last_error <- NULL
  for (ii in seq_len(max_tries)) {
    gse <- tryCatch(
      GEOquery::getGEO(
        geo_accession,
        GSEMatrix = TRUE,
        AnnotGPL = FALSE,
        destdir = geo_cache_dir
      ),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(gse)) {
      if (inherits(gse, "ExpressionSet")) gse <- list(gse)

      pheno <- purrr::map_dfr(gse, function(eset) {
        pd <- as.data.frame(Biobase::pData(eset), check.names = FALSE)

        # Some GEO series matrices already include a geo_accession column.
        # Do not call rownames_to_column(var = "geo_accession") in that case,
        # because tibble requires unique column names and will error.
        rn <- rownames(pd)
        if ("geo_accession" %in% names(pd)) {
          pd$geo_accession_from_pheno <- as.character(pd$geo_accession)
          pd$geo_accession <- rn
        } else {
          pd <- tibble::rownames_to_column(pd, var = "geo_accession")
        }

        pd$platform_id <- Biobase::annotation(eset)
        pd
      })

      names(pheno) <- make_unique_clean_names(names(pheno))
      write_csv_safely(pheno, cache_csv)
      return(pheno)
    }

    if (ii < max_tries) {
      wait_seconds <- 5 * ii
      message("GEO metadata request failed for ", geo_accession, ": ", last_error)
      message("Retrying in ", wait_seconds, " seconds (attempt ", ii + 1, " of ", max_tries, ")...")
      Sys.sleep(wait_seconds)
    }
  }

  cached <- read_cached(paste0("GEO metadata request failed after ", max_tries, " attempt(s): ", last_error, "."))
  if (!is.null(cached)) return(cached)

  warning(
    "Could not fetch GEO metadata for ", geo_accession,
    ". Continuing with minimal metadata. This is sufficient for GSE136103 sample parsing; ",
    "GSE207310 validation may need a successful metadata fetch or a metadata override file.",
    call. = FALSE
  )
  pheno <- tibble::tibble(
    geo_accession = character(),
    title = character(),
    platform_id = character()
  )
  write_csv_safely(pheno, cache_csv)
  pheno
}

download_geo_supplementary <- function(geo_accession, raw_dir = "data/raw", force_download = FALSE) {
  message_header("Preparing GEO supplementary data: ", geo_accession)

  # Canonical project raw directory remains data/raw, but many users manually
  # create a top-level raw/ folder. Search both locations and copy root-level
  # GSE*.tar / .tgz files into the canonical per-GSE folder for untarring.
  canonical_raw_dir <- raw_dir
  alternate_raw_dirs <- unique(c(canonical_raw_dir, "raw"))
  alternate_raw_dirs <- alternate_raw_dirs[!is.na(alternate_raw_dirs) & nzchar(alternate_raw_dirs)]

  geo_dir <- file.path(canonical_raw_dir, geo_accession)
  dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

  existing_raw_dirs <- alternate_raw_dirs[dir.exists(alternate_raw_dirs)]
  if (length(existing_raw_dirs) == 0) existing_raw_dirs <- canonical_raw_dir

  root_pattern <- paste0("^", geo_accession)
  supp_pattern <- "\\.(tar|tgz|tar\\.gz|mtx|txt|tsv|csv)(\\.gz)?$"

  local_files <- unique(unlist(purrr::map(existing_raw_dirs, function(rd) {
    c(
      list.files(file.path(rd, geo_accession), recursive = TRUE, full.names = TRUE),
      list.files(rd, pattern = root_pattern, recursive = FALSE, full.names = TRUE, ignore.case = TRUE)
    )
  }), use.names = FALSE))
  local_files <- unique(local_files[file.exists(local_files)])

  has_supp <- length(local_files) > 0 && any(grepl(supp_pattern, local_files, ignore.case = TRUE))
  if (!has_supp || isTRUE(force_download)) {
    message("No local supplementary files found for ", geo_accession, "; downloading from GEO into ", canonical_raw_dir, ".")
    GEOquery::getGEOSuppFiles(geo_accession, baseDir = canonical_raw_dir, makeDirectory = TRUE)
  } else {
    message("Using existing local supplementary files for ", geo_accession, ".")
  }

  # Copy root-level archives from either data/raw/ or raw/ into data/raw/GSE.../.
  root_tar_files <- unique(unlist(purrr::map(existing_raw_dirs, function(rd) {
    list.files(rd, pattern = paste0("^", geo_accession, ".*\\.(tar|tgz|tar\\.gz)$"),
               recursive = FALSE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
  root_tar_files <- root_tar_files[file.exists(root_tar_files)]
  if (length(root_tar_files) > 0) {
    purrr::walk(root_tar_files, function(f) {
      dest <- file.path(geo_dir, basename(f))
      if (!file.exists(dest)) file.copy(f, dest, overwrite = FALSE)
    })
  }

  # Untar archives found in either canonical or top-level raw/GSE... folders.
  gse_dirs_to_scan <- unique(c(geo_dir, file.path(existing_raw_dirs, geo_accession)))
  gse_dirs_to_scan <- gse_dirs_to_scan[dir.exists(gse_dirs_to_scan)]
  tar_files <- unique(unlist(purrr::map(gse_dirs_to_scan, function(gd) {
    list.files(gd, pattern = "\\.tar$|\\.tar\\.gz$|\\.tgz$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))

  if (length(tar_files) > 0) {
    untar_dir <- file.path(geo_dir, "untar")
    dir.create(untar_dir, recursive = TRUE, showWarnings = FALSE)
    purrr::walk(tar_files, function(tf) {
      message("Untarring ", basename(tf))
      utils::untar(tf, exdir = untar_dir)
    })
  }

  # Return all files from canonical data/raw/GSE... and top-level raw/GSE...
  # so downstream code can use whichever location contains the MTX/TSV files.
  all_scan_dirs <- unique(c(geo_dir, file.path(existing_raw_dirs, geo_accession)))
  all_scan_dirs <- all_scan_dirs[dir.exists(all_scan_dirs)]
  files <- unique(unlist(purrr::map(all_scan_dirs, function(gd) {
    list.files(gd, recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE))
  files <- files[file.exists(files)]

  tibble::tibble(geo_accession = geo_accession, path = files, file = basename(files))
}

parse_gse136103_sample_id <- function(path) {
  base <- basename(path)
  base <- stringr::str_remove(base, "_matrix\\.mtx\\.gz$")
  base <- stringr::str_remove(base, "_matrix\\.mtx$")
  base
}

build_gse136103_sample_table <- function(downloaded_files, geo_metadata = NULL, config = read_config()) {
  message_header("Curating GSE136103 sample metadata")

  all_files <- downloaded_files$path
  matrix_files <- all_files[grepl("matrix\\.mtx(\\.gz)?$", all_files, ignore.case = TRUE)]
  if (length(matrix_files) == 0) {
    stop("No GSE136103 matrix.mtx files found. Check the GEO supplementary download.", call. = FALSE)
  }

  sample_table <- purrr::map_dfr(matrix_files, function(mtx) {
    prefix <- parse_gse136103_sample_id(mtx)
    dir <- dirname(mtx)
    gsm <- stringr::str_extract(prefix, "GSM[0-9]+")
    sample_label <- stringr::str_remove(prefix, "^GSM[0-9]+_")

    barcodes <- list.files(dir, pattern = paste0("^", prefix, ".*barcodes.*tsv(\\.gz)?$"), full.names = TRUE, ignore.case = TRUE)
    genes <- list.files(dir, pattern = paste0("^", prefix, ".*(genes|features).*tsv(\\.gz)?$"), full.names = TRUE, ignore.case = TRUE)

    disease_status <- dplyr::case_when(
      grepl("healthy", sample_label, ignore.case = TRUE) ~ "healthy",
      grepl("cirrhotic", sample_label, ignore.case = TRUE) ~ "cirrhotic",
      TRUE ~ NA_character_
    )

    species <- dplyr::case_when(
      grepl("mouse", sample_label, ignore.case = TRUE) ~ "mouse",
      grepl("blood", sample_label, ignore.case = TRUE) ~ "human",
      grepl("healthy|cirrhotic", sample_label, ignore.case = TRUE) ~ "human",
      TRUE ~ NA_character_
    )

    tissue <- dplyr::case_when(
      grepl("blood", sample_label, ignore.case = TRUE) ~ "blood",
      grepl("mouse", sample_label, ignore.case = TRUE) ~ "mouse_liver",
      grepl("healthy|cirrhotic", sample_label, ignore.case = TRUE) ~ "liver",
      TRUE ~ NA_character_
    )

    donor_id <- stringr::str_extract(sample_label, "(?i)(healthy|cirrhotic|blood)[0-9]+")
    donor_id <- ifelse(is.na(donor_id), sample_label, donor_id)
    donor_id <- stringr::str_replace_all(tolower(donor_id), "healthy", "H")
    donor_id <- stringr::str_replace_all(donor_id, "cirrhotic", "C")
    donor_id <- stringr::str_replace_all(donor_id, "blood", "B")

    fraction <- dplyr::case_when(
      grepl("cd45\\+", sample_label, ignore.case = TRUE) ~ "CD45pos",
      grepl("cd45-a", sample_label, ignore.case = TRUE) ~ "CD45neg_A",
      grepl("cd45-b", sample_label, ignore.case = TRUE) ~ "CD45neg_B",
      grepl("cd45-", sample_label, ignore.case = TRUE) ~ "CD45neg",
      TRUE ~ NA_character_
    )

    tibble::tibble(
      geo_accession = gsm,
      sample_id = prefix,
      sample_label = sample_label,
      donor_id = donor_id,
      disease_status = disease_status,
      fibrosis_stage_harmonized = dplyr::case_when(
        disease_status == "healthy" ~ "F0_or_nonfibrotic_control",
        disease_status == "cirrhotic" ~ "F4_cirrhosis",
        TRUE ~ NA_character_
      ),
      tissue = tissue,
      species = species,
      fraction = fraction,
      matrix_path = mtx,
      barcode_path = barcodes[1] %||% NA_character_,
      feature_path = genes[1] %||% NA_character_
    )
  })

  sample_table <- sample_table %>%
    mutate(
      include_primary = species == "human" & tissue == "liver" & disease_status %in% c("healthy", "cirrhotic"),
      disease_status = factor(disease_status, levels = c("healthy", "cirrhotic"))
    )

  if (!is.null(geo_metadata) && "geo_accession" %in% names(geo_metadata)) {
    sample_table <- sample_table %>%
      left_join(
        geo_metadata %>% dplyr::select(geo_accession, title, platform_id, everything()),
        by = "geo_accession",
        suffix = c("", "_geo")
      )
  }

  write_csv_safely(sample_table, file.path(config$paths$table_dir, "gse136103_sample_metadata.csv"))
  sample_table
}

summarize_primary_dataset <- function(sample_table, config = read_config()) {
  # GEOquery can introduce list-like metadata columns. Coerce only the fields used
  # in the compact dataset summary so dplyr/vctrs never tries to group on list cols.
  summary_input <- sample_table %>%
    mutate(
      species = as.character(species),
      tissue = as.character(tissue),
      disease_status = as.character(disease_status),
      fraction = as.character(fraction),
      include_primary = as.logical(include_primary),
      donor_id = as.character(donor_id)
    )

  summary <- summary_input %>%
    dplyr::count(species, tissue, disease_status, fraction, include_primary, name = "n_libraries") %>%
    dplyr::arrange(desc(include_primary), species, tissue, disease_status, fraction)

  donor_summary <- summary_input %>%
    dplyr::filter(include_primary) %>%
    dplyr::distinct(donor_id, disease_status) %>%
    dplyr::count(disease_status, name = "n_donors")

  write_csv_safely(summary, file.path(config$paths$table_dir, "gse136103_library_summary.csv"))
  write_csv_safely(donor_summary, file.path(config$paths$table_dir, "gse136103_donor_summary.csv"))

  # Return one flat tibble rather than a list, because some targets/caches choke
  # on list-valued return objects in older package combinations.
  dplyr::bind_rows(
    summary %>% dplyr::mutate(summary_type = "libraries") %>%
      dplyr::mutate(donor_id = NA_character_, n_donors = NA_integer_),
    donor_summary %>%
      dplyr::mutate(summary_type = "donors", species = NA_character_, tissue = NA_character_, fraction = NA_character_, include_primary = NA, n_libraries = NA_integer_) %>%
      dplyr::relocate(summary_type)
  )
}
