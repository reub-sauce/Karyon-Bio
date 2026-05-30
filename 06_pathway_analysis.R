# Pathway and mechanism analysis using fgsea/MSigDB plus curated fibrosis mechanism tags.

get_msigdb_sets <- function() {
  fml <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% fml) {
    hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    reactome <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
    go_bp <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  } else {
    hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    reactome <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
    go_bp <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  }
  dplyr::bind_rows(sanitize_table(hallmark), sanitize_table(reactome), sanitize_table(go_bp)) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    mutate(gene_symbol = clean_gene_symbols(gene_symbol)) %>%
    distinct()
}

run_fgsea_one_group <- function(de_tbl, group_name, gene_sets, config = read_config()) {
  tbl <- sanitize_table(de_tbl) %>% filter(analysis_group == group_name) %>% filter(!is.na(PValue), !is.na(logFC))
  if (nrow(tbl) < 50) return(tibble::tibble())
  stats <- tbl %>%
    mutate(rank_stat = sign(logFC) * safe_neg_log10(PValue)) %>%
    arrange(desc(abs(rank_stat))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    dplyr::select(gene, rank_stat)
  ranks <- stats$rank_stat
  names(ranks) <- clean_gene_symbols(stats$gene)
  ranks <- sort(ranks, decreasing = TRUE)
  pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  fg <- suppressWarnings(
    fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = ranks,
      minSize = config$pathway$min_size,
      maxSize = config$pathway$max_size
    )
  )
  sanitize_table(fg) %>%
    mutate(
      analysis_group = group_name,
      leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
    ) %>%
    arrange(padj)
}

curated_mechanism_hits <- function(de_tbl) {
  sets <- yaml::read_yaml("resources/compartment_gene_sets.yml")$fibrosis_mechanism_sets
  tibble::tibble(
    mechanism = rep(names(sets), lengths(sets)),
    gene = clean_gene_symbols(unlist(sets, use.names = FALSE))
  ) %>%
    inner_join(de_tbl %>% mutate(gene = clean_gene_symbols(gene)), by = "gene") %>%
    mutate(mechanism_hit = TRUE)
}

plot_top_pathways <- function(fgsea_tbl, config = read_config()) {
  if (nrow(fgsea_tbl) == 0) return(invisible(NULL))
  req_groups <- unique(required_compartment_filter(fgsea_tbl %>% rename(logFC = NES, FDR = padj, gene = pathway))$analysis_group)
  tbl <- fgsea_tbl %>%
    filter(analysis_group %in% req_groups | stringr::str_detect(analysis_group, "HSC|Macrophage|Endothelial|Scar|Myofibroblast|ACKR1|PLVAP")) %>%
    group_by(analysis_group) %>%
    slice_min(order_by = padj, n = 8, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(pathway_short = stringr::str_trunc(pathway, 60), signed_score = NES * safe_neg_log10(padj))

  if (nrow(tbl) == 0) return(invisible(NULL))
  p <- ggplot(tbl, aes(x = signed_score, y = reorder(pathway_short, signed_score))) +
    geom_col() +
    facet_wrap(~ analysis_group, scales = "free_y") +
    theme_bw(base_size = 9) +
    labs(x = "NES * -log10(FDR)", y = NULL, title = "Top enriched pathways in required disease-relevant compartments")
  save_plot_safely(p, file.path(config$paths$figure_dir, "pathway_top_terms_by_compartment.pdf"), width = 13, height = 9)
}

run_pathway_analysis <- function(de_result, config = read_config()) {
  message_header("Running pathway and mechanism analysis")
  de_tbl <- sanitize_table(de_result$all_de)
  gene_sets <- get_msigdb_sets()
  groups <- unique(de_tbl$analysis_group)
  fg <- purrr::map_dfr(groups, ~ run_fgsea_one_group(de_tbl, .x, gene_sets, config)) %>% sanitize_table()
  if (nrow(fg) > 0) {
    write_csv_safely(fg, file.path(config$paths$table_dir, "fgsea_all_results.csv"))
  } else {
    write_csv_safely(tibble::tibble(), file.path(config$paths$table_dir, "fgsea_all_results.csv"))
  }

  mech <- curated_mechanism_hits(de_tbl) %>% sanitize_table()
  write_csv_safely(mech, file.path(config$paths$table_dir, "curated_fibrosis_mechanism_hits.csv"))
  plot_top_pathways(fg, config)
  list(fgsea = fg, curated_mechanism_hits = mech)
}
