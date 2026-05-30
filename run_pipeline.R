message("Starting Karyon Bio liver fibrosis workflow...")
if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Package 'targets' is required. Run: source('scripts/install_packages.R')", call. = FALSE)
}
targets::tar_make()
message("Workflow complete. See results/figures, results/tables, and results/reports.")
