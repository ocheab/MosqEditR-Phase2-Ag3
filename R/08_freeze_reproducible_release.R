# STEP 08: freeze a clean, auditable reproducible release and ZIP archive

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 08 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "zip"))

STEP <- "08"
log_step(STEP, "Freezing reproducible Manuscript 2 release")

# These outputs establish that every analytical stage completed before a release is frozen.
required_outputs <- c(
  "data_processed/00_frozen_target_sites.csv",
  "data_processed/00_frozen_gene_panel.csv",
  "data_processed/01_discovery_sample_metadata.csv",
  "data_processed/02_target_haplotype_states.rds",
  "data_processed/02_analysis_population_manifest.csv",
  "data_processed/03_population_site_metrics.csv",
  "data_processed/03_site_population_metrics.csv",
  "data_processed/04_multiplex_portfolio_metrics.csv",
  "data_processed/05_gene_population_robustness.csv",
  "data_processed/06_rank_validation_summary.csv",
  "data_processed/05b_external_validation_targets.csv",
  "data_processed/06a_ag3_external_site_validation.csv",
  "data_processed/06a_ag3_external_gene_support.csv",
  "tables/Manuscript2_results_snapshot.csv",
  "manuscript/MANUSCRIPT_RESULT_VALUES.txt",
  "tables/Figure_output_manifest.csv"
)
assert_files(required_outputs, nonempty = TRUE)

# Final structural validation before copying.
site <- data.table::fread("data_processed/03_population_site_metrics.csv", na.strings = c("", "NA", "NaN"))
gene <- data.table::fread("data_processed/05_gene_population_robustness.csv", na.strings = c("", "NA", "NaN"))
assert_unique(site, "population_site_id", "final site metrics")
assert_unique(gene, "gene_id", "final gene metrics")
if (!nrow(site) || !nrow(gene)) fail_step(STEP, "Final analysis tables are unexpectedly empty.")

release_root <- "release"
out_name <- "MosqEditR_Manuscript2_PHASE2_DISCOVERY_AG3_EXTERNAL"
out <- file.path(release_root, out_name)
zipfile <- file.path(release_root, paste0(out_name, ".zip"))

# Delete an old release directory first so stale files cannot survive a rerun.
if (dir.exists(out)) unlink(out, recursive = TRUE, force = TRUE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(zipfile)) unlink(zipfile, force = TRUE)

copy_tree <- function(src_dir, dest_root) {
  if (!dir.exists(src_dir)) return(invisible(NULL))
  files <- list.files(src_dir, full.names = TRUE, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(invisible(NULL))
  info <- file.info(files)
  files <- files[!info$isdir]
  if (!length(files)) return(invisible(NULL))
  rel <- substring(files, nchar(src_dir) + 2L)
  dest <- file.path(dest_root, src_dir, rel)
  invisible(lapply(unique(dirname(dest)), dir.create, recursive = TRUE, showWarnings = FALSE))
  ok <- file.copy(files, dest, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  if (!all(ok)) stop("Failed to copy one or more files from ", src_dir, " into release.", call. = FALSE)
  invisible(dest)
}

# data_raw is intentionally excluded to prevent accidental inclusion of remote caches or genomic chunks.
# Reproducible acquisition code and immutable manifests are included; compact analysis-ready
# products are frozen under data_processed/.
for (d in c("R", "scripts", "data_input", "data_processed", "figures", "tables", "supplement", "manuscript", "metadata", "logs", "optional")) {
  copy_tree(d, out)
}
for (f in c("analysis_config.yml", "README.md", "RUN_INSTRUCTIONS.md", "run_all.R", "ROBUSTNESS_CHANGELOG.md", "FILE_MANIFEST.csv", "SHA256SUMS.txt")) {
  if (file.exists(f)) {
    ok <- file.copy(f, file.path(out, basename(f)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    if (!ok) stop("Failed to copy release file: ", f, call. = FALSE)
  }
}

# Add a release-specific session record and manifest.
atomic_writeLines(capture.output(utils::sessionInfo()), file.path(out, "RELEASE_SESSION_INFO.txt"))
manifest_files <- list.files(out, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
manifest_files <- manifest_files[!file.info(manifest_files)$isdir]
manifest <- data.table::data.table(
  file = substring(manifest_files, nchar(out) + 2L),
  bytes = as.numeric(file.info(manifest_files)$size),
  md5 = unname(tools::md5sum(manifest_files))
)
data.table::setorder(manifest, file)
data.table::fwrite(manifest, file.path(out, "MANIFEST_MD5.tsv"), sep = "\t")

# Human-readable release QC.
release_qc <- data.table::data.table(
  metric = c(
    "Release target sites", "Release genes", "Release files before ZIP",
    "Finite site robustness scores", "Finite gene robustness scores"
  ),
  value = c(
    nrow(site), nrow(gene), nrow(manifest),
    sum(is.finite(site$population_target_site_robustness_score)),
    sum(is.finite(gene$gene_population_robustness_score))
  )
)
data.table::fwrite(release_qc, file.path(out, "RELEASE_QC.tsv"), sep = "\t")

# zip::zipr is cross-platform and does not depend on an external zip executable.
zip::zipr(
  zipfile = zipfile,
  files = out_name,
  root = release_root,
  recurse = TRUE,
  include_directories = TRUE
)
if (!file.exists(zipfile) || file.info(zipfile)$size <= 0) fail_step(STEP, "Release ZIP was not created correctly.")

zip_md5 <- unname(tools::md5sum(zipfile))
atomic_writeLines(
  c(
    paste0("Release ZIP: ", zipfile),
    paste0("Bytes: ", file.info(zipfile)$size),
    paste0("MD5: ", zip_md5),
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  file.path(release_root, paste0(out_name, "_ZIP_CHECKSUM.txt"))
)

print(release_qc)
cat("\nRelease directory: ", normalizePath(out, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Release ZIP: ", normalizePath(zipfile, winslash = "/", mustWork = TRUE), "\n", sep = "")
log_step(STEP, "Release freeze completed successfully")
