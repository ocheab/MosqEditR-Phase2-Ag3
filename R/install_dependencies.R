# Install/check R dependencies for MosqEditR Manuscript 2 hybrid workflow.
# Windows-safe revision: all remote VCF/tabix access is performed by
# Rsamtools/VariantAnnotation; native-Windows pysam is not required.

options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout", 60)))

cran <- c("data.table", "ggplot2", "scales", "zip", "yaml", "curl", "jsonlite", "digest")
bioc <- c("GenomicRanges", "IRanges", "Rsamtools", "VariantAnnotation", "S4Vectors")

repos <- getOption("repos")
if (is.null(repos) || !length(repos) || identical(unname(repos[[1]]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_cran_if_missing <- function(pkgs) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) {
    message("Installing CRAN package(s): ", paste(need, collapse = ", "))
    install.packages(need, dependencies = TRUE)
  }
  still <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still)) stop("CRAN package installation failed for: ", paste(still, collapse = ", "), call. = FALSE)
}

install_cran_if_missing(cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
need_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_bioc)) {
  message("Installing Bioconductor package(s): ", paste(need_bioc, collapse = ", "))
  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
}
still_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_bioc)) stop("Bioconductor package installation failed for: ", paste(still_bioc, collapse = ", "), call. = FALSE)

all_pkgs <- c(cran, bioc, "BiocManager")
ver_tab <- data.frame(
  package = all_pkgs,
  version = vapply(all_pkgs, function(p) as.character(utils::packageVersion(p)), character(1)),
  stringsAsFactors = FALSE
)
print(ver_tab, row.names = FALSE)
cat("\nDependencies ready.\n")
cat("Remote VCF engine: R/Bioconductor Rsamtools + VariantAnnotation.\n")
cat("Python helper is used only for Phase-2 metadata/accessibility HDF5; pysam is not used.\n")
