# LOCAL-FIRST Ag3 acquisition diagnostic.
# This diagnostic performs no network access from R.

if (!file.exists("R/helpers.R")) stop("Run this diagnostic from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "01a"
log_step(STEP, "Checking local-first Ag3 acquisition prerequisites")
cfg <- read_config()
local_dir <- as.character(cfg_get(cfg, "local_acquisition", "output_dir", default = "data_raw/ag3_local"))
release <- as.character(cfg_get(cfg, "local_acquisition", "release", default = "3.0"))

pointer_file <- file.path("metadata", "ag3_python_path.txt")
pointer_python <- if (file.exists(pointer_file)) trimws(readLines(pointer_file, warn = FALSE, n = 1L)) else character()
explicit_venv <- Sys.getenv("MOSQEDIT_AG3_VENV", unset = "")
local_appdata <- Sys.getenv("LOCALAPPDATA", unset = "")

candidates <- unique(Filter(nzchar, c(
  if (nzchar(explicit_venv)) file.path(explicit_venv, "Scripts", "python.exe") else "",
  pointer_python,
  if (nzchar(local_appdata)) file.path(local_appdata, "MosqEditR", "ag3_venv", "Scripts", "python.exe") else "",
  file.path(".venv_ag3", "Scripts", "python.exe"),
  file.path(".venv_ag3", "bin", "python")
)))
py <- candidates[file.exists(candidates)]
python_ready <- length(py) > 0L
python_version <- NA_character_
if (python_ready) {
  python_version <- tryCatch(
    paste(system2(py[[1]], "--version", stdout = TRUE, stderr = TRUE), collapse = " "),
    error = function(e) paste0("ERROR: ", conditionMessage(e))
  )
}

meta_file <- file.path(local_dir, "metadata", paste0("ag3_", release, "_sample_metadata.csv"))
complete_file <- file.path(local_dir, "ACQUISITION_COMPLETE.ok")
qc_file <- file.path(local_dir, "acquisition_qc.tsv")
access_file <- file.path(local_dir, "accessibility", "target_accessibility.tsv")

status <- data.table::data.table(
  check = c(
    "Frozen Step-00 targets", "Project Python environment", "Python version",
    "Local Ag3 metadata", "Local target accessibility", "Acquisition completion marker"
  ),
  status = c(
    if (file.exists("data_processed/00_frozen_target_sites.csv")) "PASS" else "MISSING",
    if (python_ready) "PASS" else "MISSING",
    if (python_ready) python_version else "N/A",
    if (file.exists(meta_file)) "PASS" else "MISSING",
    if (file.exists(access_file)) "PASS" else "MISSING",
    if (file.exists(complete_file)) "PASS" else "MISSING"
  )
)
print(status)

if (file.exists(qc_file)) {
  cat("\nAcquisition QC:\n")
  print(data.table::fread(qc_file, sep = "\t"))
}

cat("\nThis architecture never requires downloading a ~50-GB chromosome VCF.\n")
if (!python_ready) {
  cat("Windows: powershell -ExecutionPolicy Bypass -File scripts/setup_ag3_python.ps1 (uses a short environment path under LOCALAPPDATA by default)\n")
  cat("Bash/WSL/Linux/macOS: bash scripts/setup_ag3_python.sh\n")
}
if (!file.exists(complete_file)) {
  cat("After setup, acquire only target slices:\n")
  cat("Windows: powershell -ExecutionPolicy Bypass -File scripts/acquire_ag3_targeted.ps1\n")
  cat("Bash/WSL/Linux/macOS: bash scripts/acquire_ag3_targeted.sh\n")
}
