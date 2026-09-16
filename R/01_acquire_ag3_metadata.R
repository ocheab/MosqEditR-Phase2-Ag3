# STEP 01: prepare Ag3.0 metadata from the local targeted-acquisition cache
#
# IMPORTANT: this script performs NO network access.
# First run Step 00, then acquire the small local Ag3 handoff with:
#   Windows PowerShell: scripts/acquire_ag3_targeted.ps1
#   Git Bash/WSL/Linux: bash scripts/acquire_ag3_targeted.sh

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 01 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "01"
log_step(STEP, "Preparing Ag3.0 metadata from local targeted-acquisition files")
log_step(STEP, "Network access is disabled by design in this R step")

cfg <- read_config()
local_dir <- as.character(cfg_get(cfg, "local_acquisition", "output_dir", default = "data_raw/ag3_local"))
release <- as.character(cfg_get(cfg, "local_acquisition", "release", default = cfg_get(cfg, "ag3_release", default = "3.0")))
exclude_sets <- as.character(cfg_get(cfg, "exclude_sample_set", default = character()))
primary_taxa <- tolower(as.character(cfg_get(cfg, "primary_taxa", default = c("gambiae", "coluzzii"))))
primary_min_n <- as.integer(cfg_get(cfg, "primary_min_population_n", default = 20L))
sensitivity_ns <- as.integer(unlist(cfg_get(cfg, "sensitivity_min_population_n", default = c(10L, 30L))))

if (!length(primary_taxa)) fail_step(STEP, "No primary_taxa are configured.")
if (!is.finite(primary_min_n) || primary_min_n < 1L) fail_step(STEP, "primary_min_population_n must be a positive integer.")

meta_file <- file.path(local_dir, "metadata", paste0("ag3_", release, "_sample_metadata.csv"))
complete_file <- file.path(local_dir, "ACQUISITION_COMPLETE.ok")
manifest_file <- file.path(local_dir, "acquisition_manifest.json")
qc_file <- file.path(local_dir, "acquisition_qc.tsv")

if (!file.exists(meta_file)) {
  fail_step(STEP, paste0(
    "Local Ag3 metadata is missing: ", meta_file, "\n",
    "Run the targeted acquisition first. No 50-GB VCF is required.\n",
    "Windows PowerShell:\n  powershell -ExecutionPolicy Bypass -File scripts/acquire_ag3_targeted.ps1\n",
    "Git Bash/WSL/Linux/macOS:\n  bash scripts/acquire_ag3_targeted.sh"
  ))
}
if (!file.exists(complete_file)) {
  fail_step(STEP, paste0(
    "Local acquisition is incomplete (missing ", complete_file, "). ",
    "Rerun the targeted acquisition; validated completed contigs are cached and reused."
  ))
}

meta <- data.table::fread(meta_file, na.strings = c("", "NA", "NaN", "nan", "None"))
if (!nrow(meta)) fail_step(STEP, "Local Ag3 metadata file contains zero rows.")
assert_columns(meta, "sample_id", "local Ag3 metadata")
assert_unique(meta, "sample_id", "local Ag3 metadata")
meta[, sample_id := trimws(as.character(sample_id))]
if (any(!nonempty_string(meta$sample_id))) fail_step(STEP, "Blank sample_id detected in local Ag3 metadata.")

# A sample_set column is expected from malariagen_data.sample_metadata().
sample_set_col <- resolve_column(meta, c("sample_set", "sample_set_id"), "sample-set column", required = FALSE)
if (is.na(sample_set_col)) {
  fail_step(STEP, paste0(
    "The local metadata does not contain a sample-set column. Available columns: ",
    paste(names(meta), collapse = ", ")
  ))
}
if (sample_set_col != "sample_set") data.table::setnames(meta, sample_set_col, "sample_set")
meta[, sample_set := trimws(as.character(sample_set))]
if (length(exclude_sets)) meta <- meta[!sample_set %in% exclude_sets]
if (!nrow(meta)) fail_step(STEP, "No metadata rows remain after configured sample-set exclusions.")
assert_unique(meta, "sample_id", "Ag3 metadata after sample-set exclusions")

# Current malariagen_data metadata normally exposes 'taxon'. Keep fallbacks for
# schema drift and older locally cached tables.
taxon_col <- resolve_column(
  meta,
  c("taxon", "aim_species", "species_gambcolu_arabiensis", "species"),
  "species/taxon column",
  required = FALSE
)
if (is.na(taxon_col)) {
  fail_step(STEP, paste0(
    "No recognizable species/taxon column was found in the local Ag3 metadata. Available columns: ",
    paste(names(meta), collapse = ", ")
  ))
}
meta[, taxon_raw := as.character(get(taxon_col))]
meta[, taxon := normalize_taxon(taxon_raw)]

# Resolve location/time columns without assuming one immutable API schema.
country_col <- resolve_column(meta, c("country", "country_name"), "country column", required = FALSE)
location_col <- resolve_column(meta, c("location", "location_name", "admin2_name", "admin1_name"), "location column", required = FALSE)
year_col <- resolve_column(meta, c("year", "collection_year", "year_collected"), "collection year column", required = FALSE)

meta[, country := if (!is.na(country_col)) trimws(as.character(get(country_col))) else NA_character_]
meta[!nonempty_string(country), country := NA_character_]
meta[, location := if (!is.na(location_col)) trimws(as.character(get(location_col))) else NA_character_]
meta[!nonempty_string(location), location := NA_character_]
meta[, year_original := if (!is.na(year_col)) as.character(get(year_col)) else NA_character_]
meta[, year_numeric := suppressWarnings(as.numeric(year_original))]
meta[, primary_taxon := taxon %in% primary_taxa]
meta[is.na(primary_taxon), primary_taxon := FALSE]

# Population IDs are created only where both country and taxon are known.
meta[, population_id := NA_character_]
meta[primary_taxon == TRUE & nonempty_string(country), population_id := paste(country, taxon, sep = "__")]
meta[, country_taxon_population_defined := primary_taxon == TRUE & nonempty_string(population_id)]

# Step 01 counts are metadata-based. Step 02 recalculates them after intersection
# with the locally acquired phased-panel sample axis.
pop_counts <- meta[country_taxon_population_defined == TRUE, .(
  n_samples = data.table::uniqueN(sample_id),
  n_locations = data.table::uniqueN(location[nonempty_string(location)]),
  year_min = safe_min(year_numeric),
  year_max = safe_max(year_numeric)
), by = .(population_id, country, taxon)]
pop_counts[, primary_population_n20 := n_samples >= primary_min_n]
for (nmin in unique(sensitivity_ns)) {
  pop_counts[, (paste0("sensitivity_population_n", nmin)) := n_samples >= nmin]
}
for (nmin in c(10L, 30L)) {
  nm <- paste0("sensitivity_population_n", nmin)
  if (!nm %in% names(pop_counts)) pop_counts[, (nm) := n_samples >= nmin]
}

# Per-sample-set QC from the single local metadata table.
set_qc <- meta[, .(
  n_metadata = .N,
  n_taxon_nonmissing = sum(nonempty_string(taxon)),
  n_primary_taxon = sum(primary_taxon == TRUE, na.rm = TRUE),
  n_primary_with_population = sum(country_taxon_population_defined == TRUE, na.rm = TRUE)
), by = sample_set]

# Deterministic ordering for reproducibility.
data.table::setorder(meta, sample_set, sample_id)
data.table::setorder(pop_counts, country, taxon, population_id)
data.table::setorder(set_qc, sample_set)

atomic_fwrite(meta, "data_processed/01_ag3_sample_metadata.csv")
atomic_fwrite(pop_counts, "data_processed/01_ag3_population_manifest.csv")
atomic_fwrite(set_qc, "data_processed/01_ag3_sample_set_qc.csv")

qc <- data.table::data.table(
  metric = c(
    "Non-cross sample sets", "Samples in metadata", "Primary gambiae/coluzzii samples",
    "Primary samples with country-defined population", "Primary country-taxon populations n>=20",
    "Sensitivity populations n>=10", "Sensitivity populations n>=30",
    "Primary samples without country-defined population", "R network requests in Step 01"
  ),
  value = c(
    data.table::uniqueN(meta$sample_set), data.table::uniqueN(meta$sample_id), data.table::uniqueN(meta[primary_taxon == TRUE, sample_id]),
    data.table::uniqueN(meta[country_taxon_population_defined == TRUE, sample_id]),
    sum(pop_counts$primary_population_n20),
    sum(pop_counts$sensitivity_population_n10),
    sum(pop_counts$sensitivity_population_n30),
    data.table::uniqueN(meta[primary_taxon == TRUE & !country_taxon_population_defined, sample_id]),
    0L
  )
)
atomic_fwrite(qc, "data_processed/01_ag3_metadata_qc.csv")

checksum_inputs <- c(
  meta_file,
  if (file.exists(manifest_file)) manifest_file else character(),
  if (file.exists(qc_file)) qc_file else character(),
  "data_processed/01_ag3_sample_metadata.csv",
  "data_processed/01_ag3_population_manifest.csv",
  "data_processed/01_ag3_sample_set_qc.csv"
)
write_checksum(checksum_inputs, "logs/01_checksums.tsv")
write_session_info("logs/01_sessionInfo.txt")
print(qc)
log_step(STEP, "Local Ag3.0 metadata preparation completed successfully")
