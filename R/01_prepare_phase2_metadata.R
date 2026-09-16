# STEP 01: prepare Ag1000G Phase 2 AR1 discovery metadata from the local cache
#
# This step is deliberately OFFLINE. The public Phase-2 acquisition layer is:
#   scripts/acquire_phase2_discovery.py
# R never contacts Sanger/GCS in this step.

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
suppressPackageStartupMessages(library(data.table))

STEP <- "01"
log_step(STEP, "Preparing Ag1000G Phase 2 AR1 phased-discovery metadata from local files")
log_step(STEP, "Network access is disabled by design in this R step")

cfg <- read_config()
local_dir <- as.character(cfg_get(cfg, "phase2_acquisition", "output_dir", default = "data_raw/phase2_local"))
primary_taxa <- tolower(as.character(cfg_get(cfg, "primary_taxa", default = c("gambiae", "coluzzii"))))
primary_min_n <- as.integer(cfg_get(cfg, "primary_min_population_n", default = 20L))
sensitivity_ns <- unique(c(10L, 30L, as.integer(unlist(cfg_get(cfg, "sensitivity_min_population_n", default = c(10L, 30L))))))

meta_file <- file.path(local_dir, "metadata", "phased_sample_metadata.csv")
all_meta_file <- file.path(local_dir, "metadata", "sample_metadata.csv")
complete_file <- file.path(local_dir, "ACQUISITION_COMPLETE.ok")
manifest_file <- file.path(local_dir, "acquisition_manifest.json")

assert_files(c(meta_file, complete_file))
meta <- fread(meta_file, na.strings = c("", "NA", "NaN", "nan", "None"))
if (!nrow(meta)) fail_step(STEP, "Phase 2 phased-sample metadata contains zero rows.")
assert_columns(meta, c("sample_id", "taxon", "primary_taxon"), "Phase 2 phased-sample metadata")
assert_unique(meta, "sample_id", "Phase 2 phased-sample metadata")
meta[, sample_id := trimws(as.character(sample_id))]
if (any(!nonempty_string(meta$sample_id))) fail_step(STEP, "Blank sample_id detected in Phase 2 metadata.")

# Normalize columns produced by the Python handoff while retaining original
# release fields for auditability.
meta[, taxon := normalize_taxon(as.character(taxon))]
meta[, primary_taxon := normalize_flag(primary_taxon)]
meta[is.na(primary_taxon), primary_taxon := FALSE]
meta[, primary_taxon := primary_taxon & taxon %in% primary_taxa]

if (!"sample_set" %in% names(meta)) meta[, sample_set := "AG1000G-PHASE2-AR1"]
meta[, sample_set := trimws(as.character(sample_set))]

if (!"country" %in% names(meta)) meta[, country := NA_character_]
if (!"location" %in% names(meta)) {
  location_col <- resolve_column(meta, c("region", "site", "location_name"), "Phase 2 location", required = FALSE)
  meta[, location := if (!is.na(location_col)) as.character(get(location_col)) else NA_character_]
}
meta[, country := trimws(as.character(country))]
meta[!nonempty_string(country), country := NA_character_]
meta[, location := trimws(as.character(location))]
meta[!nonempty_string(location), location := NA_character_]

if (!"year" %in% names(meta)) meta[, year := NA_real_]
meta[, year_original := as.character(year)]
meta[, year_numeric := suppressWarnings(as.numeric(year_original))]

if (!"phase2_population" %in% names(meta)) {
  pop_col <- resolve_column(meta, c("population", "population_id"), "Phase 2 release population", required = FALSE)
  meta[, phase2_population := if (!is.na(pop_col)) trimws(as.character(get(pop_col))) else NA_character_]
} else {
  meta[, phase2_population := trimws(as.character(phase2_population))]
}
meta[!nonempty_string(phase2_population), phase2_population := NA_character_]

# AUTHORITATIVE PHASE-2 POPULATION/TAXON DEFINITION
# -------------------------------------------------
# The Phase-2 release population code is the analysis stratum.  Population codes
# ending in "gam" are explicitly An. gambiae and those ending in "col" are
# explicitly An. coluzzii.  GM, GW and KE are mixed/undetermined Phase-2 groups
# and are therefore not part of the primary gambiae/coluzzii discovery cohort.
#
# Preserve the historical M/S-marker-based handoff fields for auditability, but
# do not use them to define the primary discovery taxon.
meta[, taxon_marker_based_handoff := taxon]
meta[, primary_taxon_marker_based_handoff := primary_taxon]

pop_upper <- toupper(trimws(as.character(meta$phase2_population)))
meta[, taxon := data.table::fcase(
  grepl("GAM$", pop_upper), "gambiae",
  grepl("COL$", pop_upper), "coluzzii",
  default = NA_character_
)]
meta[, primary_taxon := !is.na(taxon) & taxon %in% primary_taxa]

meta[, population_id := phase2_population]
meta[, population_id_source := "PHASE2_RELEASE_POPULATION"]
meta[, population_defined := primary_taxon == TRUE & nonempty_string(population_id)]

# The primary population code must map to one and only one country and taxon.
descriptor_qc <- meta[population_defined == TRUE, .(
  n_countries = data.table::uniqueN(country[nonempty_string(country)]),
  n_taxa = data.table::uniqueN(taxon[nonempty_string(taxon)])
), by = population_id]
if (any(descriptor_qc$n_countries > 1L | descriptor_qc$n_taxa != 1L)) {
  bad <- descriptor_qc[n_countries > 1L | n_taxa != 1L, population_id]
  fail_step(
    STEP,
    paste0(
      "Published Phase-2 population IDs do not map uniquely to country/taxon for: ",
      paste(utils::head(bad, 20L), collapse = ", ")
    )
  )
}

# Sex is used only for the X-chromosome haploid fallback in Step 02.
sex_col <- resolve_column(meta, c("sex", "sex_call", "gender"), "sex column", required = FALSE)
meta[, sex_normalized := if (!is.na(sex_col)) tolower(trimws(as.character(get(sex_col)))) else NA_character_]
meta[, known_male_for_x_fallback := sex_normalized %in% c("m", "male", "1", "xy")]
meta[is.na(known_male_for_x_fallback), known_male_for_x_fallback := FALSE]

if (!any(meta$primary_taxon)) fail_step(STEP, "No gambiae/coluzzii Phase 2 phased samples remain after taxon normalization.")
marker_release_disagreement <- meta[
  nonempty_string(taxon_marker_based_handoff) &
    nonempty_string(taxon) &
    normalize_taxon(as.character(taxon_marker_based_handoff)) != taxon,
  .N
]
mixed_or_undetermined_n <- meta[!primary_taxon & nonempty_string(phase2_population), .N]
log_step(
  STEP,
  paste0(
    "Release-population taxonomy retained ", sum(meta$primary_taxon), " primary gambiae/coluzzii samples; ",
    mixed_or_undetermined_n, " phased samples are mixed/undetermined or outside the primary release-coded taxa; ",
    marker_release_disagreement, " sample(s) have an M/S-marker taxon differing from the release-population assignment."
  )
)

if (any(meta$primary_taxon & !meta$population_defined)) {
  log_warning(STEP, paste0(
    sum(meta$primary_taxon & !meta$population_defined),
    " primary Phase 2 sample(s) lack a usable population ID and will be excluded from population-level metrics."
  ))
}

pop_counts <- meta[population_defined == TRUE, .(
  country = {
    z <- sort(unique(country[nonempty_string(country)]))
    if (length(z)) z[[1]] else NA_character_
  },
  taxon = {
    z <- sort(unique(taxon[nonempty_string(taxon)]))
    if (length(z)) z[[1]] else NA_character_
  },
  n_samples = uniqueN(sample_id),
  n_locations = uniqueN(location[nonempty_string(location)]),
  year_min = safe_min(year_numeric),
  year_max = safe_max(year_numeric),
  population_id_source = paste(sort(unique(population_id_source)), collapse = ";")
), by = population_id]
assert_unique(pop_counts, "population_id", "Phase-2 discovery population manifest")
pop_counts[, primary_population_n20 := n_samples >= primary_min_n]
for (nmin in sensitivity_ns) pop_counts[, (paste0("sensitivity_population_n", nmin)) := n_samples >= nmin]
for (nmin in c(10L, 30L)) {
  nm <- paste0("sensitivity_population_n", nmin)
  if (!nm %in% names(pop_counts)) pop_counts[, (nm) := n_samples >= nmin]
}

set_qc <- meta[, .(
  n_metadata = .N,
  n_taxon_nonmissing = sum(nonempty_string(taxon)),
  n_primary_taxon = sum(primary_taxon == TRUE, na.rm = TRUE),
  n_primary_with_population = sum(population_defined == TRUE, na.rm = TRUE),
  n_known_males = sum(known_male_for_x_fallback == TRUE, na.rm = TRUE)
), by = sample_set]

data.table::setorder(meta, population_id, sample_id, na.last = TRUE)
data.table::setorder(pop_counts, country, taxon, population_id)
data.table::setorder(set_qc, sample_set)

atomic_fwrite(meta, "data_processed/01_discovery_sample_metadata.csv")
atomic_fwrite(pop_counts, "data_processed/01_discovery_population_manifest.csv")
atomic_fwrite(set_qc, "data_processed/01_discovery_sample_set_qc.csv")

qc <- data.table(
  metric = c(
    "Phase 2 phased samples in local handoff",
    "Primary gambiae/coluzzii phased samples",
    "Primary samples with population ID",
    "Primary populations n>=20",
    "Sensitivity populations n>=10",
    "Sensitivity populations n>=30",
    "Primary samples using country-taxon fallback population",
    "R network requests in Step 01"
  ),
  value = c(
    uniqueN(meta$sample_id),
    uniqueN(meta[primary_taxon == TRUE, sample_id]),
    uniqueN(meta[population_defined == TRUE, sample_id]),
    sum(pop_counts$primary_population_n20),
    sum(pop_counts$sensitivity_population_n10),
    sum(pop_counts$sensitivity_population_n30),
    uniqueN(meta[primary_taxon == TRUE & population_id_source == "COUNTRY_TAXON_FALLBACK", sample_id]),
    0L
  )
)
atomic_fwrite(qc, "data_processed/01_discovery_metadata_qc.csv")

checksum_inputs <- c(
  meta_file,
  if (file.exists(all_meta_file)) all_meta_file else character(),
  if (file.exists(manifest_file)) manifest_file else character(),
  "data_processed/01_discovery_sample_metadata.csv",
  "data_processed/01_discovery_population_manifest.csv",
  "data_processed/01_discovery_sample_set_qc.csv"
)
write_checksum(checksum_inputs, "logs/01_checksums.tsv")
write_session_info("logs/01_sessionInfo.txt")
print(qc)
log_step(STEP, "Phase 2 discovery metadata preparation completed successfully")
