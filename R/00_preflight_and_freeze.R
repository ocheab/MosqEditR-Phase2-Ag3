# STEP 00: preflight, frozen-input validation, and immutable analysis copies

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 00 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "00"
log_step(STEP, "Starting frozen-panel preflight and integrity checks")

cfg <- read_config()
canonical <- as.character(cfg_get(cfg, "canonical_contigs", required = TRUE))
expected <- cfg_get(cfg, "expected_frozen_counts", default = list())

required <- c(
  "analysis_config.yml",
  "data_input/13_population_site_manifest.csv",
  "data_input/13_population_validation_panel.csv",
  "data_input/12_bootstrap_rank_uncertainty.csv"
)
assert_files(required)

sites <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
panel <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))
boot <- data.table::fread(required[[4]], na.strings = c("", "NA", "NaN"))

site_required <- c(
  "population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end",
  "guide_genomic_strand", "protospacer_20nt", "pam", "population_site_rank_within_gene",
  "final_prepopulation_rank", "population_validation_role"
)
panel_required <- c(
  "gene_id", "final_prepopulation_rank", "is_external_benchmark", "population_validation_role"
)
boot_required <- c("gene_id", "bootstrap_rank_median", "bootstrap_rank_q025", "bootstrap_rank_q975")
assert_columns(sites, site_required, "site manifest")
assert_columns(panel, panel_required, "gene panel")
assert_columns(boot, boot_required, "bootstrap rank table")
assert_no_missing(sites, c("population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end", "guide_genomic_strand", "protospacer_20nt", "pam"), "site manifest")
assert_no_missing(panel, c("gene_id", "final_prepopulation_rank"), "gene panel")
assert_unique(sites, "population_site_id", "site manifest")
assert_unique(panel, "gene_id", "gene panel")
assert_unique(boot, "gene_id", "bootstrap rank table")

# Normalize booleans imported from CSV and fail if values are not interpretable.
panel[, is_external_benchmark := normalize_flag(is_external_benchmark)]
if (anyNA(panel$is_external_benchmark)) {
  fail_step(STEP, "is_external_benchmark contains values that cannot be interpreted as TRUE/FALSE.")
}

# Referential integrity.
missing_genes <- setdiff(unique(sites$gene_id), panel$gene_id)
if (length(missing_genes)) {
  fail_step(STEP, paste0("Site manifest contains gene(s) absent from gene panel: ", paste(utils::head(missing_genes, 20L), collapse = ", ")))
}
missing_boot <- setdiff(panel$gene_id, boot$gene_id)
if (length(missing_boot)) {
  fail_step(STEP, paste0("Bootstrap table is missing panel gene(s): ", paste(utils::head(missing_boot, 20L), collapse = ", ")))
}

# Genomic/site integrity.
sites[, genomic_start := as.integer(genomic_start)]
sites[, genomic_end := as.integer(genomic_end)]
if (anyNA(sites$genomic_start) || anyNA(sites$genomic_end)) fail_step(STEP, "Genomic coordinates could not be parsed as integers.")
if (any(sites$genomic_start < 1L | sites$genomic_end < sites$genomic_start)) fail_step(STEP, "Invalid genomic interval(s) detected.")
if (any((sites$genomic_end - sites$genomic_start + 1L) != 23L)) {
  bad <- sites[(genomic_end - genomic_start + 1L) != 23L, population_site_id]
  fail_step(STEP, paste0("Every frozen target must span exactly 23 bp. Invalid site(s): ", paste(utils::head(bad, 20L), collapse = ", ")))
}
if (any(!sites$guide_genomic_strand %in% c("+", "-"))) fail_step(STEP, "guide_genomic_strand must contain only '+' or '-'.")
if (any(nchar(as.character(sites$protospacer_20nt)) != 20L, na.rm = TRUE)) fail_step(STEP, "protospacer_20nt contains sequence(s) not exactly 20 nt long.")
if (any(nchar(as.character(sites$pam)) != 3L, na.rm = TRUE)) fail_step(STEP, "pam contains sequence(s) not exactly 3 nt long.")
if (any(!grepl("^[ACGTNacgtn]{20}$", as.character(sites$protospacer_20nt)))) fail_step(STEP, "protospacer_20nt contains non-DNA characters.")
if (any(!grepl("^[ACGTNacgtn]{3}$", as.character(sites$pam)))) fail_step(STEP, "pam contains non-DNA characters.")

# Each gene should map to one contig and one frozen pre-population rank.
gene_structure <- sites[, .(
  n_sites = .N,
  n_contigs = data.table::uniqueN(genomic_seqid),
  n_prepopulation_ranks = data.table::uniqueN(final_prepopulation_rank)
), by = gene_id]
if (any(gene_structure$n_contigs != 1L)) {
  fail_step(STEP, "At least one gene has frozen target sites on more than one contig; portfolio span calculations would be invalid.")
}
if (any(gene_structure$n_prepopulation_ranks != 1L)) {
  fail_step(STEP, "At least one gene has inconsistent final_prepopulation_rank values across its sites.")
}

expected_sites_per_gene <- as.integer(expected$sites_per_gene %||% NA_integer_)
if (is.finite(expected_sites_per_gene) && any(gene_structure$n_sites != expected_sites_per_gene)) {
  bad <- gene_structure[n_sites != expected_sites_per_gene]
  fail_step(STEP, paste0("Frozen panel no longer has exactly ", expected_sites_per_gene, " sites per gene. Example: ",
                         paste(utils::head(paste0(bad$gene_id, "=", bad$n_sites), 10L), collapse = ", ")))
}

# Gene ranking should be one-to-one and match between panel and sites.
panel[, final_prepopulation_rank := as.integer(final_prepopulation_rank)]
if (anyNA(panel$final_prepopulation_rank)) fail_step(STEP, "final_prepopulation_rank contains non-integer/missing values in gene panel.")
if (data.table::uniqueN(panel$final_prepopulation_rank) != nrow(panel)) fail_step(STEP, "final_prepopulation_rank is not unique across genes.")
rank_check <- merge(
  panel[, .(gene_id, panel_rank = final_prepopulation_rank)],
  sites[, .(site_rank = unique(final_prepopulation_rank)), by = gene_id],
  by = "gene_id", all = TRUE
)
if (any(rank_check$panel_rank != rank_check$site_rank, na.rm = TRUE) || anyNA(rank_check$panel_rank) || anyNA(rank_check$site_rank)) {
  fail_step(STEP, "Gene-level and site-level final_prepopulation_rank values do not agree.")
}

sites[, discovery_phased_panel_assessable := genomic_seqid %in% canonical]
sites[, discovery_unresolved_reason := data.table::fifelse(
  discovery_phased_panel_assessable,
  NA_character_,
  "CONTIG_ABSENT_FROM_PHASE2_AR1_PHASED_PANEL"
)]
# Backward-compatible aliases are retained only so older diagnostic notebooks do
# not break. New analysis scripts use the source-neutral discovery_* columns.
sites[, ag3_phased_panel_assessable := discovery_phased_panel_assessable]
sites[, ag3_unresolved_reason := discovery_unresolved_reason]

# Expected frozen counts are an integrity guardrail, not merely a note.
actual_counts <- list(
  genes = data.table::uniqueN(panel$gene_id),
  sites = nrow(sites),
  canonical_arm_sites = sum(sites$discovery_phased_panel_assessable),
  unplaced_sites = sum(!sites$discovery_phased_panel_assessable)
)
for (nm in intersect(names(expected), names(actual_counts))) {
  ex <- suppressWarnings(as.integer(expected[[nm]]))
  if (is.finite(ex) && actual_counts[[nm]] != ex) {
    fail_step(STEP, paste0("Frozen-input integrity check failed for '", nm, "': expected ", ex, ", observed ", actual_counts[[nm]], "."))
  }
}

# Preserve the frozen records in deterministic order.
data.table::setorder(panel, final_prepopulation_rank, gene_id)
data.table::setorder(sites, final_prepopulation_rank, gene_id, population_site_rank_within_gene, population_site_id)
atomic_fwrite(sites, "data_processed/00_frozen_target_sites.csv")
atomic_fwrite(panel, "data_processed/00_frozen_gene_panel.csv")

qc <- data.table::data.table(
  metric = c(
    "Frozen genes", "Frozen sites", "Canonical-arm sites", "Unplaced-contig sites",
    "Genes represented on canonical arms", "External benchmark genes", "Non-benchmark genes",
    "Sites per gene minimum", "Sites per gene maximum"
  ),
  value = c(
    data.table::uniqueN(panel$gene_id), nrow(sites), sum(sites$discovery_phased_panel_assessable),
    sum(!sites$discovery_phased_panel_assessable),
    data.table::uniqueN(sites[discovery_phased_panel_assessable == TRUE, gene_id]),
    sum(panel$is_external_benchmark), sum(!panel$is_external_benchmark),
    min(gene_structure$n_sites), max(gene_structure$n_sites)
  )
)
atomic_fwrite(qc, "data_processed/00_preflight_qc.csv")
write_checksum(
  c(required, "data_processed/00_frozen_target_sites.csv", "data_processed/00_frozen_gene_panel.csv"),
  "logs/00_checksums.tsv"
)
write_session_info("logs/00_sessionInfo.txt")
print(qc)
log_step(STEP, "Preflight completed successfully")
