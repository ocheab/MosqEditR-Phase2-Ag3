# STEP 06b: Build the sequence-aware Phase-2 vs Ag3 comparison table
#
# This script does NOT acquire data and does NOT change any discovery rank/score.
# It simply joins the corrected sequence-aware Ag3 taxon summary to the frozen
# Phase-2 site metrics using population_site_id.

if (!file.exists("R/helpers.R")) stop("Run Step 06b from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table"))
suppressPackageStartupMessages(library(data.table))

STEP <- "06b"
log_step(STEP, "Building sequence-aware Phase-2 vs Ag3 external-validation comparison")

phase2_file <- "data_processed/03_population_site_metrics.csv"
ag3_file <- "data_processed/06a_ag3_external_site_validation_by_taxon.csv"

assert_files(c(phase2_file, ag3_file), nonempty = TRUE)

p2 <- fread(phase2_file, na.strings = c("", "NA", "NaN"))
a3 <- fread(ag3_file, na.strings = c("", "NA", "NaN"))

assert_columns(
  p2,
  c(
    "population_site_id",
    "target_23bp_exact_match_fraction",
    "population_q10_target_23bp_exact",
    "pam_intact_fraction",
    "max_target_variant_alt_af"
  ),
  "Phase-2 site metrics"
)

assert_columns(
  a3,
  c(
    "population_site_id", "gene_id", "taxon",
    "strict_genotype_exact_23bp_fraction",
    "protospacer_exact_20bp_fraction",
    "pam_exact_3bp_fraction",
    "pam_ngg_intact_fraction",
    "functional_target_intact_fraction",
    "max_position_nonreference_allele_fraction_23bp"
  ),
  "sequence-aware Ag3 taxon validation"
)

assert_unique(p2, "population_site_id", "Phase-2 site metrics")
if (anyDuplicated(a3[, .(population_site_id, taxon)])) {
  stop("Duplicate population_site_id x taxon rows in Ag3 sequence-aware validation.", call. = FALSE)
}

cmp <- merge(
  a3,
  p2[, .(
    population_site_id,
    phase2_global_exact = target_23bp_exact_match_fraction,
    phase2_population_q10_exact = population_q10_target_23bp_exact,
    phase2_pam_intact = pam_intact_fraction,
    phase2_max_alt_af = max_target_variant_alt_af
  )],
  by = "population_site_id",
  all.x = TRUE,
  sort = FALSE
)

setcolorder(
  cmp,
  c(
    "population_site_id", "gene_id", "taxon",
    "phase2_global_exact",
    "phase2_population_q10_exact",
    "phase2_pam_intact",
    "phase2_max_alt_af",
    "strict_genotype_exact_23bp_fraction",
    "protospacer_exact_20bp_fraction",
    "pam_exact_3bp_fraction",
    "pam_ngg_intact_fraction",
    "functional_target_intact_fraction",
    "max_position_nonreference_allele_fraction_23bp",
    setdiff(
      names(cmp),
      c(
        "population_site_id", "gene_id", "taxon",
        "phase2_global_exact",
        "phase2_population_q10_exact",
        "phase2_pam_intact",
        "phase2_max_alt_af",
        "strict_genotype_exact_23bp_fraction",
        "protospacer_exact_20bp_fraction",
        "pam_exact_3bp_fraction",
        "pam_ngg_intact_fraction",
        "functional_target_intact_fraction",
        "max_position_nonreference_allele_fraction_23bp"
      )
    )
  )
)

# Convenience aliases for legacy exploratory scripts, while preserving the
# explicit corrected metrics alongside them.
cmp[, ag3_exact := strict_genotype_exact_23bp_fraction]
cmp[, ag3_pam := pam_ngg_intact_fraction]
cmp[, ag3_nonref_af := max_position_nonreference_allele_fraction_23bp]

out1 <- "data_processed/06b_ag3_phase2_external_validation_comparison_sequence_aware.csv"
out2 <- "ag3_phase2_external_validation_comparison_SEQUENCE_AWARE.csv"

atomic_fwrite(cmp, out1)
atomic_fwrite(cmp, out2)

qc <- data.table(
  metric = c(
    "Rows",
    "Unique locked sites",
    "Taxa",
    "Missing sequence-aware exact",
    "Missing sequence-aware functional PAM",
    "Missing maximum-position nonreference AF"
  ),
  value = as.character(c(
    nrow(cmp),
    uniqueN(cmp$population_site_id),
    paste(sort(unique(cmp$taxon)), collapse = ","),
    sum(!is.finite(cmp$strict_genotype_exact_23bp_fraction)),
    sum(!is.finite(cmp$pam_ngg_intact_fraction)),
    sum(!is.finite(cmp$max_position_nonreference_allele_fraction_23bp))
  ))
)
atomic_fwrite(qc, "data_processed/06b_ag3_phase2_comparison_qc.csv")
print(qc)

write_session_info("logs/06b_sessionInfo.txt")
log_step(STEP, paste0(
  "Sequence-aware comparison written: ", out2,
  " (", nrow(cmp), " rows; ", uniqueN(cmp$population_site_id), " sites)"
))
