# STEP 06a: Ag3.0 site-level external validation of the LOCKED Phase-2 discovery targets
#
# CRITICAL: this script reports validation only. It MUST NOT overwrite Phase-2
# discovery scores/classes/ranks and MUST NOT feed Ag3 values back into Steps 03-05.

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 06a from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))
suppressPackageStartupMessages(library(data.table))

STEP <- "06a"
log_step(STEP, "Evaluating locked Phase-2 discoveries in independent public Ag3 genotype data")

required <- c(
  "data_processed/05b_external_validation_targets.csv",
  "data_processed/03_population_site_metrics.csv",
  "data_processed/05_gene_population_robustness.csv"
)
assert_files(required, nonempty = TRUE)

cfg <- read_config()
ev <- cfg_get(cfg, "external_validation", required = TRUE)
out_base <- as.character(ev$output_dir %||% "data_raw/ag3_external")
min_called <- as.integer(ev$minimum_called_samples_per_site %||% 100L)
exact_th <- as.numeric(ev$strict_genotype_conservation_threshold_23bp %||% ev$conservation_threshold_23bp %||% 0.95)
pam_th <- as.numeric(ev$pam_ngg_threshold %||% ev$pam_threshold %||% 0.99)
maxpos_nonref_th <- as.numeric(ev$max_position_nonreference_af %||% ev$max_nonreference_allele_fraction %||% 0.05)
if (!isTRUE(as.logical(ev$require_no_discovery_retuning %||% TRUE))) {
  fail_step(STEP, "External-validation isolation must remain TRUE; Ag3 cannot retune discovery.")
}

# Prefer final validation if present. Pilot is explicitly labelled and cannot be
# promoted to final support when it does not meet the prespecified sample threshold.
final_dir <- file.path(out_base, "final")
pilot_dir <- file.path(out_base, "pilot")
if (file.exists(file.path(final_dir, "ACQUISITION_COMPLETE.ok"))) {
  run_dir <- final_dir; run_mode <- "final"
} else if (file.exists(file.path(pilot_dir, "ACQUISITION_COMPLETE.ok"))) {
  run_dir <- pilot_dir; run_mode <- "pilot"
  log_warning(STEP, "Only the Ag3 pilot is available. Results will be labelled PILOT and are not a substitute for final external validation.")
} else {
  fail_step(STEP, paste0(
    "No completed Ag3 external-validation acquisition found under ", out_base, ".\n",
    "Run scripts/acquire_ag3_external_validation.ps1 --pilot first, then --final after the pilot passes."
  ))
}

site_file <- file.path(run_dir, "site_summary.tsv")
taxon_file <- file.path(run_dir, "site_taxon_summary.tsv")
manifest_file <- file.path(run_dir, "run_manifest.json")
assert_files(c(site_file, taxon_file, manifest_file), nonempty = TRUE)

lock <- fread(required[[1]], na.strings = c("", "NA", "NaN"))
disc <- fread(required[[2]], na.strings = c("", "NA", "NaN"))
gene <- fread(required[[3]], na.strings = c("", "NA", "NaN"))
ext <- fread(site_file, sep = "\t", na.strings = c("", "NA", "NaN"))
ext_tax <- fread(taxon_file, sep = "\t", na.strings = c("", "NA", "NaN"))

assert_columns(lock, c("population_site_id", "gene_id", "selection_reason", "selection_locked_before_ag3"), "external target lock")
assert_columns(ext, c(
  "population_site_id", "n_samples_total", "n_samples_called_23bp", "n_samples_exact_23bp",
  "n_samples_pam_called", "n_samples_pam_intact", "exact_23bp_fraction", "pam_intact_fraction",
  "strict_genotype_exact_23bp_fraction", "protospacer_exact_20bp_fraction",
  "pam_ngg_intact_fraction", "pam_exact_3bp_fraction", "functional_target_intact_fraction",
  "called_alleles_23bp", "nonreference_alleles_23bp", "nonreference_allele_fraction_23bp",
  "max_position_nonreference_allele_fraction_23bp"
), "Ag3 external site summary")
assert_columns(ext_tax, c(
  "population_site_id", "taxon", "n_samples_called_23bp",
  "strict_genotype_exact_23bp_fraction", "pam_ngg_intact_fraction",
  "max_position_nonreference_allele_fraction_23bp"
), "Ag3 external taxon summary")
assert_columns(disc, c("population_site_id", "population_target_site_robustness_score", "population_robustness_class", "target_23bp_exact_match_fraction", "population_q10_target_23bp_exact", "pam_intact_fraction", "max_target_variant_alt_af"), "Phase-2 discovery site metrics")
assert_unique(lock, "population_site_id", "external target lock")
assert_unique(ext, "population_site_id", "Ag3 external site summary")
if (any(!normalize_flag(lock$selection_locked_before_ag3), na.rm = TRUE)) fail_step(STEP, "External target lock is not immutable/locked for every site.")
if (!setequal(lock$population_site_id, ext$population_site_id)) {
  fail_step(STEP, paste0(
    "Ag3 external summary does not match the locked validation target set. Missing=",
    paste(head(setdiff(lock$population_site_id, ext$population_site_id), 20L), collapse = ", "),
    "; extra=", paste(head(setdiff(ext$population_site_id, lock$population_site_id), 20L), collapse = ", ")
  ))
}

wilson_ci <- function(x, n, conf = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0 || x < 0 || x > n) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2*n)) / den
  half <- z * sqrt((p*(1-p)/n) + z^2/(4*n^2)) / den
  c(max(0, ctr-half), min(1, ctr+half))
}

ext[, c("exact_23bp_ci_low", "exact_23bp_ci_high") := {
  ci <- t(mapply(wilson_ci, n_samples_exact_23bp, n_samples_called_23bp))
  list(ci[,1], ci[,2])
}]
ext[, c("pam_intact_ci_low", "pam_intact_ci_high") := {
  ci <- t(mapply(wilson_ci, n_samples_pam_intact, n_samples_pam_called))
  list(ci[,1], ci[,2])
}]

ext[, external_sample_threshold_met := n_samples_called_23bp >= min_called]
ext[, external_stable_observed :=
      external_sample_threshold_met &
      is.finite(strict_genotype_exact_23bp_fraction) & strict_genotype_exact_23bp_fraction >= exact_th &
      is.finite(pam_ngg_intact_fraction) & pam_ngg_intact_fraction >= pam_th &
      is.finite(max_position_nonreference_allele_fraction_23bp) &
        max_position_nonreference_allele_fraction_23bp <= maxpos_nonref_th]
ext[, external_stable_conservative_ci :=
      external_sample_threshold_met &
      is.finite(exact_23bp_ci_low) & exact_23bp_ci_low >= exact_th &
      is.finite(pam_intact_ci_low) & pam_intact_ci_low >= pam_th &
      is.finite(max_position_nonreference_allele_fraction_23bp) &
        max_position_nonreference_allele_fraction_23bp <= maxpos_nonref_th]

# Taxon-stratified support. Require each configured primary taxon to contribute
# at least half of the overall prespecified minimum (minimum 20) for the strict
# taxon-balanced indicator.
primary_taxa <- tolower(as.character(ev$primary_taxa %||% c("gambiae", "coluzzii")))
taxon_min <- max(20L, ceiling(min_called / max(1L, length(primary_taxa))))
ext_tax[, taxon := tolower(as.character(taxon))]
ext_tax[, taxon_observed_stable :=
          n_samples_called_23bp >= taxon_min &
          is.finite(strict_genotype_exact_23bp_fraction) &
            strict_genotype_exact_23bp_fraction >= exact_th &
          is.finite(pam_ngg_intact_fraction) & pam_ngg_intact_fraction >= pam_th &
          is.finite(max_position_nonreference_allele_fraction_23bp) &
            max_position_nonreference_allele_fraction_23bp <= maxpos_nonref_th]
taxon_wide <- dcast(
  ext_tax[taxon %in% primary_taxa],
  population_site_id ~ taxon,
  value.var = c(
    "n_samples_called_23bp",
    "strict_genotype_exact_23bp_fraction",
    "pam_ngg_intact_fraction",
    "max_position_nonreference_allele_fraction_23bp",
    "taxon_observed_stable"
  )
)
if (!nrow(taxon_wide)) fail_step(STEP, "No configured primary taxa were found in Ag3 external taxon summary.")
for (tx in primary_taxa) {
  nm <- paste0("taxon_observed_stable_", tx)
  if (!nm %in% names(taxon_wide)) taxon_wide[, (nm) := FALSE]
  taxon_wide[is.na(get(nm)), (nm) := FALSE]
}
stable_cols <- paste0("taxon_observed_stable_", primary_taxa)
taxon_wide[, external_stable_taxon_balanced := Reduce(`&`, .SD), .SDcols = stable_cols]

res <- merge(lock, ext, by = c("population_site_id", "gene_id", "selection_reason"), all.x = TRUE, sort = FALSE)
res <- merge(res, taxon_wide, by = "population_site_id", all.x = TRUE, sort = FALSE)

# Phase-2 discovery metrics are attached from the canonical Step-03 discovery
# table exactly once.  Some frozen-lock versions already carry one or more of
# these convenience aliases (notably phase2_max_alt_af); if left in place,
# data.table::merge() creates .x/.y suffixes and the downstream concordance
# code cannot find the unsuffixed canonical column.
phase2_alias_cols <- c(
  "phase2_discovery_score",
  "phase2_discovery_class",
  "phase2_discovery_exact_23bp",
  "phase2_population_q10_exact",
  "phase2_discovery_pam_intact",
  "phase2_max_alt_af"
)
preexisting_phase2_aliases <- intersect(names(res), phase2_alias_cols)
if (length(preexisting_phase2_aliases)) {
  log_step(
    STEP,
    paste0(
      "Removing pre-existing Phase-2 alias column(s) before canonical discovery merge: ",
      paste(preexisting_phase2_aliases, collapse = ", ")
    )
  )
  res[, (preexisting_phase2_aliases) := NULL]
}

phase2_disc_join <- disc[, .(
  population_site_id,
  phase2_discovery_score = population_target_site_robustness_score,
  phase2_discovery_class = population_robustness_class,
  phase2_discovery_exact_23bp = target_23bp_exact_match_fraction,
  phase2_population_q10_exact = population_q10_target_23bp_exact,
  phase2_discovery_pam_intact = pam_intact_fraction,
  phase2_max_alt_af = max_target_variant_alt_af
)]

if (anyDuplicated(phase2_disc_join$population_site_id)) {
  fail_step(STEP, "Phase-2 discovery metrics are not unique by population_site_id.")
}

res <- merge(
  res,
  phase2_disc_join,
  by = "population_site_id",
  all.x = TRUE,
  sort = FALSE
)

suffix_collision_cols <- grep(
  "^phase2_.*\\.(x|y)$",
  names(res),
  value = TRUE
)
if (length(suffix_collision_cols)) {
  fail_step(
    STEP,
    paste0(
      "Unexpected Phase-2 merge suffix collision(s): ",
      paste(suffix_collision_cols, collapse = ", ")
    )
  )
}

assert_columns(
  res,
  c(
    "phase2_discovery_score",
    "phase2_discovery_class",
    "phase2_discovery_exact_23bp",
    "phase2_population_q10_exact",
    "phase2_discovery_pam_intact",
    "phase2_max_alt_af"
  ),
  "Step-06a merged Phase-2 discovery metrics"
)
res[, external_validation_run := toupper(run_mode)]
res[, external_validation_dataset := "Ag3.0 public Sanger per-sample all-site genotypes"]
res[, discovery_scores_retuned := FALSE]
res[, external_validation_status := fcase(
  !external_sample_threshold_met, "INSUFFICIENT_EXTERNAL_N",
  external_stable_conservative_ci & external_stable_taxon_balanced, "CONFIRMED_STRICT",
  external_stable_observed & external_stable_taxon_balanced, "CONFIRMED_OBSERVED",
  external_stable_observed, "CONFIRMED_OVERALL_TAXON_HETEROGENEITY",
  default = "NOT_CONFIRMED"
)]
setorder(res, final_prepopulation_rank, gene_id, population_site_rank_within_gene, population_site_id)
atomic_fwrite(res, "data_processed/06a_ag3_external_site_validation.csv")
atomic_fwrite(ext_tax, "data_processed/06a_ag3_external_site_validation_by_taxon.csv")

# Reproducible taxon-level Phase-2 vs Ag3 comparison table.
comparison_phase2_aliases <- c(
  "phase2_global_exact",
  "phase2_population_q10_exact",
  "phase2_pam_intact",
  "phase2_max_alt_af"
)
comparison_base <- copy(ext_tax)
preexisting_comparison_aliases <- intersect(
  names(comparison_base),
  comparison_phase2_aliases
)
if (length(preexisting_comparison_aliases)) {
  comparison_base[, (preexisting_comparison_aliases) := NULL]
}

comparison_phase2 <- disc[, .(
  population_site_id,
  phase2_global_exact = target_23bp_exact_match_fraction,
  phase2_population_q10_exact = population_q10_target_23bp_exact,
  phase2_pam_intact = pam_intact_fraction,
  phase2_max_alt_af = max_target_variant_alt_af
)]

comparison_taxon <- merge(
  comparison_base,
  comparison_phase2,
  by = "population_site_id",
  all.x = TRUE,
  sort = FALSE
)

comparison_suffix_collisions <- grep(
  "^phase2_.*\\.(x|y)$",
  names(comparison_taxon),
  value = TRUE
)
if (length(comparison_suffix_collisions)) {
  fail_step(
    STEP,
    paste0(
      "Unexpected taxon-comparison Phase-2 merge suffix collision(s): ",
      paste(comparison_suffix_collisions, collapse = ", ")
    )
  )
}
setcolorder(
  comparison_taxon,
  c(
    "population_site_id", "gene_id", "taxon",
    "phase2_global_exact", "phase2_population_q10_exact",
    "phase2_pam_intact", "phase2_max_alt_af",
    "strict_genotype_exact_23bp_fraction",
    "protospacer_exact_20bp_fraction",
    "pam_exact_3bp_fraction",
    "pam_ngg_intact_fraction",
    "functional_target_intact_fraction",
    "max_position_nonreference_allele_fraction_23bp",
    setdiff(
      names(comparison_taxon),
      c(
        "population_site_id", "gene_id", "taxon",
        "phase2_global_exact", "phase2_population_q10_exact",
        "phase2_pam_intact", "phase2_max_alt_af",
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
atomic_fwrite(
  comparison_taxon,
  "data_processed/06a_ag3_phase2_external_validation_comparison_sequence_aware.csv"
)

# Concordance analyses are descriptive validation diagnostics only.
# IMPORTANT: Phase-2 target_23bp_exact_match_fraction is a phased HAPLOTYPE metric,
# whereas the per-sample Ag3 all-sites VCFs are independently genotyped and unphased.
# Therefore the Ag3 strict-genotype exact fraction is NOT treated as numerically
# interchangeable with the Phase-2 haplotype exact-match fraction.
cc_exact <- res[
  is.finite(phase2_discovery_exact_23bp) &
  is.finite(strict_genotype_exact_23bp_fraction)
]
exact_rank_spearman <- if (nrow(cc_exact) >= 3L) suppressWarnings(
  cor(
    cc_exact$phase2_discovery_exact_23bp,
    cc_exact$strict_genotype_exact_23bp_fraction,
    method = "spearman"
  )
) else NA_real_

# The directly comparable cross-release diagnostic is the maximum per-position
# non-reference allele frequency across the locked 23-bp target.
cc_af <- res[
  is.finite(phase2_max_alt_af) &
  is.finite(max_position_nonreference_allele_fraction_23bp)
]
max_af_spearman <- if (nrow(cc_af) >= 3L) suppressWarnings(
  cor(
    cc_af$phase2_max_alt_af,
    cc_af$max_position_nonreference_allele_fraction_23bp,
    method = "spearman"
  )
) else NA_real_
max_af_mae <- if (nrow(cc_af)) mean(
  abs(
    cc_af$phase2_max_alt_af -
      cc_af$max_position_nonreference_allele_fraction_23bp
  )
) else NA_real_

# Gene-level external-support summary, preserving Phase-2 ranks/classes unchanged.
gene_ext <- res[, .(
  n_locked_sites = .N,
  n_sites_meeting_external_n = sum(external_sample_threshold_met, na.rm = TRUE),
  n_sites_confirmed_observed = sum(external_stable_observed, na.rm = TRUE),
  n_sites_confirmed_strict = sum(external_stable_conservative_ci, na.rm = TRUE),
  proportion_sites_confirmed_observed = mean(external_stable_observed, na.rm = TRUE),
  all_locked_sites_confirmed_observed = all(external_stable_observed %in% TRUE),
  all_locked_sites_taxon_balanced = all(external_stable_taxon_balanced %in% TRUE)
), by = gene_id]
gene_ext <- merge(
  gene[, .(gene_id, final_prepopulation_rank, population_filtered_prepopulation_rank, population_deployment_class, gene_population_robustness_score, best_single_site_id)],
  gene_ext, by = "gene_id", all.x = TRUE, sort = FALSE
)
gene_ext[, best_single_ag3_confirmed := res$external_stable_observed[match(best_single_site_id, res$population_site_id)]]
gene_ext[, best_single_ag3_confirmed_strict := res$external_stable_conservative_ci[match(best_single_site_id, res$population_site_id)]]
gene_ext[, ag3_external_validation_used_for_reranking := FALSE]
setorder(gene_ext, population_filtered_prepopulation_rank, final_prepopulation_rank, gene_id)
atomic_fwrite(gene_ext, "data_processed/06a_ag3_external_gene_support.csv")

qc <- data.table(
  metric = c(
    "External validation run mode",
    "Locked sites evaluated",
    "Sites meeting prespecified external sample threshold",
    "Sites confirmed by observed thresholds",
    "Sites confirmed by conservative Wilson-CI thresholds",
    "Sites with taxon-balanced observed confirmation",
    "Phase2 haplotype exact vs Ag3 strict-genotype exact Spearman rho (rank-only diagnostic)",
    "Phase2-vs-Ag3 maximum position nonreference AF Spearman rho",
    "Phase2-vs-Ag3 maximum position nonreference AF mean absolute difference",
    "Discovery ranks changed by Ag3",
    "Discovery thresholds changed by Ag3"
  ),
  value = as.character(c(
    run_mode,
    nrow(res),
    sum(res$external_sample_threshold_met, na.rm = TRUE),
    sum(res$external_stable_observed, na.rm = TRUE),
    sum(res$external_stable_conservative_ci, na.rm = TRUE),
    sum(res$external_stable_taxon_balanced, na.rm = TRUE),
    signif(exact_rank_spearman, 6),
    signif(max_af_spearman, 6),
    signif(max_af_mae, 6),
    0L,
    0L
  ))
)
atomic_fwrite(qc, "data_processed/06a_ag3_external_validation_qc.csv")

metric_semantics <- c(
  "Ag3 exact_23bp_fraction is a strict UNPHASED-GENOTYPE metric: a sample is exact only when every called allele at all 23 frozen target bases matches the frozen guide/PAM sequence.",
  "Ag3 protospacer_exact_20bp_fraction requires exact frozen-sequence agreement across the 20-nt protospacer.",
  "Ag3 pam_ngg_intact_fraction tests functional SpCas9 NGG preservation; variation at the PAM N base is allowed, whereas either guide-oriented G must remain G.",
  "Ag3 pam_exact_3bp_fraction is a stricter exact three-base PAM sequence metric and is reported separately.",
  "Ag3 max_position_nonreference_allele_fraction_23bp is the maximum mismatch-allele frequency at any one of the 23 frozen target positions and is the preferred cross-release allele-frequency concordance metric.",
  "Phase-2 target_23bp_exact_match_fraction is phased-haplotype based and is not numerically interchangeable with the Ag3 strict-genotype exact fraction; only a rank correlation is reported as a descriptive diagnostic."
)
atomic_writeLines(metric_semantics, "data_processed/06a_ag3_external_metric_semantics.txt")

write_checksum(
  c(site_file, taxon_file, manifest_file, "data_processed/05b_external_validation_targets.csv",
    "data_processed/06a_ag3_external_site_validation.csv", "data_processed/06a_ag3_external_gene_support.csv",
    "data_processed/06a_ag3_external_metric_semantics.txt",
    "data_processed/06a_ag3_phase2_external_validation_comparison_sequence_aware.csv"),
  "logs/06a_checksums.tsv"
)
write_session_info("logs/06a_sessionInfo.txt")
print(qc)
log_step(STEP, paste0("Ag3 external validation completed in ", toupper(run_mode), " mode without discovery retuning"))
