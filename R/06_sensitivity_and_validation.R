# STEP 06: sensitivity analyses, geographic heterogeneity, and internal reproducibility validation

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 06 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "06"
log_step(STEP, "Running sensitivity, heterogeneity, and validation analyses")

required <- c(
  "data_processed/02_analysis_population_manifest.csv",
  "data_processed/03_site_population_metrics.csv",
  "data_processed/03_population_site_metrics.csv",
  "data_processed/04_multiplex_portfolio_metrics.csv",
  "data_processed/04_multiplex_portfolio_by_population.csv",
  "data_processed/05_gene_population_robustness.csv",
  "data_processed/00_frozen_gene_panel.csv"
)
assert_files(required, nonempty = TRUE)

cfg <- read_config()
site_w <- unlist(cfg_get(cfg, "site_robustness_weights", required = TRUE), use.names = TRUE)
port_w <- unlist(cfg_get(cfg, "multiplex_robustness_weights", required = TRUE), use.names = TRUE)
gene_w <- unlist(cfg_get(cfg, "gene_robustness_weights", required = TRUE), use.names = TRUE)
validate_weights(site_w, "site_robustness_weights")
validate_weights(port_w, "multiplex_robustness_weights")
validate_weights(gene_w, "gene_robustness_weights")
site_w <- site_w[c("global_23bp_exact", "population_q10_23bp_exact", "population_breadth_ge_0_95", "global_pam_intact", "one_minus_max_alt_af_any_population", "accessibility_fraction_23bp")]
port_w <- port_w[c("global_at_least_one_intact", "population_q10_at_least_one_intact", "population_breadth_ge_0_95")]
gene_w <- gene_w[c("best_single_site", "best_pair", "best_triple")]
if (anyNA(site_w) || anyNA(port_w) || anyNA(gene_w)) fail_step(STEP, "One or more required sensitivity weights are missing from analysis_config.yml.")

robust_th <- cfg_get(cfg, "site_class_thresholds", "robust", required = TRUE)
intermediate_th <- cfg_get(cfg, "site_class_thresholds", "intermediate", required = TRUE)
accessibility_min <- as.numeric(cfg_get(cfg, "accessibility_min_for_classification", default = 0.90))

pops <- data.table::fread(required[[1]], na.strings = c("", "NA", "NaN"))
site_pop <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
site <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))
portfolio_base <- data.table::fread(required[[4]], na.strings = c("", "NA", "NaN"))
portfolio_pop <- data.table::fread(required[[5]], na.strings = c("", "NA", "NaN"))
gene_primary <- data.table::fread(required[[6]], na.strings = c("", "NA", "NaN"))
panel <- data.table::fread(required[[7]], na.strings = c("", "NA", "NaN"))

assert_columns(pops, c("population_id", "n_samples", "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30"), "analysis population manifest")
assert_columns(site_pop, c(
  "population_site_id", "population_id", "n_samples",
  "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30",
  "target_23bp_exact_match_fraction", "max_target_variant_alt_af"
), "site-population metrics")
assert_columns(site, c(
  "population_site_id", "gene_id", "population_metric_status", "population_site_rank_within_gene",
  "target_23bp_exact_match_fraction", "pam_intact_fraction", "accessibility_fraction_23bp",
  "population_target_site_robustness_score", "population_robustness_class"
), "site metrics")
assert_columns(gene_primary, c("gene_id", "final_prepopulation_rank", "population_filtered_prepopulation_rank", "population_filtered_rank_shift", "gene_population_robustness_score", "population_deployment_class"), "primary gene metrics")
assert_columns(portfolio_pop, c(
  "gene_id", "portfolio_size", "portfolio_site_ids", "population_id", "n_samples",
  "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30",
  "at_least_one_intact_fraction"
), "portfolio-by-population metrics")
assert_unique(site, "population_site_id", "site metrics")
assert_unique(panel, "gene_id", "frozen gene panel")
assert_unique(gene_primary, "gene_id", "primary gene metrics")

# 1) Geographic heterogeneity under the primary n>=20 population definition.
sp20 <- site_pop[primary_population_n20 == TRUE]
hetero <- sp20[, .(
  population_exact_mean = safe_mean(target_23bp_exact_match_fraction),
  population_exact_sd = safe_sd(target_23bp_exact_match_fraction),
  population_exact_min = safe_min(target_23bp_exact_match_fraction),
  population_exact_max = safe_max(target_23bp_exact_match_fraction),
  populations_assessed = sum(is.finite(target_23bp_exact_match_fraction))
), by = population_site_id]
hetero[, population_exact_range := population_exact_max - population_exact_min]
atomic_fwrite(hetero, "data_processed/06_site_geographic_heterogeneity.csv")

# 2) Threshold sensitivity for the ROBUST site definition.
thresholds <- data.table::data.table(
  scenario = c("strict", "primary", "permissive"),
  robust_global = c(0.99, as.numeric(robust_th$global_23bp_exact_min), 0.90),
  robust_q10 = c(0.95, as.numeric(robust_th$population_q10_23bp_exact_min), 0.80),
  robust_pam = c(0.995, as.numeric(robust_th$global_pam_intact_min), 0.98),
  max_af = c(0.02, as.numeric(robust_th$max_alt_af_any_population_max), 0.10)
)
class_rows <- data.table::rbindlist(lapply(seq_len(nrow(thresholds)), function(i) {
  th <- thresholds[i]
  x <- site[, .(
    population_site_id, gene_id, population_metric_status, accessibility_fraction_23bp,
    target_23bp_exact_match_fraction, population_q10_target_23bp_exact,
    pam_intact_fraction, max_target_variant_alt_af
  )]
  x[, scenario := th$scenario]
  x[, robust_under_scenario :=
      population_metric_status == "REAL_PHASE2_DISCOVERY_DATA" &
      is.finite(accessibility_fraction_23bp) & accessibility_fraction_23bp >= accessibility_min &
      is.finite(target_23bp_exact_match_fraction) & target_23bp_exact_match_fraction >= th$robust_global &
      is.finite(population_q10_target_23bp_exact) & population_q10_target_23bp_exact >= th$robust_q10 &
      is.finite(pam_intact_fraction) & pam_intact_fraction >= th$robust_pam &
      is.finite(max_target_variant_alt_af) & max_target_variant_alt_af <= th$max_af]
  x[, .(population_site_id, gene_id, scenario, robust_under_scenario)]
}), fill = TRUE)
atomic_fwrite(class_rows, "data_processed/06_site_classification_sensitivity.csv")

# Population-size thresholds. Always include the prespecified primary n=20 plus n=10/30.
requested_n <- sort(unique(c(
  10L, 20L, 30L,
  as.integer(cfg_get(cfg, "primary_min_population_n", default = 20L)),
  as.integer(unlist(cfg_get(cfg, "sensitivity_min_population_n", default = c(10L, 30L))))
)))
requested_n <- requested_n[is.finite(requested_n) & requested_n >= 1L]

row_population_flag_for_n <- function(dt, nmin) {
  exact <- if (nmin == 20L) "primary_population_n20" else paste0("sensitivity_population_n", nmin)
  if (exact %in% names(dt)) return(normalize_flag(dt[[exact]]))
  if (!"n_samples" %in% names(dt)) {
    fail_step(STEP, paste0("Cannot recompute n>=", nmin, " sensitivity because row-level n_samples is unavailable."))
  }
  as.integer(dt$n_samples) >= as.integer(nmin)
}

classify_site <- function(status, acc, ge, q10, pam, maxaf, n_pop) {
  if (identical(status, "UNRESOLVED_CONTIG_NOT_IN_PHASED_PANEL")) return("UNRESOLVED_CONTIG")
  if (!is.finite(acc) || acc < accessibility_min) return("UNRESOLVED_ACCESSIBILITY")
  if (!is.finite(n_pop) || n_pop < 1) return("UNRESOLVED_POPULATION")
  if (is.finite(ge) && ge >= as.numeric(robust_th$global_23bp_exact_min) &&
      is.finite(q10) && q10 >= as.numeric(robust_th$population_q10_23bp_exact_min) &&
      is.finite(pam) && pam >= as.numeric(robust_th$global_pam_intact_min) &&
      is.finite(maxaf) && maxaf <= as.numeric(robust_th$max_alt_af_any_population_max)) return("ROBUST")
  if (is.finite(ge) && ge >= as.numeric(intermediate_th$global_23bp_exact_min) &&
      is.finite(q10) && q10 >= as.numeric(intermediate_th$population_q10_23bp_exact_min) &&
      is.finite(pam) && pam >= as.numeric(intermediate_th$global_pam_intact_min)) return("INTERMEDIATE")
  "CONCERN"
}

classify_gene <- function(score, single_class, pair_global, pair_q10) {
  if (!is.finite(score)) return("UNRESOLVED")
  sc <- if (is.na(single_class)) "" else as.character(single_class)
  pg <- is.finite(pair_global)
  pq <- is.finite(pair_q10)
  if (identical(sc, "ROBUST") && pg && pq && pair_global >= 0.99 && pair_q10 >= 0.95) return("HIGHLY_ROBUST")
  if (nzchar(sc) && sc != "ROBUST" && pg && pq && pair_global >= 0.95 && pair_q10 >= 0.90) return("MULTIPLEX_RESCUABLE")
  if (sc %in% c("ROBUST", "INTERMEDIATE") || (pg && pair_global >= 0.90)) return("INTERMEDIATE")
  "CONCERN"
}

size_site_rows <- list()
size_port_rows <- list()
size_gene_rows <- list()
priority <- c(HIGHLY_ROBUST = 1L, MULTIPLEX_RESCUABLE = 2L, INTERMEDIATE = 3L, CONCERN = 4L, UNRESOLVED = 5L)

for (nmin in requested_n) {
  # IMPORTANT: population-size eligibility is contig-specific, not global.
  # The site/portfolio population tables already carry the denominator and
  # eligibility flags derived from the appropriate chromosome-arm phased panel.
  # This matters particularly for X, where the phased sample set is smaller.
  site_eligible <- row_population_flag_for_n(site_pop, nmin)
  if (!any(site_eligible, na.rm = TRUE)) {
    log_warning(STEP, paste0("No contig-specific site-population rows meet n>=", nmin, "; site sensitivity outputs will be unresolved."))
  }

  # Site recalculation using contig-specific eligibility.
  recal <- site_pop[site_eligible == TRUE, .(
    q10 = safe_quantile(target_23bp_exact_match_fraction, 0.10),
    breadth = if (any(is.finite(target_23bp_exact_match_fraction))) mean(target_23bp_exact_match_fraction[is.finite(target_23bp_exact_match_fraction)] >= 0.95) else NA_real_,
    maxaf = if (any(is.finite(max_target_variant_alt_af))) safe_max(max_target_variant_alt_af) else NA_real_,
    n_pop = sum(is.finite(target_23bp_exact_match_fraction))
  ), by = population_site_id]

  ss <- merge(site, recal, by = "population_site_id", all.x = TRUE, sort = FALSE)
  ss[is.na(n_pop), n_pop := 0L]
  ss[, population_q10_sensitivity := q10]
  ss[, population_breadth_sensitivity := breadth]
  ss[, max_alt_af_sensitivity := maxaf]
  ss[, site_score_sensitivity := mapply(
    function(ge, q, b, pa, ma, ac) {
      vals <- c(
        global_23bp_exact = ge,
        population_q10_23bp_exact = q,
        population_breadth_ge_0_95 = b,
        global_pam_intact = pa,
        one_minus_max_alt_af_any_population = if (is.finite(ma)) 1 - ma else NA_real_,
        accessibility_fraction_23bp = ac
      )
      renorm_weighted(vals, site_w[names(vals)])
    },
    target_23bp_exact_match_fraction, population_q10_sensitivity, population_breadth_sensitivity,
    pam_intact_fraction, max_alt_af_sensitivity, accessibility_fraction_23bp
  )]
  ss[, class_sensitivity := mapply(
    classify_site,
    population_metric_status, accessibility_fraction_23bp, target_23bp_exact_match_fraction,
    population_q10_sensitivity, pam_intact_fraction, max_alt_af_sensitivity, n_pop,
    USE.NAMES = FALSE
  )]
  ss[, min_population_n := as.integer(nmin)]
  size_site_rows[[as.character(nmin)]] <- ss[, .(
    population_site_id, gene_id, min_population_n,
    populations_assessed_sensitivity = n_pop,
    population_q10_sensitivity, population_breadth_sensitivity,
    max_alt_af_sensitivity, site_score_sensitivity, class_sensitivity
  )]

  # Portfolio recalculation using the contig-specific denominator/eligibility
  # already attached by Step 04 for each gene's chromosome arm.
  if (nrow(portfolio_base)) {
    portfolio_eligible <- row_population_flag_for_n(portfolio_pop, nmin)
    pp <- portfolio_pop[portfolio_eligible == TRUE, .(
      q10 = safe_quantile(at_least_one_intact_fraction, 0.10),
      breadth = if (any(is.finite(at_least_one_intact_fraction))) mean(at_least_one_intact_fraction[is.finite(at_least_one_intact_fraction)] >= 0.95) else NA_real_,
      n_pop = sum(is.finite(at_least_one_intact_fraction))
    ), by = .(gene_id, portfolio_size, portfolio_site_ids)]

    pb <- merge(
      portfolio_base[, .(
        gene_id, portfolio_size, portfolio_site_ids,
        global_at_least_one_intact_fraction, global_all_targets_intact_fraction,
        portfolio_span_bp
      )],
      pp,
      by = c("gene_id", "portfolio_size", "portfolio_site_ids"), all.x = TRUE, sort = FALSE
    )
    pb[is.na(n_pop), n_pop := 0L]
    pb[, portfolio_score_sensitivity := mapply(function(g, q, b) {
      vals <- c(
        global_at_least_one_intact = g,
        population_q10_at_least_one_intact = q,
        population_breadth_ge_0_95 = b
      )
      renorm_weighted(vals, port_w[names(vals)])
    }, global_at_least_one_intact_fraction, q10, breadth)]
    pb[, `:=`(
      .sort_score = data.table::fifelse(is.finite(portfolio_score_sensitivity), portfolio_score_sensitivity, -Inf),
      .sort_q10 = data.table::fifelse(is.finite(q10), q10, -Inf),
      .sort_global = data.table::fifelse(is.finite(global_at_least_one_intact_fraction), global_at_least_one_intact_fraction, -Inf),
      .sort_breadth = data.table::fifelse(is.finite(breadth), breadth, -Inf)
    )]
    data.table::setorder(
      pb, gene_id, portfolio_size,
      -.sort_score, -.sort_q10, -.sort_global, -.sort_breadth,
      portfolio_span_bp, portfolio_site_ids
    )
    pb[, rank_sensitivity := seq_len(.N), by = .(gene_id, portfolio_size)]
    pb[, c(".sort_score", ".sort_q10", ".sort_global", ".sort_breadth") := NULL]
    pb[, min_population_n := as.integer(nmin)]
    size_port_rows[[as.character(nmin)]] <- pb

    bestp <- pb[portfolio_size == 2L & rank_sensitivity == 1L, .(
      gene_id,
      best_pair_score = portfolio_score_sensitivity,
      best_pair_global = global_at_least_one_intact_fraction,
      best_pair_q10 = q10
    )]
    bestt <- pb[portfolio_size == 3L & rank_sensitivity == 1L, .(
      gene_id,
      best_triple_score = portfolio_score_sensitivity
    )]
  } else {
    pb <- data.table::data.table()
    size_port_rows[[as.character(nmin)]] <- data.table::data.table(
      gene_id = character(), portfolio_size = integer(), portfolio_site_ids = character(),
      global_at_least_one_intact_fraction = numeric(), global_all_targets_intact_fraction = numeric(),
      portfolio_span_bp = integer(), q10 = numeric(), breadth = numeric(), n_pop = integer(),
      portfolio_score_sensitivity = numeric(), rank_sensitivity = integer(), min_population_n = integer()
    )
    bestp <- data.table::data.table(
      gene_id = character(), best_pair_score = numeric(), best_pair_global = numeric(), best_pair_q10 = numeric()
    )
    bestt <- data.table::data.table(gene_id = character(), best_triple_score = numeric())
  }

  # Gene recalculation. Unresolved site classes are never eligible as a best single site.
  sg <- ss[
    is.finite(site_score_sensitivity) & class_sensitivity %in% c("ROBUST", "INTERMEDIATE", "CONCERN"),
    .(gene_id, population_site_rank_within_gene, site_score_sensitivity, class_sensitivity)
  ]
  data.table::setorder(sg, gene_id, -site_score_sensitivity, population_site_rank_within_gene)
  bs <- if (nrow(sg)) {
    sg[, .SD[1L], by = gene_id][, .(
      gene_id, best_single_score = site_score_sensitivity, best_single_class = class_sensitivity
    )]
  } else {
    data.table::data.table(gene_id = character(), best_single_score = numeric(), best_single_class = character())
  }

  gg <- merge(panel, bs, by = "gene_id", all.x = TRUE, sort = FALSE)
  gg <- merge(gg, bestp, by = "gene_id", all.x = TRUE, sort = FALSE)
  gg <- merge(gg, bestt, by = "gene_id", all.x = TRUE, sort = FALSE)
  if (nrow(gg) != nrow(panel)) fail_step(STEP, paste0("Gene sensitivity merge changed row count at n>=", nmin, "."))

  gg[, gene_score_sensitivity := mapply(
    function(s, p, t) renorm_weighted(c(s, p, t), c(gene_w[["best_single_site"]], gene_w[["best_pair"]], gene_w[["best_triple"]])),
    best_single_score, best_pair_score, best_triple_score
  )]
  gg[, class_sensitivity := mapply(
    classify_gene, gene_score_sensitivity, best_single_class, best_pair_global, best_pair_q10,
    USE.NAMES = FALSE
  )]
  gg[, class_priority := unname(priority[class_sensitivity])]
  if (anyNA(gg$class_priority)) fail_step(STEP, paste0("Unknown gene class generated at n>=", nmin, "."))
  data.table::setorder(gg, class_priority, final_prepopulation_rank, -gene_score_sensitivity, gene_id)
  gg[, population_filtered_rank_sensitivity := seq_len(.N)]
  gg[, min_population_n := as.integer(nmin)]
  size_gene_rows[[as.character(nmin)]] <- gg[, .(
    gene_id, min_population_n, gene_score_sensitivity, class_sensitivity,
    population_filtered_rank_sensitivity
  )]
}

site_size <- data.table::rbindlist(size_site_rows, fill = TRUE, use.names = TRUE)
port_size <- data.table::rbindlist(size_port_rows, fill = TRUE, use.names = TRUE)
gene_size <- data.table::rbindlist(size_gene_rows, fill = TRUE, use.names = TRUE)
atomic_fwrite(site_size, "data_processed/06_population_size_site_sensitivity.csv")
atomic_fwrite(port_size, "data_processed/06_population_size_portfolio_sensitivity.csv")
atomic_fwrite(gene_size, "data_processed/06_population_size_gene_sensitivity.csv")

log_step(STEP, "Population-size sensitivity uses contig-specific denominators carried from Steps 03/04 (including X-specific sample availability).")

# 3) Reproducibility check: n>=20 sensitivity recomputation must reproduce the primary analysis.
validation_rows <- list()
if (20L %in% requested_n) {
  s20 <- site_size[min_population_n == 20L]
  scheck <- merge(
    site[, .(population_site_id, primary_score = population_target_site_robustness_score, primary_class = population_robustness_class)],
    s20[, .(population_site_id, recalculated_score = site_score_sensitivity, recalculated_class = class_sensitivity)],
    by = "population_site_id", all = TRUE
  )
  scheck[, abs_score_difference := abs(primary_score - recalculated_score)]
  scheck[, class_matches := primary_class == recalculated_class]
  atomic_fwrite(scheck, "data_processed/06_primary_site_reproducibility_check.csv")

  g20 <- gene_size[min_population_n == 20L]
  gcheck <- merge(
    gene_primary[, .(
      gene_id, primary_gene_score = gene_population_robustness_score,
      primary_class = population_deployment_class,
      primary_rank = population_filtered_prepopulation_rank
    )],
    g20[, .(
      gene_id, recalculated_gene_score = gene_score_sensitivity,
      recalculated_class = class_sensitivity,
      recalculated_rank = population_filtered_rank_sensitivity
    )],
    by = "gene_id", all = TRUE
  )
  gcheck[, abs_score_difference := abs(primary_gene_score - recalculated_gene_score)]
  gcheck[, class_matches := primary_class == recalculated_class]
  gcheck[, rank_matches := primary_rank == recalculated_rank]
  atomic_fwrite(gcheck, "data_processed/06_primary_gene_reproducibility_check.csv")

  max_site_diff <- safe_max(scheck$abs_score_difference, default = 0)
  max_gene_diff <- safe_max(gcheck$abs_score_difference, default = 0)
  class_site_ok <- all(scheck$class_matches[!is.na(scheck$class_matches)])
  class_gene_ok <- all(gcheck$class_matches[!is.na(gcheck$class_matches)])
  rank_gene_ok <- all(gcheck$rank_matches[!is.na(gcheck$rank_matches)])

  validation_rows[["primary"]] <- data.table::data.table(
    metric = c(
      "Max absolute site score difference: primary vs n20 recomputation",
      "Max absolute gene score difference: primary vs n20 recomputation",
      "All comparable site classes reproduce at n20",
      "All comparable gene classes reproduce at n20",
      "All comparable population-filtered gene ranks reproduce at n20"
    ),
    value = c(max_site_diff, max_gene_diff, as.numeric(class_site_ok), as.numeric(class_gene_ok), as.numeric(rank_gene_ok))
  )

  # The primary analysis is defined by the same contig-specific n>=20 rule.
  # Therefore the sensitivity recomputation must reproduce it numerically, not
  # merely reproduce classes/ranks. Treat any material mismatch as a pipeline
  # consistency failure rather than silently carrying it into external validation.
  repro_tol <- 1e-10
  if (is.finite(max_site_diff) && max_site_diff > repro_tol) {
    fail_step(STEP, paste0(
      "Primary site scores do not reproduce under the contig-specific n20 recomputation; maximum absolute difference = ",
      signif(max_site_diff, 8), ". Inspect data_processed/06_primary_site_reproducibility_check.csv."
    ))
  }
  if (is.finite(max_gene_diff) && max_gene_diff > repro_tol) {
    fail_step(STEP, paste0(
      "Primary gene scores do not reproduce under the contig-specific n20 recomputation; maximum absolute difference = ",
      signif(max_gene_diff, 8), ". Inspect data_processed/06_primary_gene_reproducibility_check.csv."
    ))
  }
  if (!class_site_ok || !class_gene_ok || !rank_gene_ok) {
    fail_step(STEP, "At least one primary class/rank failed to reproduce exactly at contig-specific n>=20.")
  }
}

# 4) Rank concordance summaries.
rho_primary <- safe_cor(
  gene_primary$final_prepopulation_rank,
  gene_primary$population_filtered_prepopulation_rank,
  method = "spearman"
)
rank_sens_summary <- data.table::rbindlist(lapply(setdiff(requested_n, 20L), function(nm) {
  g20 <- gene_size[min_population_n == 20L, .(gene_id, rank20 = population_filtered_rank_sensitivity)]
  xx <- merge(
    g20,
    gene_size[min_population_n == nm, .(gene_id, rank_other = population_filtered_rank_sensitivity, class_other = class_sensitivity)],
    by = "gene_id", all = FALSE
  )
  data.table::data.table(
    metric = paste0("Spearman population-filtered rank n20 vs n", nm),
    value = safe_cor(xx$rank20, xx$rank_other, method = "spearman")
  )
}), fill = TRUE)

summary <- data.table::rbindlist(list(
  data.table::data.table(
    metric = c(
      "Spearman frozen vs population-filtered rank (primary n20)",
      "Median absolute rank shift (primary n20)",
      "Maximum absolute rank shift (primary n20)"
    ),
    value = c(
      rho_primary,
      safe_median(abs(gene_primary$population_filtered_rank_shift)),
      safe_max(abs(gene_primary$population_filtered_rank_shift))
    )
  ),
  rank_sens_summary,
  data.table::rbindlist(validation_rows, fill = TRUE)
), fill = TRUE)
atomic_fwrite(summary, "data_processed/06_rank_validation_summary.csv")

# 5) Class transitions across population-size thresholds.
wide <- data.table::dcast(
  gene_size,
  gene_id ~ min_population_n,
  value.var = "class_sensitivity",
  fun.aggregate = function(z) if (length(z)) z[[1]] else NA_character_
)
atomic_fwrite(wide, "data_processed/06_population_size_class_transitions.csv")

# 6) Benchmark-focused sensitivity table.
bench_ids <- panel[normalize_flag(is_external_benchmark) == TRUE, gene_id]
benchmark_sensitivity <- merge(
  panel[gene_id %in% bench_ids],
  gene_size[gene_id %in% bench_ids],
  by = "gene_id", all.x = TRUE, sort = FALSE
)
atomic_fwrite(benchmark_sensitivity, "data_processed/06_benchmark_population_sensitivity.csv")

write_checksum(
  c(
    "data_processed/06_site_geographic_heterogeneity.csv",
    "data_processed/06_site_classification_sensitivity.csv",
    "data_processed/06_population_size_site_sensitivity.csv",
    "data_processed/06_population_size_portfolio_sensitivity.csv",
    "data_processed/06_population_size_gene_sensitivity.csv",
    "data_processed/06_rank_validation_summary.csv",
    "data_processed/06_population_size_class_transitions.csv",
    "data_processed/06_benchmark_population_sensitivity.csv"
  ),
  "logs/06_checksums.tsv"
)
write_session_info("logs/06_sessionInfo.txt")
print(summary)
log_step(STEP, "Sensitivity and validation analyses completed successfully")
