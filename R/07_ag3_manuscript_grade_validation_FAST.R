
#!/usr/bin/env Rscript

# =============================================================================
# STEP 07: Manuscript-grade Ag3 external validation of frozen Phase-2 targets
# =============================================================================
# Inputs:
#   ag3_phase2_external_validation_comparison_SEQUENCE_AWARE.csv
#   data_processed/06a_ag3_external_site_validation.csv
#
# Scientific interpretation:
# - Phase-2 target_23bp_exact_match_fraction is a PHASED-HAPLOTYPE metric.
# - Ag3 strict_genotype_exact_23bp_fraction is an UNPHASED-GENOTYPE metric.
# - Therefore exactness comparisons are treated as rank/concordance diagnostics,
#   not numerical interchangeability/agreement tests.
# - Phase-2 max_target_variant_alt_af vs Ag3
#   max_position_nonreference_allele_fraction_23bp is the preferred quantitative
#   cross-release allele-frequency concordance comparison.
# - Ag3 never retunes Phase-2 discovery scores, ranks, thresholds or classes.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260915)

STEP <- "07"
log_msg <- function(...) {
  cat(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    "\tINFO\t", STEP, "\t", paste0(...), "\n",
    sep = ""
  )
}

fail <- function(...) stop(paste0(...), call. = FALSE)

infile <- "ag3_phase2_external_validation_comparison_SEQUENCE_AWARE.csv"
sitefile <- "data_processed/06a_ag3_external_site_validation.csv"
outdir <- "outputs/07_ag3_manuscript_validation"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) fail("Missing input: ", infile)
if (!file.exists(sitefile)) fail("Missing input: ", sitefile)

d <- fread(infile)
site06a <- fread(sitefile)

req <- c(
  "population_site_id", "gene_id", "taxon",
  "phase2_global_exact", "phase2_population_q10_exact",
  "phase2_pam_intact", "phase2_max_alt_af",
  "strict_genotype_exact_23bp_fraction",
  "protospacer_exact_20bp_fraction",
  "pam_exact_3bp_fraction",
  "pam_ngg_intact_fraction",
  "functional_target_intact_fraction",
  "max_position_nonreference_allele_fraction_23bp",
  "n_samples_total", "n_samples_called_23bp",
  "callable_sample_fraction"
)
miss <- setdiff(req, names(d))
if (length(miss)) fail("Missing required comparison columns: ", paste(miss, collapse = ", "))

if (nrow(d) != 708L) fail("Expected 708 taxon-level rows; found ", nrow(d))
if (uniqueN(d$population_site_id) != 354L) fail("Expected 354 locked sites.")
if (!setequal(unique(d$taxon), c("gambiae", "coluzzii"))) {
  fail("Expected exactly gambiae and coluzzii taxa.")
}
if (anyDuplicated(d[, .(population_site_id, taxon)])) {
  fail("population_site_id x taxon must be unique.")
}

# Confirm Phase-2 values are invariant across Ag3 taxon rows for each site.
p2_cols <- c(
  "phase2_global_exact",
  "phase2_population_q10_exact",
  "phase2_pam_intact",
  "phase2_max_alt_af"
)
p2_invariance <- d[, lapply(.SD, uniqueN), by = population_site_id, .SDcols = p2_cols]
if (any(as.matrix(p2_invariance[, ..p2_cols]) != 1L)) {
  fail("At least one Phase-2 metric varies across duplicated Ag3 taxon rows.")
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
}

rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  sqrt(mean((y[ok] - x[ok])^2))
}

mae <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  mean(abs(y[ok] - x[ok]))
}

bias <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  mean(y[ok] - x[ok])
}

ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  vx <- var(x); vy <- var(y)
  mx <- mean(x); my <- mean(y)
  sxy <- cov(x, y)
  (2 * sxy) / (vx + vy + (mx - my)^2)
}

boot_site_ci <- function(dt, stat_fun, B = 5000L, seed = 20260915L) {
  # Fast locked-site cluster bootstrap.
  #
  # Scientific design is unchanged:
  #   - resampling unit = population_site_id
  #   - both Ag3 taxon rows for each sampled site are retained together
  #   - 5,000 bootstrap replicates
  #
  # Performance improvement:
  #   Build the site -> row-index map once, rather than filtering the full
  #   data.table 354 times within every bootstrap replicate.

  set.seed(seed)

  site_rows <- split(
    seq_len(nrow(dt)),
    dt$population_site_id,
    drop = TRUE
  )
  ids <- names(site_rows)
  n_ids <- length(ids)

  if (n_ids != 354L) {
    stop(
      paste0(
        "Bootstrap expected 354 locked sites; found ",
        n_ids
      ),
      call. = FALSE
    )
  }

  cluster_sizes <- lengths(site_rows)
  if (!all(cluster_sizes == 2L)) {
    stop(
      "Bootstrap requires exactly two taxon rows per locked site.",
      call. = FALSE
    )
  }

  vals <- numeric(B)

  for (b in seq_len(B)) {
    sampled_pos <- sample.int(
      n = n_ids,
      size = n_ids,
      replace = TRUE
    )

    row_idx <- unlist(
      site_rows[sampled_pos],
      use.names = FALSE
    )

    vals[b] <- stat_fun(dt[row_idx])

    if (b %% 1000L == 0L) {
      log_msg(
        "Bootstrap progress: ",
        b,
        "/",
        B,
        " replicates"
      )
    }
  }

  c(
    estimate = stat_fun(dt),
    ci_low = unname(
      quantile(
        vals,
        0.025,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    ci_high = unname(
      quantile(
        vals,
        0.975,
        na.rm = TRUE,
        names = FALSE
      )
    )
  )
}

save_plot <- function(p, stub, width = 7.2, height = 5.6) {
  ggsave(
    filename = file.path(outdir, paste0(stub, ".pdf")),
    plot = p, width = width, height = height, units = "in",
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(outdir, paste0(stub, ".tiff")),
    plot = p, width = width, height = height, units = "in",
    dpi = 600, compression = "lzw"
  )
}

theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# -----------------------------------------------------------------------------
# 1. Pooled and taxon-specific rank concordance
# -----------------------------------------------------------------------------

pairs <- data.table(
  analysis = c(
    "Phase2 global exact vs Ag3 strict 23bp exact",
    "Phase2 population q10 exact vs Ag3 strict 23bp exact",
    "Phase2 PAM intact vs Ag3 exact 3bp PAM",
    "Phase2 PAM intact vs Ag3 functional NGG PAM",
    "Phase2 max alt AF vs Ag3 max-position nonreference AF"
  ),
  x = c(
    "phase2_global_exact",
    "phase2_population_q10_exact",
    "phase2_pam_intact",
    "phase2_pam_intact",
    "phase2_max_alt_af"
  ),
  y = c(
    "strict_genotype_exact_23bp_fraction",
    "strict_genotype_exact_23bp_fraction",
    "pam_exact_3bp_fraction",
    "pam_ngg_intact_fraction",
    "max_position_nonreference_allele_fraction_23bp"
  ),
  interpretation = c(
    "Rank diagnostic only: phased haplotype vs unphased genotype estimands",
    "Rank diagnostic only: population-tail phased metric vs unphased genotype estimand",
    "Descriptive PAM-sequence concordance",
    "Descriptive functional-PAM concordance; Ag3 allows variation at PAM N",
    "Preferred quantitative cross-release allele-frequency concordance"
  )
)

corr_rows <- list()
for (i in seq_len(nrow(pairs))) {
  xx <- pairs$x[i]
  yy <- pairs$y[i]

  corr_rows[[length(corr_rows) + 1L]] <- data.table(
    scope = "pooled_taxon_rows",
    taxon = "pooled",
    analysis = pairs$analysis[i],
    interpretation = pairs$interpretation[i],
    n = sum(is.finite(d[[xx]]) & is.finite(d[[yy]])),
    spearman_rho = safe_spearman(d[[xx]], d[[yy]]),
    pearson_r = safe_pearson(d[[xx]], d[[yy]])
  )

  for (tx in c("gambiae", "coluzzii")) {
    z <- d[taxon == tx]
    corr_rows[[length(corr_rows) + 1L]] <- data.table(
      scope = "taxon_specific",
      taxon = tx,
      analysis = pairs$analysis[i],
      interpretation = pairs$interpretation[i],
      n = sum(is.finite(z[[xx]]) & is.finite(z[[yy]])),
      spearman_rho = safe_spearman(z[[xx]], z[[yy]]),
      pearson_r = safe_pearson(z[[xx]], z[[yy]])
    )
  }
}
corr_taxon <- rbindlist(corr_rows)
fwrite(corr_taxon, file.path(outdir, "Table_S1_rank_concordance_pooled_and_taxon.csv"))

# -----------------------------------------------------------------------------
# 2. Site-level mean across taxa (354 independent locked sites)
# -----------------------------------------------------------------------------

site <- d[, .(
  gene_id = first(gene_id),
  phase2_global_exact = first(phase2_global_exact),
  phase2_population_q10_exact = first(phase2_population_q10_exact),
  phase2_pam_intact = first(phase2_pam_intact),
  phase2_max_alt_af = first(phase2_max_alt_af),
  ag3_strict_exact_mean = mean(strict_genotype_exact_23bp_fraction, na.rm = TRUE),
  ag3_protospacer_exact_mean = mean(protospacer_exact_20bp_fraction, na.rm = TRUE),
  ag3_pam_exact_mean = mean(pam_exact_3bp_fraction, na.rm = TRUE),
  ag3_pam_ngg_mean = mean(pam_ngg_intact_fraction, na.rm = TRUE),
  ag3_functional_target_mean = mean(functional_target_intact_fraction, na.rm = TRUE),
  ag3_maxpos_nonref_af_mean = mean(max_position_nonreference_allele_fraction_23bp, na.rm = TRUE),
  callable_fraction_mean = mean(callable_sample_fraction, na.rm = TRUE)
), by = population_site_id]

site_corr <- rbindlist(list(
  data.table(
    analysis = "Phase2 global exact vs mean Ag3 strict exact",
    n = nrow(site),
    spearman_rho = safe_spearman(site$phase2_global_exact, site$ag3_strict_exact_mean),
    pearson_r = safe_pearson(site$phase2_global_exact, site$ag3_strict_exact_mean),
    interpretation = "Rank diagnostic only; estimands are not numerically interchangeable"
  ),
  data.table(
    analysis = "Phase2 population q10 exact vs mean Ag3 strict exact",
    n = nrow(site),
    spearman_rho = safe_spearman(site$phase2_population_q10_exact, site$ag3_strict_exact_mean),
    pearson_r = safe_pearson(site$phase2_population_q10_exact, site$ag3_strict_exact_mean),
    interpretation = "Rank diagnostic only; estimands are not numerically interchangeable"
  ),
  data.table(
    analysis = "Phase2 PAM intact vs mean Ag3 exact 3bp PAM",
    n = nrow(site),
    spearman_rho = safe_spearman(site$phase2_pam_intact, site$ag3_pam_exact_mean),
    pearson_r = safe_pearson(site$phase2_pam_intact, site$ag3_pam_exact_mean),
    interpretation = "Descriptive exact-PAM concordance"
  ),
  data.table(
    analysis = "Phase2 PAM intact vs mean Ag3 functional NGG PAM",
    n = nrow(site),
    spearman_rho = safe_spearman(site$phase2_pam_intact, site$ag3_pam_ngg_mean),
    pearson_r = safe_pearson(site$phase2_pam_intact, site$ag3_pam_ngg_mean),
    interpretation = "Descriptive functional-PAM concordance"
  ),
  data.table(
    analysis = "Phase2 max alt AF vs mean Ag3 max-position nonreference AF",
    n = nrow(site),
    spearman_rho = safe_spearman(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    pearson_r = safe_pearson(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    interpretation = "Preferred quantitative cross-release allele-frequency concordance"
  )
))
fwrite(site_corr, file.path(outdir, "Table_S2_site_level_concordance.csv"))

# -----------------------------------------------------------------------------
# 3. Cluster bootstrap CIs by locked site
# -----------------------------------------------------------------------------

boot_defs <- list(
  list(
    name = "Phase2 global exact vs Ag3 strict exact Spearman",
    fun = function(z) safe_spearman(
      z$phase2_global_exact,
      z$strict_genotype_exact_23bp_fraction
    )
  ),
  list(
    name = "Phase2 population q10 exact vs Ag3 strict exact Spearman",
    fun = function(z) safe_spearman(
      z$phase2_population_q10_exact,
      z$strict_genotype_exact_23bp_fraction
    )
  ),
  list(
    name = "Phase2 PAM intact vs Ag3 exact 3bp PAM Spearman",
    fun = function(z) safe_spearman(
      z$phase2_pam_intact,
      z$pam_exact_3bp_fraction
    )
  ),
  list(
    name = "Phase2 PAM intact vs Ag3 functional NGG PAM Spearman",
    fun = function(z) safe_spearman(
      z$phase2_pam_intact,
      z$pam_ngg_intact_fraction
    )
  ),
  list(
    name = "Phase2 max alt AF vs Ag3 max-position AF Spearman",
    fun = function(z) safe_spearman(
      z$phase2_max_alt_af,
      z$max_position_nonreference_allele_fraction_23bp
    )
  )
)

boot_rows <- lapply(seq_along(boot_defs), function(i) {
  log_msg(
    "Starting bootstrap ",
    i,
    "/",
    length(boot_defs),
    ": ",
    boot_defs[[i]]$name
  )

  ci <- boot_site_ci(
    d,
    stat_fun = boot_defs[[i]]$fun,
    B = 5000L,
    seed = 20260915L + i
  )

  log_msg(
    "Completed bootstrap ",
    i,
    "/",
    length(boot_defs)
  )

  data.table(
    analysis = boot_defs[[i]]$name,
    estimate = unname(ci["estimate"]),
    ci_low = unname(ci["ci_low"]),
    ci_high = unname(ci["ci_high"]),
    bootstrap_unit = "locked population_site_id; both taxon rows retained",
    B = 5000L
  )
})
boot_tab <- rbindlist(boot_rows)
fwrite(boot_tab, file.path(outdir, "Table_S3_site_cluster_bootstrap_CI.csv"))

# -----------------------------------------------------------------------------
# 4. Quantitative agreement for allele frequency metric only
# -----------------------------------------------------------------------------

af_agree <- data.table(
  metric = c(
    "Spearman rho",
    "Pearson r",
    "Mean difference (Ag3 - Phase2)",
    "Mean absolute difference",
    "RMSE",
    "Lin CCC"
  ),
  value = c(
    safe_spearman(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    safe_pearson(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    bias(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    mae(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    rmse(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean),
    ccc(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean)
  )
)
fwrite(af_agree, file.path(outdir, "Table_S4_allele_frequency_agreement.csv"))

af_lm <- lm(ag3_maxpos_nonref_af_mean ~ phase2_max_alt_af, data = site)
af_lm_tab <- data.table(
  term = rownames(coef(summary(af_lm))),
  estimate = coef(summary(af_lm))[, "Estimate"],
  std_error = coef(summary(af_lm))[, "Std. Error"],
  t_value = coef(summary(af_lm))[, "t value"],
  p_value = coef(summary(af_lm))[, "Pr(>|t|)"]
)
fwrite(af_lm_tab, file.path(outdir, "Table_S5_allele_frequency_calibration_regression.csv"))

# -----------------------------------------------------------------------------
# 5. Taxon divergence
# -----------------------------------------------------------------------------

wide <- dcast(
  d,
  population_site_id + gene_id + phase2_global_exact +
    phase2_population_q10_exact + phase2_pam_intact + phase2_max_alt_af ~ taxon,
  value.var = c(
    "strict_genotype_exact_23bp_fraction",
    "protospacer_exact_20bp_fraction",
    "pam_exact_3bp_fraction",
    "pam_ngg_intact_fraction",
    "functional_target_intact_fraction",
    "max_position_nonreference_allele_fraction_23bp"
  )
)

wide[, delta_strict_exact :=
       strict_genotype_exact_23bp_fraction_gambiae -
       strict_genotype_exact_23bp_fraction_coluzzii]
wide[, abs_delta_strict_exact := abs(delta_strict_exact)]

wide[, delta_functional_target :=
       functional_target_intact_fraction_gambiae -
       functional_target_intact_fraction_coluzzii]
wide[, abs_delta_functional_target := abs(delta_functional_target)]

wide[, delta_maxpos_af :=
       max_position_nonreference_allele_fraction_23bp_gambiae -
       max_position_nonreference_allele_fraction_23bp_coluzzii]
wide[, abs_delta_maxpos_af := abs(delta_maxpos_af)]

div_summary <- rbindlist(list(
  data.table(
    metric = "Strict 23bp exact fraction absolute taxon difference",
    mean = mean(wide$abs_delta_strict_exact, na.rm = TRUE),
    median = median(wide$abs_delta_strict_exact, na.rm = TRUE),
    q90 = unname(quantile(wide$abs_delta_strict_exact, 0.90, na.rm = TRUE)),
    q95 = unname(quantile(wide$abs_delta_strict_exact, 0.95, na.rm = TRUE)),
    max = max(wide$abs_delta_strict_exact, na.rm = TRUE)
  ),
  data.table(
    metric = "Functional target intact fraction absolute taxon difference",
    mean = mean(wide$abs_delta_functional_target, na.rm = TRUE),
    median = median(wide$abs_delta_functional_target, na.rm = TRUE),
    q90 = unname(quantile(wide$abs_delta_functional_target, 0.90, na.rm = TRUE)),
    q95 = unname(quantile(wide$abs_delta_functional_target, 0.95, na.rm = TRUE)),
    max = max(wide$abs_delta_functional_target, na.rm = TRUE)
  ),
  data.table(
    metric = "Maximum-position nonreference AF absolute taxon difference",
    mean = mean(wide$abs_delta_maxpos_af, na.rm = TRUE),
    median = median(wide$abs_delta_maxpos_af, na.rm = TRUE),
    q90 = unname(quantile(wide$abs_delta_maxpos_af, 0.90, na.rm = TRUE)),
    q95 = unname(quantile(wide$abs_delta_maxpos_af, 0.95, na.rm = TRUE)),
    max = max(wide$abs_delta_maxpos_af, na.rm = TRUE)
  )
))
fwrite(div_summary, file.path(outdir, "Table_S6_taxon_divergence_summary.csv"))

top_divergent <- copy(wide)
setorder(top_divergent, -abs_delta_strict_exact, -abs_delta_functional_target)
fwrite(
  top_divergent[1:min(30L, .N)],
  file.path(outdir, "Table_S7_top30_taxon_divergent_sites.csv")
)

# -----------------------------------------------------------------------------
# 6. Mechanistic conservation descriptors
# -----------------------------------------------------------------------------
# These are descriptive continuous summaries, not new discovery thresholds.

mechanistic <- d[, .(
  n_sites = .N,
  mean_strict_exact_23bp = mean(strict_genotype_exact_23bp_fraction, na.rm = TRUE),
  median_strict_exact_23bp = median(strict_genotype_exact_23bp_fraction, na.rm = TRUE),
  mean_protospacer_exact_20bp = mean(protospacer_exact_20bp_fraction, na.rm = TRUE),
  median_protospacer_exact_20bp = median(protospacer_exact_20bp_fraction, na.rm = TRUE),
  mean_exact_pam_3bp = mean(pam_exact_3bp_fraction, na.rm = TRUE),
  mean_functional_ngg_pam = mean(pam_ngg_intact_fraction, na.rm = TRUE),
  mean_functional_target_intact = mean(functional_target_intact_fraction, na.rm = TRUE),
  mean_max_position_nonreference_af = mean(max_position_nonreference_allele_fraction_23bp, na.rm = TRUE)
), by = taxon]

fwrite(mechanistic, file.path(outdir, "Table_S8_sequence_aware_mechanistic_summary_by_taxon.csv"))

# Re-use the prespecified status from Step 06a instead of inventing new thresholds.
if ("external_validation_status" %in% names(site06a)) {
  status_tab <- site06a[, .N, by = external_validation_status][order(-N)]
  status_tab[, proportion := N / sum(N)]
  fwrite(status_tab, file.path(outdir, "Table_1_external_validation_status_counts.csv"))
}

# -----------------------------------------------------------------------------
# 7. Most discordant sites relative to Phase-2 ordering
# -----------------------------------------------------------------------------
# Because exactness estimands differ, "discordance" is defined by rank residual,
# not raw numerical difference.

site[, phase2_exact_rank := frank(-phase2_global_exact, ties.method = "average")]
site[, ag3_exact_rank := frank(-ag3_strict_exact_mean, ties.method = "average")]
site[, exact_rank_shift := ag3_exact_rank - phase2_exact_rank]
site[, abs_exact_rank_shift := abs(exact_rank_shift)]

site[, phase2_af_rank := frank(phase2_max_alt_af, ties.method = "average")]
site[, ag3_af_rank := frank(ag3_maxpos_nonref_af_mean, ties.method = "average")]
site[, af_rank_shift := ag3_af_rank - phase2_af_rank]
site[, abs_af_rank_shift := abs(af_rank_shift)]

top_rank_shift <- copy(site)
setorder(top_rank_shift, -abs_exact_rank_shift, -abs_af_rank_shift)
fwrite(
  top_rank_shift[1:min(30L, .N)],
  file.path(outdir, "Table_S9_top30_cross_release_rank_shifts.csv")
)

# -----------------------------------------------------------------------------
# 8. Publication-quality figures
# -----------------------------------------------------------------------------

rho_exact <- safe_spearman(site$phase2_global_exact, site$ag3_strict_exact_mean)
p1 <- ggplot(
  site,
  aes(x = phase2_global_exact, y = ag3_strict_exact_mean)
) +
  geom_point(alpha = 0.65, size = 1.8) +
  labs(
    title = "Cross-release rank concordance of target-sequence conservation",
    subtitle = sprintf(
      "354 locked sites; Spearman rho = %.3f. Metrics are not numerically interchangeable.",
      rho_exact
    ),
    x = "Phase-2 phased-haplotype exact-match fraction",
    y = "Ag3 mean strict-genotype exact 23-bp fraction"
  ) +
  theme_pub
save_plot(p1, "Figure_1_phase2_vs_ag3_exact_rank_concordance")

rho_af <- safe_spearman(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean)
p2 <- ggplot(
  site,
  aes(x = phase2_max_alt_af, y = ag3_maxpos_nonref_af_mean)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.65, size = 1.8) +
  labs(
    title = "Cross-release allele-frequency concordance",
    subtitle = sprintf(
      "354 locked sites; Spearman rho = %.3f; dashed line = identity",
      rho_af
    ),
    x = "Phase-2 maximum target-variant alternate-allele frequency",
    y = "Ag3 mean maximum-position nonreference allele frequency"
  ) +
  theme_pub
save_plot(p2, "Figure_2_phase2_vs_ag3_max_position_AF")

p3 <- ggplot(
  wide,
  aes(
    x = strict_genotype_exact_23bp_fraction_gambiae,
    y = strict_genotype_exact_23bp_fraction_coluzzii
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.65, size = 1.8) +
  coord_equal() +
  labs(
    title = "Taxon concordance of Ag3 strict 23-bp target conservation",
    x = "An. gambiae strict-genotype exact fraction",
    y = "An. coluzzii strict-genotype exact fraction"
  ) +
  theme_pub
save_plot(p3, "Figure_3_ag3_gambiae_vs_coluzzii_strict_exact")

long_fun <- melt(
  d[, .(
    population_site_id,
    taxon,
    strict_23bp = strict_genotype_exact_23bp_fraction,
    protospacer_20bp = protospacer_exact_20bp_fraction,
    functional_ngg_pam = pam_ngg_intact_fraction,
    functional_target = functional_target_intact_fraction
  )],
  id.vars = c("population_site_id", "taxon"),
  variable.name = "metric",
  value.name = "fraction"
)

metric_labels <- c(
  strict_23bp = "Strict 23-bp exact",
  protospacer_20bp = "20-bp protospacer exact",
  functional_ngg_pam = "Functional NGG PAM",
  functional_target = "Functional target intact"
)
long_fun[, metric := factor(metric, levels = names(metric_labels), labels = metric_labels)]

p4 <- ggplot(long_fun, aes(x = fraction, colour = taxon)) +
  stat_ecdf(linewidth = 0.9) +
  facet_wrap(~ metric, ncol = 2) +
  labs(
    title = "Sequence-aware conservation of locked targets in Ag3",
    x = "Fraction of samples",
    y = "Empirical cumulative probability",
    colour = "Taxon"
  ) +
  theme_pub
save_plot(p4, "Figure_4_sequence_aware_conservation_ECDF", width = 8.2, height = 6.4)

if ("external_validation_status" %in% names(site06a)) {
  status_plot <- site06a[, .N, by = external_validation_status]
  p5 <- ggplot(
    status_plot,
    aes(x = reorder(external_validation_status, N), y = N)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Prespecified external-validation status of locked Phase-2 sites",
      x = NULL,
      y = "Number of locked sites"
    ) +
    theme_pub
  save_plot(p5, "Figure_5_external_validation_status_counts", width = 7.4, height = 4.8)
}

top20 <- top_rank_shift[1:min(20L, .N)]
top20[, population_site_id := factor(
  population_site_id,
  levels = rev(population_site_id)
)]
p6 <- ggplot(top20, aes(x = population_site_id, y = exact_rank_shift)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Largest cross-release shifts in target-conservation rank",
    subtitle = "Positive values indicate a less favourable Ag3 rank than Phase-2",
    x = "Locked target site",
    y = "Ag3 rank - Phase-2 rank"
  ) +
  theme_pub
save_plot(p6, "Figure_6_top20_exact_rank_shifts", width = 8.2, height = 6.2)

# -----------------------------------------------------------------------------
# 9. Compact manuscript-ready QC summary
# -----------------------------------------------------------------------------

qc <- data.table(
  metric = c(
    "Taxon-level rows",
    "Unique locked sites",
    "Unique genes",
    "Taxa",
    "Mean callable sample fraction",
    "Phase2 global exact vs Ag3 strict exact Spearman rho",
    "Phase2 q10 exact vs Ag3 strict exact Spearman rho",
    "Phase2 max AF vs Ag3 max-position AF Spearman rho",
    "Phase2 max AF vs Ag3 max-position AF MAE",
    "Phase2 max AF vs Ag3 max-position AF RMSE",
    "Discovery ranks changed by Ag3",
    "Discovery thresholds changed by Ag3"
  ),
  value = c(
    nrow(d),
    uniqueN(d$population_site_id),
    uniqueN(d$gene_id),
    paste(sort(unique(d$taxon)), collapse = ","),
    signif(mean(d$callable_sample_fraction, na.rm = TRUE), 6),
    signif(safe_spearman(d$phase2_global_exact, d$strict_genotype_exact_23bp_fraction), 6),
    signif(safe_spearman(d$phase2_population_q10_exact, d$strict_genotype_exact_23bp_fraction), 6),
    signif(safe_spearman(d$phase2_max_alt_af, d$max_position_nonreference_allele_fraction_23bp), 6),
    signif(mae(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean), 6),
    signif(rmse(site$phase2_max_alt_af, site$ag3_maxpos_nonref_af_mean), 6),
    0,
    0
  )
)
fwrite(qc, file.path(outdir, "07_ag3_manuscript_validation_qc.csv"))

# Save site-level derived dataset for reproducibility.
fwrite(site, file.path(outdir, "07_site_level_mean_across_taxa.csv"))
fwrite(wide, file.path(outdir, "07_site_level_taxon_wide.csv"))

# Session information.
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))

log_msg(
  "Manuscript-grade Ag3 validation completed. Outputs written to ",
  outdir
)

print(qc)
