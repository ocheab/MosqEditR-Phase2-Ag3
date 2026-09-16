# STEP 07: generate publication-quality figures and manuscript/supplement tables

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 07 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "ggplot2", "scales", "yaml"))

STEP <- "07"
log_step(STEP, "Building publication figures and tables")

required <- c(
  "data_processed/03_population_site_metrics.csv",
  "data_processed/03_site_population_metrics.csv",
  "data_processed/04_best_multiplex_portfolios.csv",
  "data_processed/05_gene_population_robustness.csv",
  "data_processed/06_site_geographic_heterogeneity.csv",
  "data_processed/06_rank_validation_summary.csv",
  "data_processed/06a_ag3_external_site_validation.csv",
  "data_processed/06a_ag3_external_gene_support.csv"
)
assert_files(required, nonempty = TRUE)

cfg <- read_config()
png_dpi <- as.integer(cfg_get(cfg, "plotting", "png_dpi", default = 400L))
tiff_dpi <- as.integer(cfg_get(cfg, "plotting", "tiff_dpi", default = 600L))
if (!is.finite(png_dpi) || png_dpi < 72L) png_dpi <- 400L
if (!is.finite(tiff_dpi) || tiff_dpi < 150L) tiff_dpi <- 600L

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("supplement", recursive = TRUE, showWarnings = FALSE)

site <- data.table::fread(required[[1]], na.strings = c("", "NA", "NaN"))
site_pop <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
portfolio <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))
gene <- data.table::fread(required[[4]], na.strings = c("", "NA", "NaN"))
hetero <- data.table::fread(required[[5]], na.strings = c("", "NA", "NaN"))
rank_summary <- data.table::fread(required[[6]], na.strings = c("", "NA", "NaN"))
ext_site <- data.table::fread(required[[7]], na.strings = c("", "NA", "NaN"))
ext_gene <- data.table::fread(required[[8]], na.strings = c("", "NA", "NaN"))

assert_columns(site, c("population_site_id", "gene_id", "population_target_site_robustness_score", "population_robustness_class", "final_prepopulation_rank"), "site metrics")
assert_columns(site_pop, c("population_site_id", "population_id", "target_23bp_exact_match_fraction", "primary_population_n20"), "site-population metrics")
assert_columns(portfolio, c("gene_id", "portfolio_size", "portfolio_site_ids", "multiplex_robustness_score", "portfolio_rank_within_gene_size"), "best portfolios")
assert_columns(gene, c("gene_id", "is_external_benchmark", "final_prepopulation_rank", "population_filtered_prepopulation_rank", "gene_population_robustness_score", "best_single_site_score"), "gene metrics")
assert_columns(hetero, c("population_site_id", "population_exact_mean", "population_exact_range"), "heterogeneity metrics")
assert_columns(ext_site, c("population_site_id", "phase2_discovery_exact_23bp", "exact_23bp_fraction", "external_validation_status", "external_stable_observed"), "Ag3 external site validation")
assert_columns(ext_gene, c("gene_id", "best_single_ag3_confirmed", "ag3_external_validation_used_for_reranking"), "Ag3 external gene support")

gene[, is_external_benchmark := normalize_flag(is_external_benchmark)]
if (anyNA(gene$is_external_benchmark)) fail_step(STEP, "is_external_benchmark contains unparseable values in gene metrics.")

fig_manifest <- list()
fig_i <- 0L

theme_pub <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
    plot.title = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "grey95")
  )

placeholder_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
}

save_plot <- function(p, stem, w = 8, h = 5.5, note = "") {
  pdf_path <- file.path("figures", paste0(stem, ".pdf"))
  png_path <- file.path("figures", paste0(stem, ".png"))
  tif_path <- file.path("figures", paste0(stem, ".tif"))

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  tryCatch(
    ggplot2::ggsave(pdf_path, plot = p, width = w, height = h, units = "in", device = pdf_device, limitsize = FALSE),
    error = function(e) {
      log_warning(STEP, paste0("Cairo/PDF save failed for ", stem, "; retrying with base PDF. Error: ", conditionMessage(e)))
      ggplot2::ggsave(pdf_path, plot = p, width = w, height = h, units = "in", device = "pdf", limitsize = FALSE)
    }
  )
  ggplot2::ggsave(png_path, plot = p, width = w, height = h, units = "in", dpi = png_dpi, limitsize = FALSE)
  tif_ok <- FALSE
  if (isTRUE(capabilities("tiff"))) {
    tif_ok <- tryCatch({
      ggplot2::ggsave(tif_path, plot = p, width = w, height = h, units = "in", dpi = tiff_dpi, device = "tiff", compression = "lzw", limitsize = FALSE)
      TRUE
    }, error = function(e1) {
      log_warning(STEP, paste0("TIFF save with LZW failed for ", stem, "; retrying without compression. Error: ", conditionMessage(e1)))
      tryCatch({
        ggplot2::ggsave(tif_path, plot = p, width = w, height = h, units = "in", dpi = tiff_dpi, device = "tiff", limitsize = FALSE)
        TRUE
      }, error = function(e2) {
        log_warning(STEP, paste0("TIFF output is unavailable for ", stem, ". PDF and PNG were created. Error: ", conditionMessage(e2)))
        FALSE
      })
    })
  } else {
    log_warning(STEP, paste0("This R build has no TIFF capability; skipping TIFF for ", stem, ". PDF and PNG are still required."))
  }

  for (path in c(pdf_path, png_path)) {
    if (!file.exists(path) || file.info(path)$size <= 0) fail_step(STEP, paste0("Required figure output was not created correctly: ", path))
  }
  if (tif_ok && (!file.exists(tif_path) || file.info(tif_path)$size <= 0)) {
    log_warning(STEP, paste0("TIFF writer returned without a valid file for ", stem, "; recording TIFF as unavailable."))
    tif_ok <- FALSE
  }
  if (!tif_ok && file.exists(tif_path)) unlink(tif_path, force = TRUE)

  fig_i <<- fig_i + 1L
  fig_manifest[[fig_i]] <<- data.table::data.table(
    figure_stem = stem,
    width_in = w,
    height_in = h,
    pdf = pdf_path,
    png = png_path,
    tif = if (tif_ok) tif_path else NA_character_,
    note = note
  )
  invisible(c(pdf_path, png_path, if (tif_ok) tif_path else character()))
}

# Figure 2: distribution of site robustness scores by classification.
f2 <- site[is.finite(population_target_site_robustness_score)]
if (nrow(f2)) {
  p2 <- ggplot2::ggplot(f2, ggplot2::aes(x = population_target_site_robustness_score)) +
    ggplot2::geom_histogram(bins = 30, boundary = 0) +
    ggplot2::facet_wrap(~population_robustness_class, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Population target-site robustness score",
      y = "Number of target sites"
    ) + theme_pub
  save_plot(p2, "Fig2_site_robustness_distribution")
} else {
  save_plot(placeholder_plot("No finite site robustness scores were available."), "Fig2_site_robustness_distribution", note = "Placeholder: no finite site robustness scores")
}

# Figure 3: population heatmap for the best assessable site of top-ranked novel genes.
site_for_best <- site[
  is.finite(population_target_site_robustness_score) &
    population_robustness_class %in% c("ROBUST", "INTERMEDIATE", "CONCERN")
]
data.table::setorder(site_for_best, gene_id, -population_target_site_robustness_score, population_site_rank_within_gene, population_site_id)
best_site <- if (nrow(site_for_best)) site_for_best[, .SD[1L], by = gene_id] else site_for_best
novel_gene <- gene[is_external_benchmark == FALSE][order(population_filtered_prepopulation_rank)]
topgenes <- utils::head(novel_gene$gene_id, 30L)
heat_ids <- best_site[gene_id %in% topgenes, population_site_id]
h <- site_pop[population_site_id %in% heat_ids & primary_population_n20 == TRUE]
if (nrow(h)) {
  h <- merge(h, best_site[, .(population_site_id, gene_id)], by = "population_site_id", all.x = TRUE, sort = FALSE)
  plotted_genes <- topgenes[topgenes %in% h$gene_id]
  h[, gene_id := factor(gene_id, levels = rev(plotted_genes))]
  pop_order <- unique(h$population_id)
  h[, population_id := factor(population_id, levels = pop_order)]
  p3 <- ggplot2::ggplot(h, ggplot2::aes(x = population_id, y = gene_id, fill = target_23bp_exact_match_fraction)) +
    ggplot2::geom_tile(na.rm = FALSE) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), labels = scales::percent_format(accuracy = 1), na.value = "grey90") +
    ggplot2::labs(x = "Country × taxon population", y = "Gene", fill = "23-bp exact") +
    theme_pub +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1), legend.position = "right")
  save_plot(p3, "Fig3_population_exact_match_heatmap", 10, 8)
} else {
  save_plot(placeholder_plot("No n>=20 site-by-population heatmap data were available."), "Fig3_population_exact_match_heatmap", 10, 8, note = "Placeholder: no heatmap data")
}

# Figure 4: single-site robustness versus best two-site portfolio robustness.
pair <- portfolio[portfolio_size == 2L]
x4 <- merge(
  gene[, .(gene_id, is_external_benchmark, best_single_site_score)],
  pair[, .(gene_id, best_pair_score = multiplex_robustness_score)],
  by = "gene_id", all = FALSE
)
x4 <- x4[is.finite(best_single_site_score) & is.finite(best_pair_score)]
if (nrow(x4)) {
  p4 <- ggplot2::ggplot(x4, ggplot2::aes(x = best_single_site_score, y = best_pair_score, shape = is_external_benchmark)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Best single-site robustness",
      y = "Best two-site portfolio robustness",
      shape = "External benchmark"
    ) + theme_pub
  save_plot(p4, "Fig4_multiplex_rescue")
} else {
  save_plot(placeholder_plot("No genes had both finite single-site and pair robustness scores."), "Fig4_multiplex_rescue", note = "Placeholder: no paired robustness data")
}

# Figure 5: rank-shift dumbbell plot for the top 50 frozen pre-population genes.
r <- gene[final_prepopulation_rank <= 50L][order(final_prepopulation_rank)]
r <- r[is.finite(final_prepopulation_rank) & is.finite(population_filtered_prepopulation_rank)]
if (nrow(r)) {
  r[, gene_factor := factor(gene_id, levels = rev(gene_id))]
  pts <- data.table::rbindlist(list(
    r[, .(gene_factor, rank = final_prepopulation_rank, ranking = "Frozen pre-population")],
    r[, .(gene_factor, rank = population_filtered_prepopulation_rank, ranking = "Population-filtered")]
  ))
  p5 <- ggplot2::ggplot(r, ggplot2::aes(y = gene_factor)) +
    ggplot2::geom_segment(ggplot2::aes(x = final_prepopulation_rank, xend = population_filtered_prepopulation_rank, yend = gene_factor), linewidth = 0.4) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = rank, y = gene_factor, shape = ranking), size = 1.8, inherit.aes = FALSE) +
    ggplot2::labs(x = "Rank", y = "Gene", shape = NULL) +
    theme_pub +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  save_plot(p5, "Fig5_population_rank_shift_top50", 9, 10)
} else {
  save_plot(placeholder_plot("No top-50 rank-shift data were available."), "Fig5_population_rank_shift_top50", 9, 10, note = "Placeholder: no rank-shift data")
}

# Figure 6: external benchmark population robustness.
b <- gene[is_external_benchmark == TRUE]
if (nrow(b)) {
  if (!"benchmark_name" %in% names(b)) b[, benchmark_name := NA_character_]
  b[, benchmark_label := ifelse(nonempty_string(benchmark_name), benchmark_name, gene_id)]
  b <- b[is.finite(gene_population_robustness_score)]
}
if (nrow(b)) {
  p6 <- ggplot2::ggplot(b, ggplot2::aes(x = reorder(benchmark_label, gene_population_robustness_score), y = gene_population_robustness_score)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Gene-level population robustness score") + theme_pub
  save_plot(p6, "Fig6_benchmark_population_robustness", 8, 5)
} else {
  save_plot(placeholder_plot("No external benchmarks had finite population robustness scores."), "Fig6_benchmark_population_robustness", 8, 5, note = "Placeholder: no finite benchmark scores")
}

# Figure 7: geographic heterogeneity of exact target match.
x7 <- merge(
  site[, .(population_site_id, gene_id, final_prepopulation_rank, population_robustness_class)],
  hetero,
  by = "population_site_id", all.x = TRUE, sort = FALSE
)
x7 <- x7[is.finite(population_exact_mean) & is.finite(population_exact_range)]
if (nrow(x7)) {
  p7 <- ggplot2::ggplot(x7, ggplot2::aes(x = population_exact_mean, y = population_exact_range, shape = population_robustness_class)) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Mean population exact-match fraction",
      y = "Range across populations",
      shape = "Site class"
    ) + theme_pub
  save_plot(p7, "Fig7_geographic_heterogeneity")
} else {
  save_plot(placeholder_plot("No finite geographic heterogeneity metrics were available."), "Fig7_geographic_heterogeneity", note = "Placeholder: no heterogeneity data")
}


# Figure 8: independent Ag3 external validation of frozen Phase-2 discoveries.
x8 <- ext_site[is.finite(phase2_discovery_exact_23bp) & is.finite(exact_23bp_fraction)]
if (nrow(x8)) {
  p8 <- ggplot2::ggplot(x8, ggplot2::aes(x = phase2_discovery_exact_23bp, y = exact_23bp_fraction, shape = external_validation_status)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.80, size = 2) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Phase-2 discovery: 23-bp exact-match fraction",
      y = "Ag3 external validation: 23-bp exact-match fraction",
      shape = "External validation"
    ) + theme_pub
  save_plot(p8, "Fig8_phase2_ag3_external_validation", 8, 6)
} else {
  save_plot(placeholder_plot("No paired Phase-2/Ag3 external validation metrics were available."), "Fig8_phase2_ag3_external_validation", 8, 6, note = "Placeholder: no external-validation concordance data")
}

# Main manuscript tables.
novel_sorted <- gene[is_external_benchmark == FALSE][order(population_filtered_prepopulation_rank)]
atomic_fwrite(utils::head(novel_sorted, 20L), "tables/Table2_top20_population_validated_candidates.csv")
atomic_fwrite(gene[is_external_benchmark == TRUE][order(population_filtered_prepopulation_rank)], "tables/Table3_benchmark_population_validation.csv")
atomic_fwrite(portfolio[portfolio_rank_within_gene_size == 1L][order(gene_id, portfolio_size)], "tables/Table4_best_multiplex_portfolios.csv")

# Compact ranking/sensitivity summary table for manuscript methods/results cross-checking.
atomic_fwrite(rank_summary, "tables/Table5_rank_and_sensitivity_summary.csv")
atomic_fwrite(ext_site, "tables/Table6_ag3_external_site_validation.csv")

# Full supplementary tables.
atomic_fwrite(site, "supplement/TableS1_all_site_population_metrics.csv")
atomic_fwrite(site_pop, "supplement/TableS2_site_by_population_metrics.csv")
atomic_fwrite(gene, "supplement/TableS3_all_gene_population_robustness.csv")
atomic_fwrite(portfolio, "supplement/TableS4_best_gene_portfolios.csv")
atomic_fwrite(ext_gene, "supplement/TableS5_ag3_external_gene_support.csv")

fig_manifest_dt <- data.table::rbindlist(fig_manifest, fill = TRUE, use.names = TRUE)
atomic_fwrite(fig_manifest_dt, "tables/Figure_output_manifest.csv")
write_checksum(
  c(
    unlist(fig_manifest_dt[, .(pdf, png, tif)]),
    "tables/Table2_top20_population_validated_candidates.csv",
    "tables/Table3_benchmark_population_validation.csv",
    "tables/Table4_best_multiplex_portfolios.csv",
    "tables/Table5_rank_and_sensitivity_summary.csv",
    "tables/Table6_ag3_external_site_validation.csv",
    "supplement/TableS1_all_site_population_metrics.csv",
    "supplement/TableS2_site_by_population_metrics.csv",
    "supplement/TableS3_all_gene_population_robustness.csv",
    "supplement/TableS4_best_gene_portfolios.csv",
    "supplement/TableS5_ag3_external_gene_support.csv"
  ),
  "logs/07_checksums.tsv"
)
write_session_info("logs/07_sessionInfo.txt")
print(fig_manifest_dt)
log_step(STEP, "Publication figures and tables completed successfully")
