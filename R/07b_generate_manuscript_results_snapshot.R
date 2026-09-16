# STEP 07b: generate a compact manuscript-results snapshot from completed real outputs

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 07b from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages("data.table")

STEP <- "07b"
log_step(STEP, "Generating manuscript results snapshot")

required <- c(
  "data_processed/02_analysis_sample_metadata.csv",
  "data_processed/02_analysis_population_manifest.csv",
  "data_processed/03_population_site_metrics.csv",
  "data_processed/04_multiplex_portfolio_metrics.csv",
  "data_processed/05_gene_population_robustness.csv",
  "data_processed/06_rank_validation_summary.csv",
  "data_processed/06a_ag3_external_validation_qc.csv"
)
assert_files(required, nonempty = TRUE)

meta <- data.table::fread(required[[1]], na.strings = c("", "NA", "NaN"))
pops <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
site <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))
portfolio <- data.table::fread(required[[4]], na.strings = c("", "NA", "NaN"))
gene <- data.table::fread(required[[5]], na.strings = c("", "NA", "NaN"))
rankv <- data.table::fread(required[[6]], na.strings = c("", "NA", "NaN"))
extv <- data.table::fread(required[[7]], na.strings = c("", "NA", "NaN"))

assert_columns(meta, c("sample_id", "primary_taxon"), "analysis sample metadata")
assert_columns(pops, c("population_id", "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30"), "analysis population manifest")
assert_columns(site, c("population_robustness_class", "population_target_site_robustness_score", "target_23bp_exact_match_fraction", "population_q10_target_23bp_exact"), "site metrics")
assert_columns(gene, c("gene_id", "population_deployment_class", "gene_population_robustness_score", "population_filtered_prepopulation_rank", "is_external_benchmark"), "gene metrics")
assert_columns(rankv, c("metric", "value"), "rank validation summary")

meta[, primary_taxon := normalize_flag(primary_taxon)]
primary_meta <- meta[primary_taxon == TRUE]
site_class <- site[, .N, by = population_robustness_class][order(population_robustness_class)]
gene_class <- gene[, .N, by = population_deployment_class][order(population_deployment_class)]
finite_site <- site[is.finite(population_target_site_robustness_score)]
finite_gene <- gene[is.finite(gene_population_robustness_score)]

core <- data.table::data.table(
  section = "core",
  metric = c(
    "Primary gambiae/coluzzii samples in phased panel",
    "Primary populations n>=20",
    "Sensitivity populations n>=10",
    "Sensitivity populations n>=30",
    "Assessable sites with finite robustness score",
    "Median site robustness score",
    "Median global 23-bp exact-match fraction",
    "Median population q10 23-bp exact-match fraction",
    "Evaluated multiplex portfolios",
    "Genes with finite population robustness score",
    "Median gene population robustness score"
  ),
  value = as.character(c(
    data.table::uniqueN(primary_meta$sample_id),
    sum(pops$primary_population_n20 == TRUE, na.rm = TRUE),
    sum(pops$sensitivity_population_n10 == TRUE, na.rm = TRUE),
    sum(pops$sensitivity_population_n30 == TRUE, na.rm = TRUE),
    nrow(finite_site),
    safe_median(finite_site$population_target_site_robustness_score),
    safe_median(finite_site$target_23bp_exact_match_fraction),
    safe_median(finite_site$population_q10_target_23bp_exact),
    nrow(portfolio),
    nrow(finite_gene),
    safe_median(finite_gene$gene_population_robustness_score)
  ))
)

snapshot <- data.table::rbindlist(list(
  core,
  site_class[, .(section = "site_class", metric = population_robustness_class, value = as.character(N))],
  gene_class[, .(section = "gene_class", metric = population_deployment_class, value = as.character(N))],
  rankv[, .(section = "rank_validation", metric = as.character(metric), value = as.character(value))],
  extv[, .(section = "ag3_external_validation", metric = as.character(metric), value = as.character(value))]
), fill = TRUE, use.names = TRUE)

# Add top-ranked novel candidates and benchmarks as explicit manuscript-ready rows.
gene[, is_external_benchmark := normalize_flag(is_external_benchmark)]
novel <- gene[is_external_benchmark == FALSE][order(population_filtered_prepopulation_rank)]
if (nrow(novel)) {
  top_n <- utils::head(novel, 10L)
  snapshot <- data.table::rbindlist(list(
    snapshot,
    top_n[, .(
      section = "top_novel_gene",
      metric = paste0("rank_", population_filtered_prepopulation_rank, "_", gene_id),
      value = paste0("score=", signif(gene_population_robustness_score, 6), ";class=", population_deployment_class)
    )]
  ), fill = TRUE)
}

bench <- gene[is_external_benchmark == TRUE][order(population_filtered_prepopulation_rank)]
if (nrow(bench)) {
  if (!"benchmark_name" %in% names(bench)) bench[, benchmark_name := NA_character_]
  snapshot <- data.table::rbindlist(list(
    snapshot,
    bench[, .(
      section = "benchmark_gene",
      metric = ifelse(nonempty_string(benchmark_name), benchmark_name, gene_id),
      value = paste0("gene_id=", gene_id, ";score=", signif(gene_population_robustness_score, 6), ";class=", population_deployment_class)
    )]
  ), fill = TRUE)
}

atomic_fwrite(snapshot, "tables/Manuscript2_results_snapshot.csv")

lines <- c(
  "MOSQEDIT-R MANUSCRIPT 2 — REAL RESULTS SNAPSHOT",
  "Generated from Phase-2 phased discovery plus locked Ag3 site-level external-validation outputs. Do not edit values manually.",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste(snapshot$section, snapshot$metric, snapshot$value, sep = "\t")
)
atomic_writeLines(lines, "manuscript/MANUSCRIPT_RESULT_VALUES.txt")
write_checksum(
  c("tables/Manuscript2_results_snapshot.csv", "manuscript/MANUSCRIPT_RESULT_VALUES.txt"),
  "logs/07b_checksums.tsv"
)
write_session_info("logs/07b_sessionInfo.txt")
print(snapshot)
log_step(STEP, "Manuscript results snapshot completed successfully")
