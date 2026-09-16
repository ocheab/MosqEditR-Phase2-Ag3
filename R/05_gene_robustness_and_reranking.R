# STEP 05: gene-level robustness, deployment classes, and population-validated reranking

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 05 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "05"
log_step(STEP, "Building gene-level population robustness and reranking")

required <- c(
  "data_processed/00_frozen_gene_panel.csv",
  "data_processed/03_population_site_metrics.csv",
  "data_processed/04_best_multiplex_portfolios.csv"
)
assert_files(required, nonempty = TRUE)

cfg <- read_config()
gene_weights <- unlist(cfg_get(cfg, "gene_robustness_weights", required = TRUE), use.names = TRUE)
validate_weights(gene_weights, "gene_robustness_weights")
needed_w <- c("best_single_site", "best_pair", "best_triple")
if (!all(needed_w %in% names(gene_weights))) fail_step(STEP, paste0("gene_robustness_weights must define: ", paste(needed_w, collapse = ", ")))
gene_weights <- gene_weights[needed_w]

panel <- data.table::fread(required[[1]], na.strings = c("", "NA", "NaN"))
site <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
best <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))

assert_columns(panel, c("gene_id", "final_prepopulation_rank", "is_external_benchmark"), "frozen gene panel")
assert_unique(panel, "gene_id", "frozen gene panel")
assert_columns(site, c(
  "gene_id", "population_site_id", "population_site_rank_within_gene",
  "population_target_site_robustness_score", "population_robustness_class",
  "target_23bp_exact_match_fraction", "population_q10_target_23bp_exact",
  "pam_intact_fraction", "max_target_variant_alt_af"
), "population site metrics")
assert_unique(site, "population_site_id", "population site metrics")

# Step 04 may legitimately produce zero portfolios; fread will still retain headers from the robust script.
expected_best_cols <- c(
  "gene_id", "portfolio_size", "portfolio_site_ids", "portfolio_span_bp",
  "global_at_least_one_intact_fraction", "population_q10_at_least_one_intact",
  "population_breadth_ge_0_95", "multiplex_robustness_score"
)
assert_columns(best, expected_best_cols, "best multiplex portfolios")

panel[, is_external_benchmark := normalize_flag(is_external_benchmark)]
if (anyNA(panel$is_external_benchmark)) fail_step(STEP, "is_external_benchmark contains an unparseable value.")

# Select ONLY the needed site columns before merge. This avoids suffixing or silently
# replacing gene-level fields such as final_prepopulation_rank.
single <- site[
  is.finite(population_target_site_robustness_score) &
    population_robustness_class %in% c("ROBUST", "INTERMEDIATE", "CONCERN"),
  .(
    gene_id,
    population_site_id,
    population_site_rank_within_gene,
    population_target_site_robustness_score,
    population_robustness_class,
    target_23bp_exact_match_fraction,
    population_q10_target_23bp_exact,
    pam_intact_fraction,
    max_target_variant_alt_af
  )
]
data.table::setorder(single, gene_id, -population_target_site_robustness_score, population_site_rank_within_gene, population_site_id)

if (nrow(single)) {
  best_single <- single[, .SD[1L], by = gene_id]
  data.table::setnames(
    best_single,
    old = c(
      "population_site_id", "population_target_site_robustness_score", "population_robustness_class",
      "target_23bp_exact_match_fraction", "population_q10_target_23bp_exact",
      "pam_intact_fraction", "max_target_variant_alt_af"
    ),
    new = c(
      "best_single_site_id", "best_single_site_score", "best_single_site_class",
      "best_single_global_exact", "best_single_population_q10_exact",
      "best_single_pam_intact", "best_single_max_alt_af"
    )
  )
} else {
  best_single <- data.table::data.table(
    gene_id = character(), best_single_site_id = character(), population_site_rank_within_gene = integer(),
    best_single_site_score = numeric(), best_single_site_class = character(),
    best_single_global_exact = numeric(), best_single_population_q10_exact = numeric(),
    best_single_pam_intact = numeric(), best_single_max_alt_af = numeric()
  )
}

pair <- best[portfolio_size == 2L, .(
  gene_id,
  best_pair_sites = portfolio_site_ids,
  best_pair_score = multiplex_robustness_score,
  best_pair_global_atleast1 = global_at_least_one_intact_fraction,
  best_pair_q10_atleast1 = population_q10_at_least_one_intact,
  best_pair_breadth_ge95 = population_breadth_ge_0_95,
  best_pair_span_bp = portfolio_span_bp
)]
triple <- best[portfolio_size == 3L, .(
  gene_id,
  best_triple_sites = portfolio_site_ids,
  best_triple_score = multiplex_robustness_score,
  best_triple_global_atleast1 = global_at_least_one_intact_fraction,
  best_triple_q10_atleast1 = population_q10_at_least_one_intact,
  best_triple_breadth_ge95 = population_breadth_ge_0_95,
  best_triple_span_bp = portfolio_span_bp
)]

# Explicit typed empty tables guarantee downstream columns exist even when no pair/triple was possible.
if (!nrow(pair)) {
  pair <- data.table::data.table(
    gene_id = character(), best_pair_sites = character(), best_pair_score = numeric(),
    best_pair_global_atleast1 = numeric(), best_pair_q10_atleast1 = numeric(),
    best_pair_breadth_ge95 = numeric(), best_pair_span_bp = integer()
  )
}
if (!nrow(triple)) {
  triple <- data.table::data.table(
    gene_id = character(), best_triple_sites = character(), best_triple_score = numeric(),
    best_triple_global_atleast1 = numeric(), best_triple_q10_atleast1 = numeric(),
    best_triple_breadth_ge95 = numeric(), best_triple_span_bp = integer()
  )
}

# Best table must contain at most one row per gene per portfolio size.
if (nrow(pair) && anyDuplicated(pair$gene_id)) fail_step(STEP, "More than one 'best pair' row exists for at least one gene.")
if (nrow(triple) && anyDuplicated(triple$gene_id)) fail_step(STEP, "More than one 'best triple' row exists for at least one gene.")

# Merges are now collision-safe because only new columns are carried from site/portfolio tables.
g <- merge(panel, best_single, by = "gene_id", all.x = TRUE, sort = FALSE)
g <- merge(g, pair, by = "gene_id", all.x = TRUE, sort = FALSE)
g <- merge(g, triple, by = "gene_id", all.x = TRUE, sort = FALSE)
if (nrow(g) != nrow(panel) || data.table::uniqueN(g$gene_id) != nrow(panel)) fail_step(STEP, "Gene-level merges changed the number/uniqueness of panel genes.")

w_single <- unname(gene_weights[["best_single_site"]])
w_pair <- unname(gene_weights[["best_pair"]])
w_triple <- unname(gene_weights[["best_triple"]])
g[, gene_population_robustness_score := mapply(
  function(s, p, t) renorm_weighted(c(s, p, t), c(w_single, w_pair, w_triple)),
  best_single_site_score, best_pair_score, best_triple_score
)]

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

g[, population_deployment_class := mapply(
  classify_gene,
  gene_population_robustness_score,
  best_single_site_class,
  best_pair_global_atleast1,
  best_pair_q10_atleast1,
  USE.NAMES = FALSE
)]

priority <- c(HIGHLY_ROBUST = 1L, MULTIPLEX_RESCUABLE = 2L, INTERMEDIATE = 3L, CONCERN = 4L, UNRESOLVED = 5L)
g[, deployment_class_priority := unname(priority[population_deployment_class])]
if (anyNA(g$deployment_class_priority)) fail_step(STEP, "An unknown population_deployment_class was generated.")

# Pure population deployment rank: score can reorder genes within deployment class.
data.table::setorder(g, deployment_class_priority, -gene_population_robustness_score, final_prepopulation_rank, gene_id)
g[, population_deployment_rank := seq_len(.N)]

# Conservative population-filtered rank: frozen biology-led order is preserved within deployment class.
data.table::setorder(g, deployment_class_priority, final_prepopulation_rank, -gene_population_robustness_score, gene_id)
g[, population_filtered_prepopulation_rank := seq_len(.N)]
# Positive values indicate an upward movement (e.g., frozen rank 30 -> filtered rank 10 gives +20).
g[, population_filtered_rank_shift := final_prepopulation_rank - population_filtered_prepopulation_rank]
g[, population_filtered_rank_direction := data.table::fcase(
  population_filtered_rank_shift > 0, "UP",
  population_filtered_rank_shift < 0, "DOWN",
  default = "UNCHANGED"
)]

atomic_fwrite(g, "data_processed/05_gene_population_robustness.csv")
atomic_fwrite(g[is_external_benchmark == TRUE], "data_processed/05_benchmark_population_robustness.csv")
novel <- g[is_external_benchmark == FALSE][order(population_filtered_prepopulation_rank)]
atomic_fwrite(novel[seq_len(min(50L, nrow(novel)))], "data_processed/05_top50_population_validated_novel.csv")

qc <- g[, .N, by = population_deployment_class]
qc[, deployment_class_priority := unname(priority[population_deployment_class])]
data.table::setorder(qc, deployment_class_priority)
atomic_fwrite(qc, "data_processed/05_gene_robustness_qc.csv")
write_checksum(
  c(
    "data_processed/05_gene_population_robustness.csv",
    "data_processed/05_benchmark_population_robustness.csv",
    "data_processed/05_top50_population_validated_novel.csv"
  ),
  "logs/05_checksums.tsv"
)
write_session_info("logs/05_sessionInfo.txt")
print(qc)
log_step(STEP, "Gene-level population robustness completed successfully")
