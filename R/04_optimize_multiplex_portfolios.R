# STEP 04: optimize within-gene two- and three-site target portfolios using phased haplotypes

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 04 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "04"
log_step(STEP, "Optimizing phased multiplex target portfolios")

required <- c(
  "data_processed/00_frozen_target_sites.csv",
  "data_processed/02_analysis_sample_metadata.csv",
  "data_processed/02_analysis_population_manifest.csv",
  "data_processed/02_contig_population_manifest.csv",
  "data_processed/02_target_haplotype_states.rds",
  "data_processed/03_population_site_metrics.csv"
)
assert_files(required)

cfg <- read_config()
portfolio_sizes <- sort(unique(as.integer(unlist(cfg_get(cfg, "multiplex_portfolio_sizes", default = c(2L, 3L))))))
portfolio_sizes <- portfolio_sizes[is.finite(portfolio_sizes) & portfolio_sizes >= 2L]
if (!length(portfolio_sizes)) fail_step(STEP, "No valid multiplex_portfolio_sizes >=2 are configured.")
weights <- unlist(cfg_get(cfg, "multiplex_robustness_weights", required = TRUE), use.names = TRUE)
validate_weights(weights, "multiplex_robustness_weights")
needed_w <- c("global_at_least_one_intact", "population_q10_at_least_one_intact", "population_breadth_ge_0_95")
if (!all(needed_w %in% names(weights))) fail_step(STEP, paste0("multiplex_robustness_weights must define: ", paste(needed_w, collapse = ", ")))
weights <- weights[needed_w]

sites <- data.table::fread(required[[1]], na.strings = c("", "NA", "NaN"))
meta <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
pops <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))
contig_pops <- data.table::fread(required[[4]], na.strings = c("", "NA", "NaN"))
states <- readRDS(required[[5]])
site_metrics <- data.table::fread(required[[6]], na.strings = c("", "NA", "NaN"))

assert_columns(sites, c("population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end", "discovery_phased_panel_assessable"), "frozen target sites")
assert_columns(meta, c("sample_id", "primary_taxon", "population_id"), "analysis sample metadata")
assert_columns(pops, c("population_id", "n_samples", "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30"), "analysis population manifest")
assert_columns(contig_pops, c("contig", "population_id", "n_samples", "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30"), "contig-specific population manifest")
assert_columns(site_metrics, c("population_site_id", "gene_id", "population_robustness_class"), "site metrics")
assert_unique(site_metrics, "population_site_id", "site metrics")
if (anyDuplicated(contig_pops[, .(contig, population_id)])) {
  fail_step(STEP, "Contig-specific population manifest is not unique by (contig, population_id).")
}


if (is.null(states$exact_23bp) || !is.matrix(states$exact_23bp)) fail_step(STEP, "states$exact_23bp is missing or malformed.")
state_ids <- rownames(states$exact_23bp)
hap_names <- colnames(states$exact_23bp)
if (is.null(state_ids) || is.null(hap_names)) fail_step(STEP, "Haplotype-state matrix requires row and column names.")

meta[, primary_taxon := normalize_flag(primary_taxon)]
analysis_meta <- meta[primary_taxon == TRUE]
hap_samples <- haplotype_to_sample(hap_names)
if (any(!hap_samples %in% analysis_meta$sample_id)) fail_step(STEP, "Haplotype-state columns contain samples absent from analysis metadata.")
hap_pops <- analysis_meta$population_id[match(hap_samples, analysis_meta$sample_id)]
global_hap <- seq_along(hap_names)
all_pops_n10 <- unique(contig_pops[sensitivity_population_n10 == TRUE, population_id])
if (!length(all_pops_n10)) fail_step(STEP, "No contig-specific n>=10 populations are available for portfolio sensitivity summaries.")

eligible_site_ids <- site_metrics[
  population_robustness_class %in% c("ROBUST", "INTERMEDIATE", "CONCERN"),
  population_site_id
]
eligible_site_ids <- intersect(eligible_site_ids, state_ids)

portfolio_state <- function(mat) apply(mat, 2L, at_least_one_true)
all_intact_state <- function(mat) apply(mat, 2L, all_true_state)

rows <- list()
pop_rows <- list()
ri <- 0L
pi <- 0L

for (gene in unique(sites$gene_id)) {
  gene_sites <- sites[
    gene_id == gene & discovery_phased_panel_assessable == TRUE & population_site_id %in% eligible_site_ids
  ]
  sids <- gene_sites$population_site_id[gene_sites$population_site_id %in% state_ids]
  sids <- unique(sids)
  if (length(sids) < 2L) next

  # All target sites within one gene must be on one contig; Step 00 checks this, but re-check locally.
  if (data.table::uniqueN(gene_sites$genomic_seqid) != 1L) {
    fail_step(STEP, paste0("Gene ", gene, " has eligible target sites on multiple contigs; genomic span is undefined."))
  }
  gene_contig <- unique(as.character(gene_sites$genomic_seqid))[[1]]
  gene_pops_n10 <- unique(contig_pops[contig == gene_contig & sensitivity_population_n10 == TRUE, population_id])
  if (!length(gene_pops_n10)) next
  gene_pop_manifest <- contig_pops[contig == gene_contig, .(
    population_id, n_samples, primary_population_n20, sensitivity_population_n10, sensitivity_population_n30
  )]

  for (k in portfolio_sizes[portfolio_sizes <= length(sids)]) {
    combos <- utils::combn(sids, k, simplify = FALSE)
    for (combo in combos) {
      m <- states$exact_23bp[combo, global_hap, drop = FALSE]
      atleast <- portfolio_state(m)
      allint <- all_intact_state(m)
      combo_id <- paste(combo, collapse = ";")

      per_pop_list <- lapply(gene_pops_n10, function(pop) {
        idx <- which(hap_pops == pop)
        x <- atleast[idx]
        y <- allint[idx]
        data.table::data.table(
          gene_id = gene,
          portfolio_size = as.integer(k),
          portfolio_site_ids = combo_id,
          population_id = pop,
          at_least_one_intact_fraction = if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
          all_targets_intact_fraction = if (!length(y) || all(is.na(y))) NA_real_ else mean(y, na.rm = TRUE),
          n_haplotypes_total = length(idx),
          n_haplotypes_callable_atleast1 = sum(!is.na(x)),
          n_haplotypes_callable_all = sum(!is.na(y))
        )
      })
      per_pop <- data.table::rbindlist(per_pop_list, fill = TRUE, use.names = TRUE)
      per_pop <- merge(
        per_pop,
        gene_pop_manifest,
        by = "population_id", all.x = TRUE, sort = FALSE
      )
      if (anyNA(per_pop$n_samples)) fail_step(STEP, paste0("Portfolio population rows failed to match population manifest for gene ", gene, "."))

      pp20 <- per_pop[primary_population_n20 == TRUE]
      popv <- pp20$at_least_one_intact_fraction
      global_atleast <- if (all(is.na(atleast))) NA_real_ else mean(atleast, na.rm = TRUE)
      global_all <- if (all(is.na(allint))) NA_real_ else mean(allint, na.rm = TRUE)
      q10 <- safe_quantile(popv, 0.10)
      breadth <- if (any(is.finite(popv))) mean(popv[is.finite(popv)] >= 0.95) else NA_real_
      score_components <- c(
        global_at_least_one_intact = global_atleast,
        population_q10_at_least_one_intact = q10,
        population_breadth_ge_0_95 = breadth
      )
      score <- renorm_weighted(score_components, weights[names(score_components)])

      coords <- sites[match(combo, population_site_id)]
      if (nrow(coords) != length(combo) || anyNA(coords$genomic_start) || anyNA(coords$genomic_end)) {
        fail_step(STEP, paste0("Could not resolve genomic coordinates for portfolio ", combo_id, "."))
      }
      span <- max(coords$genomic_end) - min(coords$genomic_start) + 1L

      ri <- ri + 1L
      rows[[ri]] <- data.table::data.table(
        gene_id = gene,
        portfolio_size = as.integer(k),
        portfolio_site_ids = combo_id,
        n_sites_callable = nrow(m),
        portfolio_span_bp = as.integer(span),
        global_at_least_one_intact_fraction = global_atleast,
        population_q10_at_least_one_intact = q10,
        population_breadth_ge_0_95 = breadth,
        global_all_targets_intact_fraction = global_all,
        populations_assessed = sum(is.finite(popv)),
        multiplex_robustness_score = score
      )
      pi <- pi + 1L
      pop_rows[[pi]] <- per_pop
    }
  }
}

portfolio_schema <- list(
  gene_id = character(), portfolio_size = integer(), portfolio_site_ids = character(), n_sites_callable = integer(),
  portfolio_span_bp = integer(), global_at_least_one_intact_fraction = numeric(),
  population_q10_at_least_one_intact = numeric(), population_breadth_ge_0_95 = numeric(),
  global_all_targets_intact_fraction = numeric(), populations_assessed = integer(),
  multiplex_robustness_score = numeric(), portfolio_rank_within_gene_size = integer()
)
pop_schema <- list(
  population_id = character(), gene_id = character(), portfolio_size = integer(), portfolio_site_ids = character(),
  at_least_one_intact_fraction = numeric(), all_targets_intact_fraction = numeric(),
  n_haplotypes_total = integer(), n_haplotypes_callable_atleast1 = integer(), n_haplotypes_callable_all = integer(),
  n_samples = integer(), primary_population_n20 = logical(), sensitivity_population_n10 = logical(), sensitivity_population_n30 = logical()
)

portfolio <- if (length(rows)) data.table::rbindlist(rows, fill = TRUE, use.names = TRUE) else empty_table(portfolio_schema[names(portfolio_schema) != "portfolio_rank_within_gene_size"])
portfolio_pop <- if (length(pop_rows)) data.table::rbindlist(pop_rows, fill = TRUE, use.names = TRUE) else empty_table(pop_schema)

if (nrow(portfolio)) {
  portfolio[, `:=`(
    .sort_score = data.table::fifelse(is.finite(multiplex_robustness_score), multiplex_robustness_score, -Inf),
    .sort_q10 = data.table::fifelse(is.finite(population_q10_at_least_one_intact), population_q10_at_least_one_intact, -Inf),
    .sort_global = data.table::fifelse(is.finite(global_at_least_one_intact_fraction), global_at_least_one_intact_fraction, -Inf),
    .sort_breadth = data.table::fifelse(is.finite(population_breadth_ge_0_95), population_breadth_ge_0_95, -Inf)
  )]
  data.table::setorder(
    portfolio,
    gene_id, portfolio_size,
    -.sort_score, -.sort_q10, -.sort_global, -.sort_breadth,
    portfolio_span_bp, portfolio_site_ids
  )
  portfolio[, portfolio_rank_within_gene_size := seq_len(.N), by = .(gene_id, portfolio_size)]
  portfolio[, c(".sort_score", ".sort_q10", ".sort_global", ".sort_breadth") := NULL]
} else {
  portfolio[, portfolio_rank_within_gene_size := integer()]
}

atomic_fwrite(portfolio, "data_processed/04_multiplex_portfolio_metrics.csv")
atomic_fwrite(portfolio_pop, "data_processed/04_multiplex_portfolio_by_population.csv")
best <- if (nrow(portfolio)) portfolio[portfolio_rank_within_gene_size == 1L] else portfolio
atomic_fwrite(best, "data_processed/04_best_multiplex_portfolios.csv")

qc <- data.table::data.table(
  metric = c(
    "Eligible single sites", "Evaluated portfolios", "Genes with pair portfolio", "Genes with triple portfolio"
  ),
  value = c(
    length(eligible_site_ids), nrow(portfolio),
    data.table::uniqueN(portfolio[portfolio_size == 2L, gene_id]),
    data.table::uniqueN(portfolio[portfolio_size == 3L, gene_id])
  )
)
atomic_fwrite(qc, "data_processed/04_multiplex_qc.csv")
write_checksum(
  c("data_processed/04_multiplex_portfolio_metrics.csv", "data_processed/04_multiplex_portfolio_by_population.csv", "data_processed/04_best_multiplex_portfolios.csv"),
  "logs/04_checksums.tsv"
)
write_session_info("logs/04_sessionInfo.txt")
print(qc)
log_step(STEP, "Multiplex optimization completed successfully")
