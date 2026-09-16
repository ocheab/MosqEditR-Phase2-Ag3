# STEP 05b: freeze the Ag3 external-validation target set BEFORE external data are read.
#
# This script enforces validation isolation. Selection is based only on frozen
# inputs and Phase-2 discovery outputs. Ag3 data must never alter this target set.

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 05b from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))
suppressPackageStartupMessages(library(data.table))

STEP <- "05b"
log_step(STEP, "Freezing external-validation target set from Phase-2 discovery outputs")

required <- c(
  "data_processed/00_frozen_target_sites.csv",
  "data_processed/00_frozen_gene_panel.csv",
  "data_processed/03_population_site_metrics.csv",
  "data_processed/04_best_multiplex_portfolios.csv",
  "data_processed/05_gene_population_robustness.csv"
)
assert_files(required, nonempty = TRUE)

cfg <- read_config()
ev <- cfg_get(cfg, "external_validation", default = list())
include_single <- isTRUE(as.logical(ev$include_best_single_per_gene %||% TRUE))
include_pair <- isTRUE(as.logical(ev$include_best_pair_members %||% TRUE))
include_triple <- isTRUE(as.logical(ev$include_best_triple_members %||% TRUE))
include_bench <- isTRUE(as.logical(ev$include_external_benchmark_genes %||% TRUE))

sites <- fread(required[[1]], na.strings = c("", "NA", "NaN"))
panel <- fread(required[[2]], na.strings = c("", "NA", "NaN"))
site_metrics <- fread(required[[3]], na.strings = c("", "NA", "NaN"))
best <- fread(required[[4]], na.strings = c("", "NA", "NaN"))
gene <- fread(required[[5]], na.strings = c("", "NA", "NaN"))

assert_columns(sites, c(
  "population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end",
  "guide_genomic_strand", "discovery_phased_panel_assessable"
), "frozen target sites")
assert_columns(panel, c("gene_id", "is_external_benchmark"), "frozen gene panel")
assert_columns(site_metrics, c("population_site_id", "gene_id", "population_target_site_robustness_score"), "Phase-2 site metrics")
assert_columns(gene, c("gene_id", "best_single_site_id", "best_pair_sites", "best_triple_sites"), "Phase-2 gene metrics")
assert_unique(sites, "population_site_id", "frozen target sites")
assert_unique(gene, "gene_id", "Phase-2 gene metrics")
panel[, is_external_benchmark := normalize_flag(is_external_benchmark)]
panel[is.na(is_external_benchmark), is_external_benchmark := FALSE]

selection <- data.table(population_site_id = character(), selection_reason = character())
add_selection <- function(ids, reason) {
  ids <- unique(trimws(as.character(ids)))
  ids <- ids[nonempty_string(ids)]
  if (!length(ids)) return(invisible(NULL))
  selection <<- rbind(selection, data.table(population_site_id = ids, selection_reason = reason), fill = TRUE)
  invisible(NULL)
}

if (include_single) add_selection(gene$best_single_site_id, "PHASE2_BEST_SINGLE")

split_portfolios <- function(x) {
  x <- as.character(x)
  unlist(lapply(x[nonempty_string(x)], function(z) strsplit(z, ";", fixed = TRUE)[[1L]]), use.names = FALSE)
}
if (include_pair) add_selection(split_portfolios(gene$best_pair_sites), "PHASE2_BEST_PAIR_MEMBER")
if (include_triple) add_selection(split_portfolios(gene$best_triple_sites), "PHASE2_BEST_TRIPLE_MEMBER")
if (include_bench) {
  bg <- panel[is_external_benchmark == TRUE, gene_id]
  add_selection(sites[gene_id %in% bg, population_site_id], "EXTERNAL_BENCHMARK_GENE_SITE")
}

if (!nrow(selection)) fail_step(STEP, "External-validation target selection produced zero sites.")
selection <- selection[, .(selection_reason = paste(sort(unique(selection_reason)), collapse = ";")), by = population_site_id]
lock <- merge(sites, selection, by = "population_site_id", all = FALSE, sort = FALSE)
if (!nrow(lock)) fail_step(STEP, "No selected external-validation IDs matched the frozen target panel.")
if (nrow(lock) != uniqueN(lock$population_site_id)) fail_step(STEP, "Duplicate external-validation target IDs after merge.")
if (any(lock$discovery_phased_panel_assessable != TRUE)) {
  bad <- lock[discovery_phased_panel_assessable != TRUE, population_site_id]
  fail_step(STEP, paste0("External-validation set contains unresolved/noncanonical target(s): ", paste(head(bad, 20L), collapse = ", ")))
}

# Add frozen discovery annotations for transparent, audit-ready validation.
ann <- site_metrics[, .(
  population_site_id,
  phase2_site_score = population_target_site_robustness_score,
  phase2_site_class = population_robustness_class,
  phase2_global_exact = target_23bp_exact_match_fraction,
  phase2_population_q10_exact = population_q10_target_23bp_exact,
  phase2_pam_intact = pam_intact_fraction,
  phase2_max_alt_af = max_target_variant_alt_af
)]
lock <- merge(lock, ann, by = "population_site_id", all.x = TRUE, sort = FALSE)
lock[, `:=`(
  validation_dataset = "Ag3.0 public Sanger per-sample genotypes",
  discovery_dataset = "Ag1000G Phase 2 AR1 phased haplotypes",
  selection_locked_before_ag3 = TRUE,
  external_validation_can_retune_discovery = FALSE
)]
setorder(lock, final_prepopulation_rank, gene_id, population_site_rank_within_gene, population_site_id)

out <- "data_processed/05b_external_validation_targets.csv"
atomic_fwrite(lock, out)
# The target lock is immutable evidence that external data were not used in selection.
write_sha256(out, "metadata/05b_external_validation_target_lock.sha256")

qc <- data.table(
  metric = c(
    "External-validation frozen sites",
    "Genes represented",
    "Benchmark genes represented",
    "Sites selected as best single",
    "Sites selected as pair member",
    "Sites selected as triple member",
    "Sites selected because benchmark gene",
    "Ag3 rows read during target freezing"
  ),
  value = c(
    nrow(lock),
    uniqueN(lock$gene_id),
    uniqueN(lock[gene_id %in% panel[is_external_benchmark == TRUE, gene_id], gene_id]),
    sum(grepl("PHASE2_BEST_SINGLE", lock$selection_reason, fixed = TRUE)),
    sum(grepl("PHASE2_BEST_PAIR_MEMBER", lock$selection_reason, fixed = TRUE)),
    sum(grepl("PHASE2_BEST_TRIPLE_MEMBER", lock$selection_reason, fixed = TRUE)),
    sum(grepl("EXTERNAL_BENCHMARK_GENE_SITE", lock$selection_reason, fixed = TRUE)),
    0L
  )
)
atomic_fwrite(qc, "data_processed/05b_external_validation_target_qc.csv")
write_session_info("logs/05b_sessionInfo.txt")
print(qc)
log_step(STEP, "External-validation target lock completed; Ag3 has not been accessed")
