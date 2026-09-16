# STEP 03: compute site-level and site-by-population robustness metrics

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 03 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))

STEP <- "03"
log_step(STEP, "Computing site-level population genomic robustness metrics")

required <- c(
  "data_processed/00_frozen_target_sites.csv",
  "data_processed/02_analysis_sample_metadata.csv",
  "data_processed/02_analysis_population_manifest.csv",
  "data_processed/02_contig_population_manifest.csv",
  "data_processed/02_target_haplotype_states.rds",
  "data_processed/02_target_variant_genotypes.rds",
  "data_processed/02_target_accessibility.csv"
)
assert_files(required)

cfg <- read_config()
weights <- unlist(cfg_get(cfg, "site_robustness_weights", required = TRUE), use.names = TRUE)
validate_weights(weights, "site_robustness_weights")
required_weight_names <- c(
  "global_23bp_exact", "population_q10_23bp_exact", "population_breadth_ge_0_95",
  "global_pam_intact", "one_minus_max_alt_af_any_population", "accessibility_fraction_23bp"
)
if (!all(required_weight_names %in% names(weights))) {
  fail_step(STEP, paste0("site_robustness_weights must define: ", paste(required_weight_names, collapse = ", ")))
}
weights <- weights[required_weight_names]

robust_th <- cfg_get(cfg, "site_class_thresholds", "robust", required = TRUE)
intermediate_th <- cfg_get(cfg, "site_class_thresholds", "intermediate", required = TRUE)
accessibility_min <- as.numeric(cfg_get(cfg, "accessibility_min_for_classification", default = 0.90))

sites <- data.table::fread(required[[1]], na.strings = c("", "NA", "NaN"))
meta <- data.table::fread(required[[2]], na.strings = c("", "NA", "NaN"))
pops <- data.table::fread(required[[3]], na.strings = c("", "NA", "NaN"))
contig_pops <- data.table::fread(required[[4]], na.strings = c("", "NA", "NaN"))
states <- readRDS(required[[5]])
variants <- readRDS(required[[6]])
access <- data.table::fread(required[[7]], na.strings = c("", "NA", "NaN"))

assert_columns(sites, c("population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end", "discovery_phased_panel_assessable"), "frozen target sites")
assert_columns(meta, c("sample_id", "primary_taxon", "population_id"), "analysis sample metadata")
assert_columns(pops, c("population_id", "n_samples", "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30"), "analysis population manifest")
assert_columns(contig_pops, c("contig", "population_id", "n_samples", "primary_population_n20", "sensitivity_population_n10", "sensitivity_population_n30"), "contig-specific population manifest")
assert_columns(access, c("population_site_id", "accessibility_fraction_23bp", "accessibility_complete"), "target accessibility table")
assert_unique(access, "population_site_id", "target accessibility table")

# Defensive normalization of the contig-specific population manifest.
# A contig/population pair must have exactly one set of denominator/eligibility values.
# Earlier X-aware handoff versions could contain exact duplicate rows; those duplicates
# must not be allowed to multiply site-by-population rows in later joins.
manifest_flag_cols <- c(
  "primary_population_n20",
  "sensitivity_population_n10",
  "sensitivity_population_n30"
)
for (nm in manifest_flag_cols) {
  contig_pops[, (nm) := normalize_flag(get(nm))]
}

dup_manifest <- contig_pops[
  ,
  .(
    rows = .N,
    n_n_samples = data.table::uniqueN(n_samples, na.rm = FALSE),
    n_primary_population_n20 = data.table::uniqueN(primary_population_n20, na.rm = FALSE),
    n_sensitivity_population_n10 = data.table::uniqueN(sensitivity_population_n10, na.rm = FALSE),
    n_sensitivity_population_n30 = data.table::uniqueN(sensitivity_population_n30, na.rm = FALSE)
  ),
  by = .(contig, population_id)
][rows > 1L]

if (nrow(dup_manifest)) {
  conflict_manifest <- dup_manifest[
    n_n_samples > 1L |
      n_primary_population_n20 > 1L |
      n_sensitivity_population_n10 > 1L |
      n_sensitivity_population_n30 > 1L
  ]

  if (nrow(conflict_manifest)) {
    preview <- paste(
      utils::head(
        paste0(conflict_manifest$contig, "/", conflict_manifest$population_id),
        20L
      ),
      collapse = ", "
    )
    fail_step(
      STEP,
      paste0(
        "Conflicting duplicate rows were found in the contig-specific population manifest for: ",
        preview,
        ". Step 03 will not choose between inconsistent population denominators."
      )
    )
  }

  n_before <- nrow(contig_pops)
  data.table::setorder(contig_pops, contig, population_id)
  contig_pops <- unique(contig_pops, by = c("contig", "population_id"))
  n_removed <- n_before - nrow(contig_pops)

  log_step(
    STEP,
    paste0(
      "Collapsed ", n_removed,
      " exactly consistent duplicate contig/population manifest row(s); ",
      nrow(contig_pops), " unique contig/population rows remain."
    )
  )
}

if (anyDuplicated(contig_pops[, .(contig, population_id)])) {
  fail_step(STEP, "Contig-specific population manifest is not unique by contig + population_id after normalization.")
}

for (nm in c("exact_23bp", "pam_intact", "protospacer_exact")) {
  if (is.null(states[[nm]]) || !is.matrix(states[[nm]])) fail_step(STEP, paste0("states$", nm, " is missing or not a matrix."))
}
if (!identical(dim(states$exact_23bp), dim(states$pam_intact)) || !identical(dim(states$exact_23bp), dim(states$protospacer_exact))) {
  fail_step(STEP, "Haplotype-state matrices do not have identical dimensions.")
}
if (!identical(rownames(states$exact_23bp), rownames(states$pam_intact)) ||
    !identical(rownames(states$exact_23bp), rownames(states$protospacer_exact)) ||
    !identical(colnames(states$exact_23bp), colnames(states$pam_intact)) ||
    !identical(colnames(states$exact_23bp), colnames(states$protospacer_exact))) {
  fail_step(STEP, "Haplotype-state matrices do not have identical row/column identities.")
}
assert_unique(meta, "sample_id", "analysis sample metadata")
meta[, primary_taxon := normalize_flag(primary_taxon)]
analysis_meta <- meta[primary_taxon == TRUE]

hap_names <- colnames(states$exact_23bp)
if (is.null(hap_names) || !length(hap_names)) fail_step(STEP, "Haplotype-state matrix has no named columns.")
hap_to_sample <- haplotype_to_sample(hap_names)
if (any(!hap_to_sample %in% analysis_meta$sample_id)) {
  bad <- unique(hap_to_sample[!hap_to_sample %in% analysis_meta$sample_id])
  fail_step(STEP, paste0("Haplotype columns contain sample IDs absent from analysis metadata: ", paste(utils::head(bad, 20L), collapse = ", ")))
}
hap_to_pop <- analysis_meta$population_id[match(hap_to_sample, analysis_meta$sample_id)]
global_hap_idx <- seq_along(hap_names)

sensitivity_pops_n10 <- unique(contig_pops[sensitivity_population_n10 == TRUE, population_id])
primary_pops_n20 <- unique(contig_pops[primary_population_n20 == TRUE, population_id])
if (!length(sensitivity_pops_n10)) fail_step(STEP, "No contig-specific n>=10 populations are available after phased-panel intersection.")
if (!length(primary_pops_n20)) fail_step(STEP, "No contig-specific n>=20 primary populations are available after phased-panel intersection.")

# Build a site-to-variant-position map so an invariant target (no VCF records) can be
# distinguished from a variant-bearing target with no called alleles in a specific population.
site_variant_rows <- list()
svi <- 0L
for (contig in names(variants)) {
  obj <- variants[[contig]]
  if (is.null(obj$variants) || !nrow(obj$variants)) next
  ss <- sites[genomic_seqid == contig & discovery_phased_panel_assessable == TRUE]
  if (!nrow(ss)) next
  for (v in seq_len(nrow(obj$variants))) {
    pos <- as.integer(obj$variants$position[[v]])
    hit <- ss[genomic_start <= pos & genomic_end >= pos, population_site_id]
    if (!length(hit)) next
    for (sid in hit) {
      svi <- svi + 1L
      site_variant_rows[[svi]] <- data.table::data.table(population_site_id = sid, contig = contig, position = pos)
    }
  }
}
site_variant_map <- if (length(site_variant_rows)) {
  unique(data.table::rbindlist(site_variant_rows, use.names = TRUE, fill = TRUE), by = c("population_site_id", "contig", "position"))
} else {
  data.table::data.table(population_site_id = character(), contig = character(), position = integer())
}
atomic_fwrite(site_variant_map, "data_processed/03_target_variant_position_map.csv")

# Variant allele frequencies by site and population for every n>=10 group.
variant_population_rows <- list()
vp_i <- 0L
for (contig in names(variants)) {
  obj <- variants[[contig]]
  if (is.null(obj$variants) || !nrow(obj$variants)) next
  if (is.null(obj$gt) || !is.matrix(obj$gt)) fail_step(STEP, paste0("Variant genotype object for ", contig, " is malformed."))
  if (nrow(obj$gt) != nrow(obj$variants)) fail_step(STEP, paste0("Variant table and genotype matrix row counts differ for ", contig, "."))

  gt <- obj$gt
  if (is.null(colnames(gt))) fail_step(STEP, paste0("Variant genotype matrix for ", contig, " has no sample names."))
  sample_pop <- analysis_meta$population_id[match(colnames(gt), analysis_meta$sample_id)]
  if (any(is.na(match(colnames(gt), analysis_meta$sample_id)))) {
    fail_step(STEP, paste0("Variant genotype matrix for ", contig, " contains sample IDs absent from analysis metadata."))
  }
  ss <- sites[genomic_seqid == contig & discovery_phased_panel_assessable == TRUE]
  ctg <- contig
  eligible_pops_contig_n10 <- unique(contig_pops[contig == ctg & sensitivity_population_n10 == TRUE, population_id])

  for (v in seq_len(nrow(obj$variants))) {
    pos <- as.integer(obj$variants$position[[v]])
    hit_sites <- ss[genomic_start <= pos & genomic_end >= pos, population_site_id]
    if (!length(hit_sites)) next

    parsed <- parse_gt_haplotypes(gt[v, ], require_phased = TRUE)
    alleles <- c(parsed$a1, parsed$a2)
    allele_pop <- c(sample_pop, sample_pop)
    alt_indicator <- rep(NA_integer_, length(alleles))
    called <- !is.na(alleles)
    alt_indicator[called] <- as.integer(alleles[called] != "0")

    for (pop in eligible_pops_contig_n10) {
      x <- alt_indicator[allele_pop == pop]
      x <- x[!is.na(x)]
      if (!length(x)) next
      af <- mean(x)
      for (sid in hit_sites) {
        vp_i <- vp_i + 1L
        variant_population_rows[[vp_i]] <- data.table::data.table(
          population_site_id = sid,
          population_id = pop,
          contig = contig,
          position = pos,
          alt_af = af,
          heterozygosity_proxy = 2 * af * (1 - af),
          n_called_haplotypes = length(x)
        )
      }
    }
  }
}
variant_pop <- if (length(variant_population_rows)) {
  data.table::rbindlist(variant_population_rows, fill = TRUE, use.names = TRUE)
} else {
  data.table::data.table(
    population_site_id = character(), population_id = character(), contig = character(),
    position = integer(), alt_af = numeric(), heterozygosity_proxy = numeric(), n_called_haplotypes = integer()
  )
}
atomic_fwrite(variant_pop, "data_processed/03_target_variant_population_af.csv")

# Site x population metrics for all n>=10 populations.
site_pop_rows <- list()
sp_i <- 0L
state_site_ids <- rownames(states$exact_23bp)
for (sid in state_site_ids) {
  ctg <- sites[population_site_id == sid, genomic_seqid][[1]]
  eligible_pops_site_n10 <- unique(contig_pops[contig == ctg & sensitivity_population_n10 == TRUE, population_id])
  ex_all <- states$exact_23bp[sid, , drop = TRUE]
  pa_all <- states$pam_intact[sid, , drop = TRUE]
  pr_all <- states$protospacer_exact[sid, , drop = TRUE]
  n_variant_positions_site <- data.table::uniqueN(site_variant_map[population_site_id == sid, position])

  for (pop in eligible_pops_site_n10) {
    idx <- which(hap_to_pop == pop)
    if (!length(idx)) next
    ex <- ex_all[idx]
    pa <- pa_all[idx]
    pr <- pr_all[idx]
    vp <- variant_pop[population_site_id == sid & population_id == pop]

    max_af <- if (n_variant_positions_site == 0L) {
      0
    } else if (nrow(vp) && any(is.finite(vp$alt_af))) {
      safe_max(vp$alt_af)
    } else {
      NA_real_
    }
    hetero <- if (n_variant_positions_site == 0L) {
      0
    } else if (nrow(vp) && any(is.finite(vp$heterozygosity_proxy))) {
      sum(vp$heterozygosity_proxy[is.finite(vp$heterozygosity_proxy)]) / 23
    } else {
      NA_real_
    }

    sp_i <- sp_i + 1L
    site_pop_rows[[sp_i]] <- data.table::data.table(
      population_site_id = sid,
      contig = ctg,
      population_id = pop,
      n_haplotypes_total = length(idx),
      n_haplotypes_called_23bp = sum(!is.na(ex)),
      callable_fraction_23bp = mean(!is.na(ex)),
      target_23bp_exact_match_fraction = if (all(is.na(ex))) NA_real_ else mean(ex, na.rm = TRUE),
      pam_intact_fraction = if (all(is.na(pa))) NA_real_ else mean(pa, na.rm = TRUE),
      protospacer_exact_match_fraction = if (all(is.na(pr))) NA_real_ else mean(pr, na.rm = TRUE),
      max_target_variant_alt_af = max_af,
      mean_target_heterozygosity_proxy = hetero,
      n_variable_positions = data.table::uniqueN(vp$position),
      n_variant_positions_in_panel = n_variant_positions_site
    )
  }
}
if (!length(site_pop_rows)) fail_step(STEP, "No site-by-population metrics could be computed.")
site_pop <- data.table::rbindlist(site_pop_rows, fill = TRUE, use.names = TRUE)

dup_site_pop <- site_pop[
  ,
  .N,
  by = .(population_site_id, contig, population_id)
][N > 1L]
if (nrow(dup_site_pop)) {
  preview <- paste(
    utils::head(
      paste0(
        dup_site_pop$population_site_id, "/",
        dup_site_pop$contig, "/",
        dup_site_pop$population_id
      ),
      20L
    ),
    collapse = ", "
  )
  fail_step(
    STEP,
    paste0(
      "Duplicate site-by-population rows were generated before manifest merge: ",
      preview,
      ". This indicates duplicated eligibility groups upstream."
    )
  )
}

site_pop <- merge(
  site_pop,
  contig_pops[, .(contig, population_id, n_samples, primary_population_n20, sensitivity_population_n10, sensitivity_population_n30)],
  by = c("contig", "population_id"), all.x = TRUE, sort = FALSE
)
if (anyNA(site_pop$n_samples)) fail_step(STEP, "Some site-by-population rows could not be matched to the population manifest.")
if (anyDuplicated(site_pop[, .(population_site_id, contig, population_id)])) {
  fail_step(STEP, "Manifest merge produced duplicate site-by-population rows.")
}
atomic_fwrite(site_pop, "data_processed/03_site_population_metrics.csv")

# Primary site summaries: global metrics use all phased-panel primary samples;
# geographic summaries use prespecified country x taxon populations with n>=20.
summary_rows <- lapply(state_site_ids, function(sid) {
  ctg <- sites[population_site_id == sid, genomic_seqid][[1]]
  ex <- states$exact_23bp[sid, global_hap_idx, drop = TRUE]
  pa <- states$pam_intact[sid, global_hap_idx, drop = TRUE]
  pr <- states$protospacer_exact[sid, global_hap_idx, drop = TRUE]
  sp <- site_pop[population_site_id == sid & primary_population_n20 == TRUE]
  a <- access[population_site_id == sid]

  pop_exact <- sp$target_23bp_exact_match_fraction
  global_exact <- if (all(is.na(ex))) NA_real_ else mean(ex, na.rm = TRUE)
  global_pam <- if (all(is.na(pa))) NA_real_ else mean(pa, na.rm = TRUE)
  global_proto <- if (all(is.na(pr))) NA_real_ else mean(pr, na.rm = TRUE)
  q10 <- safe_quantile(pop_exact, 0.10)
  breadth <- if (any(is.finite(pop_exact))) mean(pop_exact[is.finite(pop_exact)] >= 0.95) else NA_real_
  max_af <- if (nrow(sp) && any(is.finite(sp$max_target_variant_alt_af))) safe_max(sp$max_target_variant_alt_af) else NA_real_
  acc <- if (nrow(a) == 1L && isTRUE(a$accessibility_complete[[1]])) as.numeric(a$accessibility_fraction_23bp[[1]]) else NA_real_

  score_components <- c(
    global_23bp_exact = global_exact,
    population_q10_23bp_exact = q10,
    population_breadth_ge_0_95 = breadth,
    global_pam_intact = global_pam,
    one_minus_max_alt_af_any_population = if (is.finite(max_af)) 1 - max_af else NA_real_,
    accessibility_fraction_23bp = acc
  )
  score <- renorm_weighted(score_components, weights[names(score_components)])

  data.table::data.table(
    population_site_id = sid,
    population_n_samples = sum(contig_pops[contig == ctg, n_samples], na.rm = TRUE),
    population_n_populations_eligible_n20 = nrow(sp),
    population_callable_fraction = mean(!is.na(ex)),
    pam_intact_fraction = global_pam,
    protospacer_exact_match_fraction = global_proto,
    target_23bp_exact_match_fraction = global_exact,
    population_q10_target_23bp_exact = q10,
    population_breadth_ge_0_95 = breadth,
    max_target_variant_alt_af = max_af,
    one_minus_max_target_variant_alt_af = if (is.finite(max_af)) 1 - max_af else NA_real_,
    mean_target_heterozygosity_proxy = if (nrow(sp) && any(is.finite(sp$mean_target_heterozygosity_proxy))) safe_mean(sp$mean_target_heterozygosity_proxy) else NA_real_,
    populations_with_target_23bp_exact_match_ge_0_95 = sum(pop_exact >= 0.95, na.rm = TRUE),
    populations_assessed = sum(is.finite(pop_exact)),
    accessibility_fraction_23bp = acc,
    population_target_site_robustness_score = score
  )
})
site_summary <- data.table::rbindlist(summary_rows, fill = TRUE, use.names = TRUE)

site_global <- merge(sites, site_summary, by = "population_site_id", all.x = TRUE, sort = FALSE)
site_global <- merge(
  site_global,
  access[, .(population_site_id, accessibility_records_expected, accessibility_records_found, accessibility_complete)],
  by = "population_site_id", all.x = TRUE, sort = FALSE
)

# Explicit classification prevents missing geographic evidence from being mislabeled as CONCERN.
robust_global_min <- as.numeric(robust_th$global_23bp_exact_min)
robust_q10_min <- as.numeric(robust_th$population_q10_23bp_exact_min)
robust_pam_min <- as.numeric(robust_th$global_pam_intact_min)
robust_maxaf_max <- as.numeric(robust_th$max_alt_af_any_population_max)
int_global_min <- as.numeric(intermediate_th$global_23bp_exact_min)
int_q10_min <- as.numeric(intermediate_th$population_q10_23bp_exact_min)
int_pam_min <- as.numeric(intermediate_th$global_pam_intact_min)

site_global[, population_robustness_class := "CONCERN"]
site_global[discovery_phased_panel_assessable != TRUE, population_robustness_class := "UNRESOLVED_CONTIG"]
site_global[
  discovery_phased_panel_assessable == TRUE &
    (is.na(accessibility_complete) | accessibility_complete != TRUE | !is.finite(accessibility_fraction_23bp)),
  population_robustness_class := "UNRESOLVED_ACCESSIBILITY"
]
site_global[
  discovery_phased_panel_assessable == TRUE & !is.na(accessibility_complete) & accessibility_complete &
    is.finite(accessibility_fraction_23bp) & accessibility_fraction_23bp < accessibility_min,
  population_robustness_class := "UNRESOLVED_ACCESSIBILITY"
]
site_global[
  discovery_phased_panel_assessable == TRUE &
    population_robustness_class == "CONCERN" &
    populations_assessed < 1L,
  population_robustness_class := "UNRESOLVED_POPULATION"
]
site_global[
  population_robustness_class == "CONCERN" &
    is.finite(target_23bp_exact_match_fraction) & target_23bp_exact_match_fraction >= robust_global_min &
    is.finite(population_q10_target_23bp_exact) & population_q10_target_23bp_exact >= robust_q10_min &
    is.finite(pam_intact_fraction) & pam_intact_fraction >= robust_pam_min &
    is.finite(max_target_variant_alt_af) & max_target_variant_alt_af <= robust_maxaf_max,
  population_robustness_class := "ROBUST"
]
site_global[
  population_robustness_class == "CONCERN" &
    is.finite(target_23bp_exact_match_fraction) & target_23bp_exact_match_fraction >= int_global_min &
    is.finite(population_q10_target_23bp_exact) & population_q10_target_23bp_exact >= int_q10_min &
    is.finite(pam_intact_fraction) & pam_intact_fraction >= int_pam_min,
  population_robustness_class := "INTERMEDIATE"
]

site_global[, `:=`(
  population_source_release = "Ag1000G Phase 2 AR1 phased discovery",
  population_taxa = paste(as.character(cfg_get(cfg, "primary_taxa", default = c("gambiae", "coluzzii"))), collapse = ";"),
  population_accessibility_mask = "Ag1000G Phase 2 AR1 accessibility mask",
  population_metric_status = data.table::fcase(
    population_robustness_class == "UNRESOLVED_CONTIG", "UNRESOLVED_CONTIG_NOT_IN_PHASED_PANEL",
    population_robustness_class == "UNRESOLVED_ACCESSIBILITY", "UNRESOLVED_ACCESSIBILITY",
    population_robustness_class == "UNRESOLVED_POPULATION", "UNRESOLVED_POPULATION_SUMMARY",
    default = "REAL_PHASE2_DISCOVERY_DATA"
  )
)]

if (nrow(site_global) != nrow(sites) || data.table::uniqueN(site_global$population_site_id) != nrow(sites)) {
  fail_step(STEP, "Final site metric table is not one row per frozen target site.")
}
data.table::setorder(site_global, final_prepopulation_rank, gene_id, population_site_rank_within_gene, population_site_id)
atomic_fwrite(site_global, "data_processed/03_population_site_metrics.csv")

qc <- site_global[, .N, by = population_robustness_class][order(population_robustness_class)]
atomic_fwrite(qc, "data_processed/03_population_site_metrics_qc.csv")
write_checksum(
  c(
    "data_processed/03_target_variant_position_map.csv", "data_processed/03_target_variant_population_af.csv",
    "data_processed/03_site_population_metrics.csv", "data_processed/03_population_site_metrics.csv"
  ),
  "logs/03_checksums.tsv"
)
write_session_info("logs/03_sessionInfo.txt")
print(qc)
log_step(STEP, "Population site metrics completed successfully")
