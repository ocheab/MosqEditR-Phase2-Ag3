# STEP 02: convert locally acquired targeted Ag3 phased data into analysis objects
#
# This R step is deliberately OFFLINE. It reads the small local handoff created by
# scripts/acquire_ag3_targeted.py and never opens a URL, gs:// path, or remote VCF.

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run Step 02 from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml"))
suppressPackageStartupMessages(library(data.table))

STEP <- "02"
log_step(STEP, "Starting local-only conversion of targeted Ag3.0 phased-panel data")
log_step(STEP, "No network access and no chromosome-scale VCF download are performed in this R step")

assert_files(c(
  "data_processed/00_frozen_target_sites.csv",
  "data_processed/01_ag3_sample_metadata.csv",
  "data_processed/01_ag3_population_manifest.csv"
))

cfg <- read_config()
canonical <- as.character(cfg_get(cfg, "canonical_contigs", required = TRUE))
local_dir <- as.character(cfg_get(cfg, "local_acquisition", "output_dir", default = "data_raw/ag3_local"))
release <- as.character(cfg_get(cfg, "local_acquisition", "release", default = cfg_get(cfg, "ag3_release", default = "3.0")))
phasing_panel <- as.character(cfg_get(cfg, "local_acquisition", "phasing_analysis", default = cfg_get(cfg, "phasing_panel", default = "gamb_colu")))
primary_min_n <- as.integer(cfg_get(cfg, "primary_min_population_n", default = 20L))
sensitivity_ns <- unique(c(10L, 30L, as.integer(unlist(cfg_get(cfg, "sensitivity_min_population_n", default = c(10L, 30L))))))

complete_file <- file.path(local_dir, "ACQUISITION_COMPLETE.ok")
access_file <- file.path(local_dir, "accessibility", "target_accessibility.tsv")
acq_qc_file <- file.path(local_dir, "acquisition_qc.tsv")
assert_files(c(complete_file, access_file))

sites <- fread("data_processed/00_frozen_target_sites.csv", na.strings = c("", "NA", "NaN"))
meta <- fread("data_processed/01_ag3_sample_metadata.csv", na.strings = c("", "NA", "NaN"))
access_dt <- fread(access_file, sep = "\t", na.strings = c("", "NA", "NaN"))

assert_columns(sites, c(
  "population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end",
  "guide_genomic_strand", "ag3_phased_panel_assessable"
), "frozen target sites")
assert_columns(meta, c("sample_id", "primary_taxon", "population_id", "country", "taxon"), "Ag3 sample metadata")
assert_columns(access_dt, c(
  "population_site_id", "accessibility_records_expected", "accessibility_records_found",
  "accessibility_complete", "accessibility_fraction_23bp", "pam_accessibility_fraction",
  "protospacer_accessibility_fraction"
), "local target accessibility")
assert_unique(access_dt, "population_site_id", "local target accessibility")

meta[, primary_taxon := normalize_flag(primary_taxon)]
primary_meta <- meta[primary_taxon == TRUE]
if (!nrow(primary_meta)) fail_step(STEP, "No primary-taxon samples are available from Step 01.")
assert_unique(primary_meta, "sample_id", "primary-taxon metadata")

query_sites <- sites[genomic_seqid %in% canonical & ag3_phased_panel_assessable == TRUE]
if (!nrow(query_sites)) fail_step(STEP, "No frozen target sites are assessable on configured canonical contigs.")

# The local accessibility handoff must cover exactly the canonical frozen targets.
missing_access <- setdiff(query_sites$population_site_id, access_dt$population_site_id)
extra_access <- setdiff(access_dt$population_site_id, query_sites$population_site_id)
if (length(missing_access) || length(extra_access)) {
  fail_step(STEP, paste0(
    "Local accessibility target IDs do not match the frozen canonical panel. Missing=",
    paste(utils::head(missing_access, 20L), collapse = ", "),
    "; extra=", paste(utils::head(extra_access, 20L), collapse = ", ")
  ))
}
access_dt[, accessibility_complete := normalize_flag(accessibility_complete)]
if (anyNA(access_dt$accessibility_complete)) fail_step(STEP, "accessibility_complete contains non-boolean values.")
if (any(access_dt$accessibility_records_expected != 23L, na.rm = TRUE)) {
  fail_step(STEP, "Every targeted accessibility record must expect exactly 23 positions.")
}

# Optional sex metadata is used only as a fallback when an X-chromosome sample
# has no variant calls in the very small queried target windows.
sex_col <- resolve_column(primary_meta, c("sex_call", "sex", "sex_calling"), "sex column", required = FALSE)
known_male <- rep(FALSE, nrow(primary_meta))
if (!is.na(sex_col)) {
  sx <- tolower(trimws(as.character(primary_meta[[sex_col]])))
  known_male <- sx %in% c("m", "male", "1", "xy")
}
primary_meta[, known_male_for_x_fallback := known_male]

fixed_gt_cols <- c("variant_id", "contig", "position", "ref", "alt", "filter")
state_exact <- list()
state_pam <- list()
state_proto <- list()
variant_store <- list()
contig_qc <- list()
master_samples <- NULL

for (contig in canonical) {
  ss <- query_sites[genomic_seqid == contig]
  if (!nrow(ss)) next

  gt_file <- file.path(local_dir, "haplotypes", paste0(contig, ".phased_genotypes.tsv"))
  if (!file.exists(gt_file)) {
    fail_step(STEP, paste0(
      "Missing targeted phased-genotype cache for ", contig, ": ", gt_file, "\n",
      "Rerun scripts/acquire_ag3_targeted.py. It will reuse completed cached contigs."
    ))
  }

  log_step(STEP, paste0("Reading local targeted phased data for ", contig, " (", nrow(ss), " frozen target sites)"))

  # Read header first. This preserves the sample axis even when a target set has
  # zero variant records and the file therefore contains only the header row.
  hdr <- fread(gt_file, sep = "\t", nrows = 0L, check.names = FALSE)
  header_names <- names(hdr)
  if (!all(fixed_gt_cols %in% header_names)) {
    fail_step(STEP, paste0(
      "Malformed targeted genotype file for ", contig, ". Missing fixed columns: ",
      paste(setdiff(fixed_gt_cols, header_names), collapse = ", ")
    ))
  }
  panel_samples <- setdiff(header_names, fixed_gt_cols)
  if (!length(panel_samples)) fail_step(STEP, paste0("No phased-panel sample columns found in ", gt_file, "."))
  if (anyDuplicated(panel_samples)) fail_step(STEP, paste0("Duplicate sample columns found in ", gt_file, "."))

  available <- panel_samples[panel_samples %in% primary_meta$sample_id]
  if (!length(available)) fail_step(STEP, paste0("No primary Ag3 samples matched the local phased-panel file for ", contig, "."))

  if (is.null(master_samples)) {
    master_samples <- available
  } else {
    if (!setequal(master_samples, available)) {
      fail_step(STEP, paste0(
        "Local phased-panel sample membership differs across contigs. On ", contig,
        ", missing relative to first contig: ", paste(utils::head(setdiff(master_samples, available), 10L), collapse = ", "),
        "; extra: ", paste(utils::head(setdiff(available, master_samples), 10L), collapse = ", ")
      ))
    }
    available <- master_samples
  }

  x <- fread(gt_file, sep = "\t", check.names = FALSE)
  if (!all(header_names %in% names(x)) || !identical(names(x), header_names)) {
    fail_step(STEP, paste0("Header changed while reading local genotype file for ", contig, "."))
  }

  if (nrow(x)) {
    x[, position := as.integer(position)]
    if (anyNA(x$position) || any(x$position < 1L)) fail_step(STEP, paste0("Invalid variant position in ", gt_file, "."))
    if (any(as.character(x$contig) != contig)) fail_step(STEP, paste0("Contig mismatch inside ", gt_file, "."))
    if (anyDuplicated(x$variant_id)) fail_step(STEP, paste0("Duplicate variant_id in ", gt_file, "."))
    if (any(nchar(as.character(x$ref)) != 1L)) fail_step(STEP, paste0("Non-SNP REF allele in local phased data for ", contig, "."))
    alt_bad <- vapply(strsplit(as.character(x$alt), ",", fixed = TRUE), function(z) any(nchar(z) != 1L), logical(1))
    if (any(alt_bad)) fail_step(STEP, paste0("Non-SNP ALT allele in local phased data for ", contig, "."))
  }

  if (nrow(x)) {
    gt <- as.matrix(x[, ..available])
    storage.mode(gt) <- "character"
    colnames(gt) <- available
  } else {
    gt <- matrix(character(0), nrow = 0L, ncol = length(available), dimnames = list(NULL, available))
  }
  if (!identical(colnames(gt), master_samples)) gt <- gt[, master_samples, drop = FALSE]

  var_dt <- data.table::data.table(
    variant_index = seq_len(nrow(x)),
    contig = if (nrow(x)) as.character(x$contig) else character(),
    position = if (nrow(x)) as.integer(x$position) else integer(),
    ref = if (nrow(x)) as.character(x$ref) else character(),
    alt = if (nrow(x)) as.character(x$alt) else character(),
    filter = if (nrow(x)) as.character(x$filter) else character()
  )
  pos <- var_dt$position
  if (nrow(gt) != nrow(var_dt)) fail_step(STEP, paste0("Variant table and genotype matrix row counts differ for ", contig, "."))
  variant_store[[contig]] <- list(variants = var_dt, gt = gt, sample_ids = colnames(gt))

  # Determine whether each sample has a real second haplotype on this contig.
  # We inspect the acquired phased calls. When no target variant gives ploidy
  # evidence on X, known male metadata is used as a conservative fallback.
  second_haplotype_exists <- rep(TRUE, length(master_samples))
  names(second_haplotype_exists) <- master_samples
  mixed_ploidy_samples <- character()

  if (nrow(gt)) {
    for (j in seq_along(master_samples)) {
      g <- as.character(gt[, j])
      parsed <- parse_gt_haplotypes(g, require_phased = TRUE)
      any_h1 <- any(!is.na(parsed$a1))
      any_h2 <- any(!is.na(parsed$a2))
      if (any_h1 && !any_h2) second_haplotype_exists[[j]] <- FALSE
      if (any_h2 && any(is.na(parsed$a2) & !is.na(parsed$a1))) {
        mixed_ploidy_samples <- c(mixed_ploidy_samples, master_samples[[j]])
      }
    }
  }

  if (identical(contig, "X")) {
    no_ploidy_evidence <- rep(nrow(gt) == 0L, length(master_samples))
    names(no_ploidy_evidence) <- master_samples
    if (nrow(gt)) {
      for (j in seq_along(master_samples)) {
        parsed <- parse_gt_haplotypes(gt[, j], require_phased = TRUE)
        no_ploidy_evidence[[j]] <- !any(!is.na(parsed$a1)) && !any(!is.na(parsed$a2))
      }
    }
    male_map <- primary_meta$known_male_for_x_fallback[match(master_samples, primary_meta$sample_id)]
    male_map[is.na(male_map)] <- FALSE
    fallback_haploid <- no_ploidy_evidence & male_map
    second_haplotype_exists[fallback_haploid] <- FALSE
  }

  if (length(mixed_ploidy_samples)) {
    log_warning(STEP, paste0(
      "Mixed one-/two-haplotype call representation detected on ", contig, " for ",
      length(unique(mixed_ploidy_samples)), " sample(s). A second haplotype is retained where any called h2 evidence exists; missing h2 calls remain NA."
    ))
  }

  hap_names <- make_haplotype_names(master_samples)
  exact_mat <- matrix(TRUE, nrow = nrow(ss), ncol = length(hap_names), dimnames = list(ss$population_site_id, hap_names))
  pam_mat <- exact_mat
  proto_mat <- exact_mat

  haploid_samples <- names(second_haplotype_exists)[!second_haplotype_exists]
  if (length(haploid_samples)) {
    h2_cols <- paste0(haploid_samples, "|h2")
    exact_mat[, h2_cols] <- NA
    pam_mat[, h2_cols] <- NA
    proto_mat[, h2_cols] <- NA
  }

  for (k in seq_len(nrow(ss))) {
    s <- ss[k]
    idx <- which(pos >= as.integer(s$genomic_start) & pos <= as.integer(s$genomic_end))
    if (!length(idx)) next

    pam_pos <- if (identical(as.character(s$guide_genomic_strand), "+")) {
      (as.integer(s$genomic_end) - 2L):as.integer(s$genomic_end)
    } else {
      as.integer(s$genomic_start):(as.integer(s$genomic_start) + 2L)
    }

    for (j in idx) {
      parsed <- parse_gt_haplotypes(gt[j, ], require_phased = TRUE)
      r1 <- allele_is_reference(parsed$a1)
      r2 <- allele_is_reference(parsed$a2)
      new_state <- interleave_haplotypes(r1, r2)
      exact_mat[k, ] <- combine_exact_state(exact_mat[k, ], new_state)
      if (pos[[j]] %in% pam_pos) {
        pam_mat[k, ] <- combine_exact_state(pam_mat[k, ], new_state)
      } else {
        proto_mat[k, ] <- combine_exact_state(proto_mat[k, ], new_state)
      }
    }
  }

  state_exact[[contig]] <- exact_mat
  state_pam[[contig]] <- pam_mat
  state_proto[[contig]] <- proto_mat

  a <- access_dt[population_site_id %in% ss$population_site_id]
  if (nrow(a) != nrow(ss)) fail_step(STEP, paste0("Accessibility row count mismatch on ", contig, "."))
  contig_qc[[contig]] <- data.table::data.table(
    contig = contig,
    n_target_sites = nrow(ss),
    n_phased_panel_samples = length(master_samples),
    n_variant_records_in_targets = nrow(var_dt),
    n_inferred_haploid_samples = length(haploid_samples),
    n_sites_with_complete_accessibility = sum(a$accessibility_complete == TRUE, na.rm = TRUE),
    n_sites_with_incomplete_accessibility = sum(a$accessibility_complete != TRUE | is.na(a$accessibility_complete))
  )
}

if (is.null(master_samples) || !length(master_samples)) fail_step(STEP, "No phased-panel samples were recovered from local targeted files.")
if (!length(state_exact)) fail_step(STEP, "No target haplotype-state matrices were produced.")

# Enforce one identical haplotype axis across all contigs before row-binding.
expected_haps <- make_haplotype_names(master_samples)
for (nm in names(state_exact)) {
  if (!identical(colnames(state_exact[[nm]]), expected_haps)) stop("Haplotype column mismatch in exact-state matrix for ", nm, call. = FALSE)
  if (!identical(colnames(state_pam[[nm]]), expected_haps)) stop("Haplotype column mismatch in PAM-state matrix for ", nm, call. = FALSE)
  if (!identical(colnames(state_proto[[nm]]), expected_haps)) stop("Haplotype column mismatch in protospacer-state matrix for ", nm, call. = FALSE)
}

exact_all <- do.call(rbind, state_exact)
pam_all <- do.call(rbind, state_pam)
proto_all <- do.call(rbind, state_proto)
row_order <- match(query_sites$population_site_id, rownames(exact_all))
if (anyNA(row_order)) {
  miss <- query_sites$population_site_id[is.na(row_order)]
  fail_step(STEP, paste0("State matrices are missing queried target site(s): ", paste(utils::head(miss, 20L), collapse = ", ")))
}
exact_all <- exact_all[row_order, , drop = FALSE]
pam_all <- pam_all[row_order, , drop = FALSE]
proto_all <- proto_all[row_order, , drop = FALSE]

states <- list(
  exact_23bp = exact_all,
  pam_intact = pam_all,
  protospacer_exact = proto_all,
  sample_ids = master_samples,
  haplotype_names = expected_haps,
  reference_genome = "AgamP4",
  phasing_panel = phasing_panel,
  source_mode = "LOCAL_TARGETED_MALARIAGEN_DATA_API",
  note = paste(
    "Rows are population_site_id; columns are phased sample haplotypes sample_id|h1 or sample_id|h2.",
    "TRUE means reference-identical at all observed biallelic SNPs in the queried component;",
    "FALSE means at least one non-reference allele; NA means state cannot be established from called haplotypes."
  )
)
atomic_saveRDS(states, "data_processed/02_target_haplotype_states.rds", compress = "xz")
atomic_saveRDS(variant_store, "data_processed/02_target_variant_genotypes.rds", compress = "xz")

# Reorder local accessibility to the frozen target order and persist the same
# downstream interface used by the original workflow.
access_dt <- access_dt[match(query_sites$population_site_id, population_site_id)]
if (anyNA(access_dt$population_site_id)) fail_step(STEP, "Failed to deterministically reorder accessibility rows.")
atomic_fwrite(access_dt, "data_processed/02_target_accessibility.csv")

# Re-freeze metadata/population manifest to exact phased-panel sample membership.
analysis_meta <- meta[sample_id %in% master_samples]
analysis_meta[, vcf_sample_order := match(sample_id, master_samples)]
data.table::setorder(analysis_meta, vcf_sample_order)
if (nrow(analysis_meta) != length(master_samples)) {
  missing_meta <- setdiff(master_samples, analysis_meta$sample_id)
  fail_step(STEP, paste0("Not every phased-panel sample maps to Step 01 metadata. Missing: ", paste(utils::head(missing_meta, 20L), collapse = ", ")))
}
if (!identical(analysis_meta$sample_id, master_samples)) fail_step(STEP, "Analysis metadata could not be aligned exactly to phased-panel sample order.")

analysis_pops <- analysis_meta[primary_taxon == TRUE & nonempty_string(population_id), .(
  n_samples = data.table::uniqueN(sample_id),
  n_locations = data.table::uniqueN(location[nonempty_string(location)]),
  year_min = safe_min(year_numeric),
  year_max = safe_max(year_numeric)
), by = .(population_id, country, taxon)]
analysis_pops[, primary_population_n20 := n_samples >= primary_min_n]
for (nmin in sensitivity_ns) analysis_pops[, (paste0("sensitivity_population_n", nmin)) := n_samples >= nmin]
for (nmin in c(10L, 30L)) {
  nm <- paste0("sensitivity_population_n", nmin)
  if (!nm %in% names(analysis_pops)) analysis_pops[, (nm) := n_samples >= nmin]
}
data.table::setorder(analysis_pops, country, taxon, population_id)
atomic_fwrite(analysis_meta, "data_processed/02_analysis_sample_metadata.csv")
atomic_fwrite(analysis_pops, "data_processed/02_analysis_population_manifest.csv")

unresolved <- sites[!genomic_seqid %in% canonical | ag3_phased_panel_assessable == FALSE, .(
  population_site_id, gene_id, genomic_seqid, genomic_start, genomic_end,
  status = "UNRESOLVED_CONTIG_NOT_IN_AG3_GAMB_COLU_PHASED_PANEL"
)]
atomic_fwrite(unresolved, "data_processed/02_unresolved_unplaced_targets.csv")

contig_qc_dt <- data.table::rbindlist(contig_qc, fill = TRUE, use.names = TRUE)
atomic_fwrite(contig_qc_dt, "data_processed/02_contig_extraction_qc.csv")
qc <- data.table::data.table(
  metric = c(
    "Canonical sites queried", "Unplaced sites unresolved", "Genes with unresolved unplaced sites",
    "Primary samples in phased panel", "Country-taxon populations n>=20 after phased-panel intersection",
    "Sites with complete 23-bp accessibility", "Sites with incomplete accessibility",
    "R network requests in Step 02", "Chromosome-scale VCFs downloaded by Step 02"
  ),
  value = c(
    nrow(query_sites), nrow(unresolved), data.table::uniqueN(unresolved$gene_id),
    length(master_samples), sum(analysis_pops$primary_population_n20),
    sum(access_dt$accessibility_complete == TRUE, na.rm = TRUE),
    sum(access_dt$accessibility_complete != TRUE | is.na(access_dt$accessibility_complete)),
    0L, 0L
  )
)
atomic_fwrite(qc, "data_processed/02_target_extraction_qc.csv")

checksum_inputs <- c(
  complete_file, access_file,
  unlist(lapply(names(state_exact), function(contig) file.path(local_dir, "haplotypes", paste0(contig, ".phased_genotypes.tsv")))),
  if (file.exists(acq_qc_file)) acq_qc_file else character(),
  "data_processed/02_target_haplotype_states.rds", "data_processed/02_target_variant_genotypes.rds",
  "data_processed/02_target_accessibility.csv", "data_processed/02_analysis_sample_metadata.csv",
  "data_processed/02_analysis_population_manifest.csv", "data_processed/02_unresolved_unplaced_targets.csv"
)
write_checksum(checksum_inputs, "logs/02_checksums.tsv")
write_session_info("logs/02_sessionInfo.txt")
print(qc)
log_step(STEP, "Local targeted Ag3.0 phased-data conversion completed successfully")
