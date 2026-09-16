# WINDOWS-SAFE PHASE-2 ACQUISITION B
# Query the public Sanger Phase-2 SHAPEIT VCFs directly with Bioconductor
# Rsamtools/VariantAnnotation. Only the frozen small target windows are read.
# No chromosome-scale VCF is downloaded and no Python/pysam is required.

if (!file.exists("R/helpers.R")) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    p <- sub("^--file=", "", arg[[1]])
    candidate <- dirname(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(candidate, "R", "helpers.R"))) setwd(candidate)
  }
}
if (!file.exists("R/helpers.R")) stop("Run this script from the project root.", call. = FALSE)
source("R/helpers.R")
set_project_root()
require_packages(c("data.table", "yaml", "curl", "jsonlite"), bioc = c("Rsamtools", "VariantAnnotation", "GenomicRanges", "IRanges"))
suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(IRanges)
})

STEP <- "00b"
log_step(STEP, "Starting Windows-native Phase-2 phased VCF acquisition via remote Rsamtools TabixFile")
log_step(STEP, "Only frozen target intervals are queried; whole chromosome VCF download is prohibited")

cfg <- read_config()
local_dir <- as.character(cfg_get(cfg, "phase2_acquisition", "output_dir", default = "data_raw/phase2_local"))
canonical <- as.character(cfg_get(cfg, "canonical_contigs", required = TRUE))
prep_marker <- file.path(local_dir, "PREPARATION_COMPLETE.ok")
manifest_file <- file.path(local_dir, "haplotypes", "remote_vcf_manifest.csv")
meta_file <- file.path(local_dir, "metadata", "sample_metadata.csv")
access_file <- file.path(local_dir, "accessibility", "target_accessibility.tsv")
plan_file <- file.path(local_dir, "request_plan.csv")
assert_files(c(
  "data_processed/00_frozen_target_sites.csv",
  prep_marker, manifest_file, meta_file, access_file, plan_file
))

sites <- fread("data_processed/00_frozen_target_sites.csv", na.strings = c("", "NA", "NaN"))
meta <- fread(meta_file, na.strings = c("", "NA", "NaN", "nan", "None"))
remote <- fread(manifest_file, na.strings = c("", "NA", "NaN"))
plan <- fread(plan_file, na.strings = c("", "NA", "NaN"))
assert_columns(sites, c("population_site_id", "genomic_seqid", "genomic_start", "genomic_end"), "frozen sites")
assert_columns(meta, c("sample_id", "primary_taxon"), "Phase-2 raw metadata")
assert_columns(remote, c("contig", "vcf_url", "index_url"), "Phase-2 remote VCF manifest")
assert_columns(plan, c("kind", "contig", "start", "end"), "Phase-2 request plan")
assert_unique(remote, "contig", "Phase-2 remote VCF manifest")

meta[, primary_taxon := normalize_flag(primary_taxon)]
meta[is.na(primary_taxon), primary_taxon := FALSE]
primary_meta <- meta[primary_taxon == TRUE]
if (!nrow(primary_meta)) fail_step(STEP, "No primary gambiae/coluzzii samples in prepared Phase-2 metadata.")
assert_unique(primary_meta, "sample_id", "Phase-2 primary metadata")
primary_ids <- as.character(primary_meta$sample_id)

# Explicit one-byte HTTP range probe. If the server does not return 206, abort
# rather than risk a full-object transfer.
range_probe <- function(url) {
  # HEAD-only fail-closed check. We never issue a one-byte GET probe here because
  # a server that ignores Range could otherwise return the entire multi-GB VCF.
  h <- curl::new_handle()
  curl::handle_setopt(h, nobody = TRUE, followlocation = TRUE)
  res <- tryCatch(curl::curl_fetch_memory(url, handle = h), error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, status = NA_integer_, message = conditionMessage(res), msg = conditionMessage(res)))
  hdr <- tryCatch(curl::parse_headers_list(res$headers), error = function(e) list())
  ar <- tolower(as.character(hdr[["accept-ranges"]] %||% hdr[["Accept-Ranges"]] %||% ""))
  ok <- as.integer(res$status_code) < 400L && grepl("bytes", ar, fixed = TRUE)
  list(ok = ok, status = as.integer(res$status_code), message = paste0("Accept-Ranges=", ar), msg = paste0("Accept-Ranges=", ar))
}

retry <- function(fun, label, attempts = 3L, sleep_sec = c(0, 2, 5)) {
  errs <- character()
  for (i in seq_len(attempts)) {
    if (sleep_sec[[min(i, length(sleep_sec))]] > 0) Sys.sleep(sleep_sec[[min(i, length(sleep_sec))]])
    ans <- tryCatch(fun(), error = function(e) e)
    if (!inherits(ans, "error")) return(ans)
    errs <- c(errs, paste0("attempt ", i, ": ", conditionMessage(ans)))
  }
  fail_step(STEP, paste0(label, " failed after ", attempts, " attempts. ", paste(errs, collapse = " | ")))
}

fixed_cols <- c("variant_id", "contig", "position", "ref", "alt", "filter")
master_samples <- NULL
contig_qc <- list()
out_hap <- file.path(local_dir, "haplotypes")
dir.create(out_hap, recursive = TRUE, showWarnings = FALSE)

for (contig in canonical) {
  ss <- sites[genomic_seqid == contig]
  if (!nrow(ss)) next
  rr <- remote[as.character(remote$contig) == contig]
  if (nrow(rr) != 1L) fail_step(STEP, paste0("Expected exactly one remote Phase-2 VCF for ", contig, "."))
  url <- as.character(rr$vcf_url[[1]])
  idx <- as.character(rr$index_url[[1]])

  probe <- range_probe(url)
  if (!isTRUE(probe$ok)) {
    fail_step(STEP, paste0(
      "Public Phase-2 VCF does not confirm HTTP byte-range serving for ", contig,
      " (HTTP ", probe$status %||% "NA", "). Refusing any possible whole-file transfer. ", probe$message
    ))
  }

  tf <- Rsamtools::TabixFile(url, index = idx)
  header <- retry(function() VariantAnnotation::scanVcfHeader(tf), paste0("Remote VCF header for ", contig))
  panel_samples <- as.character(VariantAnnotation::samples(header))
  if (!length(panel_samples)) fail_step(STEP, paste0("No sample IDs in Phase-2 phased VCF header for ", contig, "."))
  available <- panel_samples[panel_samples %in% primary_ids]
  if (length(available) < 100L) fail_step(STEP, paste0(contig, ": only ", length(available), " primary samples intersect the phased VCF header."))

  if (is.null(master_samples)) {
    master_samples <- available
  } else {
    if (!setequal(master_samples, available)) {
      fail_step(STEP, paste0(
        "Phase-2 phased sample membership differs across contigs at ", contig,
        ". Missing=", paste(head(setdiff(master_samples, available), 10L), collapse = ","),
        "; extra=", paste(head(setdiff(available, master_samples), 10L), collapse = ",")
      ))
    }
    available <- master_samples
  }

  wp <- plan[kind == "haplotype" & as.character(plan$contig) == contig]
  if (!nrow(wp)) fail_step(STEP, paste0("No haplotype request windows for ", contig, "."))
  gr <- GenomicRanges::GRanges(
    seqnames = rep(contig, nrow(wp)),
    ranges = IRanges::IRanges(start = as.integer(wp$start), end = as.integer(wp$end))
  )
  param <- VariantAnnotation::ScanVcfParam(
    which = gr,
    samples = available,
    fixed = c("ALT", "FILTER"),
    info = NA_character_,
    geno = "GT"
  )

  log_step(STEP, paste0(contig, ": querying ", nrow(wp), " small remote target window(s) for ", length(available), " phased samples"))
  vcf <- retry(function() VariantAnnotation::readVcf(tf, genome = "AgamP4", param = param), paste0("Remote tabix query for ", contig))

  if (nrow(vcf)) {
    rrng <- VariantAnnotation::rowRanges(vcf)
    refv <- as.character(VariantAnnotation::ref(vcf))
    altv <- vapply(VariantAnnotation::alt(vcf), function(z) {
      if (length(z) == 1L) as.character(z[[1]]) else paste(as.character(z), collapse = ",")
    }, character(1))
    filt <- as.character(VariantAnnotation::fixed(vcf)$FILTER)
    posv <- as.integer(start(rrng))
    contigv <- as.character(seqnames(rrng))
    gt <- VariantAnnotation::geno(vcf)$GT
    if (is.null(gt)) fail_step(STEP, paste0("GT field missing from Phase-2 phased VCF for ", contig, "."))
    gt <- as.matrix(gt)
    if (!identical(colnames(gt), available)) {
      m <- match(available, colnames(gt))
      if (anyNA(m)) fail_step(STEP, paste0("Requested sample columns missing from returned GT matrix for ", contig, "."))
      gt <- gt[, m, drop = FALSE]
    }

    # Discovery is intentionally restricted to biallelic SNPs on the phased scaffold.
    keep <- nchar(refv) == 1L & !grepl(",", altv, fixed = TRUE) & nchar(altv) == 1L
    refv <- refv[keep]; altv <- altv[keep]; filt <- filt[keep]
    posv <- posv[keep]; contigv <- contigv[keep]; gt <- gt[keep, , drop = FALSE]

    key <- paste(contigv, posv, refv, altv, sep = ":")
    if (anyDuplicated(key)) {
      # Overlapping query windows can cause duplicate records. Deduplicate only
      # if every duplicate row carries identical GT values.
      dupkeys <- unique(key[duplicated(key)])
      for (kk in dupkeys) {
        ii <- which(key == kk)
        if (length(ii) > 1L && !all(vapply(ii[-1L], function(j) identical(gt[ii[[1]], ], gt[j, ]), logical(1)))) {
          fail_step(STEP, paste0("Non-identical duplicate VCF record returned for ", kk, "."))
        }
      }
      keep1 <- !duplicated(key)
      key <- key[keep1]; refv <- refv[keep1]; altv <- altv[keep1]; filt <- filt[keep1]
      posv <- posv[keep1]; contigv <- contigv[keep1]; gt <- gt[keep1, , drop = FALSE]
    }

    ord <- order(posv, refv, altv)
    key <- key[ord]; refv <- refv[ord]; altv <- altv[ord]; filt <- filt[ord]
    posv <- posv[ord]; contigv <- contigv[ord]; gt <- gt[ord, , drop = FALSE]

    fixed <- data.table(
      variant_id = key,
      contig = contigv,
      position = posv,
      ref = refv,
      alt = altv,
      filter = filt
    )
    gtdt <- as.data.table(gt)
    setnames(gtdt, available)
    out <- cbind(fixed, gtdt)
  } else {
    out <- as.data.table(setNames(replicate(length(fixed_cols) + length(available), character(), simplify = FALSE), c(fixed_cols, available)))
  }

  outfile <- file.path(out_hap, paste0(contig, ".phased_targets.tsv"))
  atomic_fwrite(out, outfile, sep = "\t")
  contig_qc[[contig]] <- data.table(
    contig = contig,
    n_target_sites = nrow(ss),
    n_request_windows = nrow(wp),
    n_phased_samples = length(available),
    n_biallelic_snp_records_in_targets = nrow(out),
    local_cache_bytes = file.info(outfile)$size,
    whole_vcf_downloaded = FALSE
  )
  log_step(STEP, paste0(contig, ": cached ", nrow(out), " biallelic target-overlap SNP record(s); whole VCF downloaded = FALSE"))
}

if (is.null(master_samples) || !length(master_samples)) fail_step(STEP, "No Phase-2 phased samples were acquired.")
phase_meta <- meta[match(master_samples, sample_id)]
if (nrow(phase_meta) != length(master_samples) || anyNA(phase_meta$sample_id)) {
  fail_step(STEP, "Not every phased VCF sample maps to prepared Phase-2 metadata.")
}
if (!identical(as.character(phase_meta$sample_id), master_samples)) fail_step(STEP, "Phase-2 phased metadata order mismatch.")
atomic_fwrite(phase_meta, file.path(local_dir, "metadata", "phased_sample_metadata.csv"))
qc <- rbindlist(contig_qc, use.names = TRUE, fill = TRUE)
atomic_fwrite(qc, file.path(local_dir, "acquisition_qc.tsv"), sep = "\t")

manifest <- list(
  schema_version = "phase2-sanger-windows-rsamtools-v3",
  discovery_dataset = "Ag1000G Phase 2 AR1",
  access_engine = "R/Bioconductor Rsamtools + VariantAnnotation remote TabixFile",
  n_frozen_targets = nrow(sites),
  n_canonical_targets = sites[genomic_seqid %in% canonical, .N],
  n_phased_primary_samples = length(master_samples),
  whole_chromosome_vcf_downloads = 0L,
  contig_qc = as.data.frame(qc),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(manifest, file.path(local_dir, "acquisition_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
atomic_writeLines("phase2-sanger-windows-rsamtools-v3", file.path(local_dir, "ACQUISITION_COMPLETE.ok"))
write_session_info("logs/00b_phase2_acquisition_sessionInfo.txt")
log_step(STEP, paste0("Phase-2 remote phased acquisition completed successfully for ", length(master_samples), " primary samples"))
