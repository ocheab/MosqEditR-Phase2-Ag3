# WINDOWS-SAFE AG3 EXTERNAL-VALIDATION ACQUISITION
# Public Sanger per-sample *all-sites* VCFs are queried remotely with
# Bioconductor Rsamtools/VariantAnnotation. Whole ~GB VCFs are never downloaded.
# Scientific role: external validation only, after the Phase-2 target lock exists.

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
require_packages(c("data.table", "yaml", "curl", "jsonlite", "digest"), bioc = c("Rsamtools", "VariantAnnotation", "SummarizedExperiment", "GenomicRanges", "IRanges"))
suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(IRanges)
})

STEP <- "05c"
args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
final_mode <- "--final" %in% args
pilot_mode <- "--pilot" %in% args || !final_mode
if (final_mode && "--pilot" %in% args) fail_step(STEP, "Use only one of --pilot or --final.")
run_mode <- if (final_mode) "final" else "pilot"

log_step(STEP, paste0("Starting Ag3 public-Sanger external acquisition in ", toupper(run_mode), " mode via Rsamtools"))
log_step(STEP, "Ag3 is external validation only; Phase-2 scores/ranks/thresholds are immutable")

cfg <- read_config()
ev <- cfg_get(cfg, "external_validation", default = list())
if (!isTRUE(as.logical(ev$enabled %||% TRUE))) fail_step(STEP, "external_validation.enabled is FALSE.")
if (!isTRUE(as.logical(ev$require_no_discovery_retuning %||% TRUE))) fail_step(STEP, "Validation isolation must remain enabled.")
lock_file <- "data_processed/05b_external_validation_targets.csv"
assert_file(lock_file)
targets <- fread(lock_file, na.strings = c("", "NA", "NaN"))
assert_columns(targets, c(
  "population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end",
  "guide_genomic_strand", "selection_reason", "selection_locked_before_ag3"
), "external target lock")
assert_unique(targets, "population_site_id", "external target lock")
targets[, selection_locked_before_ag3 := normalize_flag(selection_locked_before_ag3)]
if (any(targets$selection_locked_before_ag3 != TRUE | is.na(targets$selection_locked_before_ag3))) {
  fail_step(STEP, "Target lock is not immutable for every external-validation site.")
}
canonical <- as.character(cfg_get(cfg, "canonical_contigs", required = TRUE))
if (any(!targets$genomic_seqid %in% canonical)) fail_step(STEP, "External target lock contains noncanonical contigs.")
if (any((targets$genomic_end - targets$genomic_start + 1L) != 23L)) fail_step(STEP, "External target lock contains a non-23-bp target.")
lock_hash <- digest::digest(file = lock_file, algo = "sha256", serialize = FALSE)

curl_text_lines <- function(url) {
  res <- curl::curl_fetch_memory(url)
  if (as.integer(res$status_code) >= 400L) stop("HTTP ", res$status_code, " for ", url, call. = FALSE)
  x <- strsplit(rawToChar(res$content), "\\r?\\n", perl = TRUE)[[1L]]
  trimws(x[nzchar(trimws(x))])
}
normalize_taxon_local <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- sub("^an\\.\\s*", "", z)
  z <- sub("^anopheles\\s+", "", z)
  z[z %in% c("gambiae s.s.", "gambiae_ss", "s")] <- "gambiae"
  z[z %in% c("m")] <- "coluzzii"
  z
}
sample_id_from_feature <- function(x) {
  b <- basename(sub("/+$", "", x))
  b <- sub("\\.gatk\\.zarr\\.zip$", "", b, ignore.case = TRUE)
  b <- sub("\\.zarr\\.zip$", "", b, ignore.case = TRUE)
  b <- sub("\\.vcf\\.gz$", "", b, ignore.case = TRUE)
  sub("\\..*$", "", b)
}
features <- curl_text_lines(as.character(ev$features_url))
labels <- curl_text_lines(as.character(ev$labels_url))
if (length(features) != length(labels)) fail_step(STEP, paste0("Public Ag3 feature/label manifest lengths differ: ", length(features), " vs ", length(labels)))
if (length(features) < 1000L) fail_step(STEP, paste0("Ag3 public manifest unexpectedly contains only ", length(features), " rows."))
md <- data.table(feature_url = features, taxon_raw = labels)
md[, sample_id := vapply(feature_url, sample_id_from_feature, character(1))]
md[, taxon := normalize_taxon_local(taxon_raw)]
if (anyDuplicated(md$sample_id)) fail_step(STEP, "Duplicate sample IDs in Ag3 public manifest.")
base_url <- sub("/+$", "", as.character(ev$public_vcf_base_url %||% "https://cog.sanger.ac.uk/vo_agam_output"))
md[, vcf_url := paste0(base_url, "/", sample_id, ".vcf.gz")]
primary_taxa <- tolower(as.character(ev$primary_taxa %||% c("gambiae", "coluzzii")))
md <- md[taxon %in% primary_taxa]
if (!nrow(md)) fail_step(STEP, "No configured primary taxa remain in Ag3 public manifest.")

per_taxon <- if (final_mode) as.integer(ev$final_max_samples_per_taxon %||% 250L) else as.integer(ev$pilot_samples_per_taxon %||% 25L)
selected <- rbindlist(lapply(primary_taxa, function(tx) {
  d <- md[taxon == tx][order(sample_id)]
  if (!nrow(d)) fail_step(STEP, paste0("No Ag3 samples labelled ", tx, "."))
  if (per_taxon > 0L) d <- head(d, per_taxon)
  d
}), use.names = TRUE, fill = TRUE)

out_base <- as.character(ev$output_dir %||% "data_raw/ag3_external")
out_dir <- file.path(out_base, run_mode)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (final_mode && !file.exists(file.path(out_base, "pilot", "ACQUISITION_COMPLETE.ok"))) {
  fail_step(STEP, "--final is blocked until the small Ag3 pilot completes successfully.")
}
atomic_fwrite(selected, file.path(out_dir, "selected_samples.tsv"), sep = "\t")
atomic_fwrite(targets, file.path(out_dir, "locked_targets.tsv"), sep = "\t")

# For all-sites VCFs, do NOT merge across gaps: every extra base would create a
# record and unnecessary transfer. Merge only overlapping/adjacent target spans.
merge_exact_windows <- function(d) {
  ans <- list(); k <- 0L
  for (cc in unique(as.character(d$genomic_seqid))) {
    z <- d[genomic_seqid == cc][order(genomic_start, genomic_end)]
    s0 <- as.integer(z$genomic_start[[1L]]); e0 <- as.integer(z$genomic_end[[1L]])
    if (nrow(z) > 1L) for (i in 2:nrow(z)) {
      s <- as.integer(z$genomic_start[[i]]); e <- as.integer(z$genomic_end[[i]])
      if (s <= e0 + 1L) e0 <- max(e0, e) else {
        k <- k + 1L; ans[[k]] <- data.table(contig = cc, start = s0, end = e0)
        s0 <- s; e0 <- e
      }
    }
    k <- k + 1L; ans[[k]] <- data.table(contig = cc, start = s0, end = e0)
  }
  rbindlist(ans)
}
windows <- merge_exact_windows(targets)
window_bp <- sum(windows$end - windows$start + 1L)
log_step(STEP, paste0("Locked external targets: ", nrow(targets), " sites / ", uniqueN(targets$gene_id), " genes; SHA256=", substr(lock_hash, 1L, 16L), "..."))
log_step(STEP, paste0("Ag3 ", run_mode, ": selected ", nrow(selected), " samples; ", paste(selected[, .N, by = taxon][, paste0(taxon, "=", N)], collapse = ", ")))
log_step(STEP, paste0("Per-sample exact-window plan: ", nrow(windows), " window(s) / ", format(window_bp, big.mark = ","), " logical bp"))
log_step(STEP, "Bandwidth guard: public all-sites VCFs remain remote; only tabix HTTP ranges are read")
if (dry_run) {
  log_step(STEP, "DRY RUN complete: small public manifests read; no genomic VCF opened")
  quit(save = "no", status = 0L)
}

range_probe <- function(url) {
  # Fail-safe byte-range verification.
  #
  # HEAD is not reliable for this Sanger object endpoint on all Windows/libcurl
  # combinations, and may omit Accept-Ranges even when byte-range GETs work.
  # Instead request exactly one byte and require HTTP 206 + Content-Range.
  # maxfilesize=1024 makes this fail closed: if the server ignores Range and
  # advertises the multi-GB object length, libcurl aborts before downloading it.
  h <- curl::new_handle()
  curl::handle_setopt(
    h,
    range = "0-0",
    followlocation = TRUE,
    maxfilesize = 1024,
    connecttimeout = 20,
    timeout = 30
  )
  res <- tryCatch(curl::curl_fetch_memory(url, handle = h), error = function(e) e)
  if (inherits(res, "error")) {
    return(list(
      ok = FALSE,
      status = NA_integer_,
      message = conditionMessage(res),
      msg = conditionMessage(res),
      content_range = NA_character_
    ))
  }

  hdr <- tryCatch(curl::parse_headers_list(res$headers), error = function(e) list())
  # Header names returned by parse_headers_list() are normally lower-case, but
  # tolerate either representation.
  cr <- as.character(hdr[["content-range"]] %||% hdr[["Content-Range"]] %||% "")
  status <- as.integer(res$status_code)
  one_byte <- length(res$content) <= 1L
  ok <- identical(status, 206L) && nzchar(cr) && grepl("^bytes[[:space:]]+0-0/", cr, ignore.case = TRUE) && one_byte

  list(
    ok = ok,
    status = status,
    message = paste0("Content-Range=", cr, "; bytes_received=", length(res$content)),
    msg = paste0("Content-Range=", cr, "; bytes_received=", length(res$content)),
    content_range = cr
  )
}
retry <- function(fun, label, attempts = 3L) {
  errs <- character()
  for (i in seq_len(attempts)) {
    if (i > 1L) Sys.sleep(c(2, 5)[min(i - 1L, 2L)])
    ans <- tryCatch(fun(), error = function(e) e)
    if (!inherits(ans, "error")) return(ans)
    errs <- c(errs, paste0("attempt ", i, ": ", conditionMessage(ans)))
  }
  stop(label, " failed: ", paste(errs, collapse = " | "), call. = FALSE)
}
parse_gt <- function(gt) {
  if (is.na(gt) || !nzchar(gt) || gt %in% c(".", "./.", ".|.")) return(c(called = 0, refonly = 0, ncall = 0, nnonref = 0))
  p <- strsplit(gt, "[|/]", perl = TRUE)[[1L]]
  p <- p[p != "." & nzchar(p)]
  if (!length(p)) return(c(called = 0, refonly = 0, ncall = 0, nnonref = 0))
  a <- suppressWarnings(as.integer(p))
  if (anyNA(a)) return(c(called = 0, refonly = 0, ncall = 0, nnonref = 0))
  c(called = 1, refonly = as.integer(all(a == 0L)), ncall = length(a), nnonref = sum(a != 0L))
}
pam_positions <- function(st, en, strand) {
  if (strand == "+") return((en - 2L):en)
  if (strand == "-") return(st:(st + 2L))
  stop("Invalid guide strand: ", strand, call. = FALSE)
}

sample_states <- function(url, sid, taxon) {
  pr <- range_probe(url)
  if (!isTRUE(pr$ok)) {
    stop(
      "VCF byte-range probe failed (HTTP ", pr$status, "; ", pr$message, "): ",
      url,
      call. = FALSE
    )
  }
  ir <- range_probe(paste0(url, ".tbi"))
  if (!isTRUE(ir$ok)) {
    stop(
      "VCF index byte-range probe failed (HTTP ", ir$status, "; ", ir$message, "): ",
      paste0(url, ".tbi"),
      call. = FALSE
    )
  }
  log_step(STEP, paste0("Range access confirmed for ", sid, ": ", pr$content_range, "; index ", ir$content_range))
  tf <- Rsamtools::TabixFile(url, index = paste0(url, ".tbi"))
  hdr <- retry(function() VariantAnnotation::scanVcfHeader(tf), paste0("header for ", sid))
  hs <- as.character(VariantAnnotation::samples(hdr))
  if (length(hs) != 1L) stop("Expected one VCF sample for ", sid, "; header has ", length(hs), call. = FALSE)
  if (!identical(hs[[1L]], sid)) stop("VCF sample/header mismatch: selected ", sid, " but header is ", hs[[1L]], call. = FALSE)
  gr <- GenomicRanges::GRanges(seqnames = windows$contig, ranges = IRanges::IRanges(windows$start, windows$end))
  param <- VariantAnnotation::ScanVcfParam(which = gr, samples = hs, fixed = NA_character_, info = NA_character_, geno = "GT")
  vcf <- retry(function() VariantAnnotation::readVcf(tf, genome = "AgamP4", param = param), paste0("targeted all-sites query for ", sid))
  if (!nrow(vcf)) stop("No records returned for locked target windows in ", sid, call. = FALSE)
  rr <- SummarizedExperiment::rowRanges(vcf)
  keys <- paste(as.character(seqnames(rr)), as.integer(start(rr)), sep = ":")
  if (anyDuplicated(keys)) stop("Duplicate genomic positions returned in all-sites VCF for ", sid, call. = FALSE)
  gt <- VariantAnnotation::geno(vcf)$GT
  if (is.null(gt) || ncol(gt) != 1L) stop("Expected one GT column in ", sid, call. = FALSE)
  gtv <- as.character(gt[, 1L])
  names(gtv) <- keys

  rows <- vector("list", nrow(targets))
  for (i in seq_len(nrow(targets))) {
    t <- targets[i]
    cc <- as.character(t$genomic_seqid); st <- as.integer(t$genomic_start); en <- as.integer(t$genomic_end)
    pos <- st:en; kk <- paste(cc, pos, sep = ":"); vals <- unname(gtv[kk])
    present <- sum(!is.na(vals))
    ppam <- pam_positions(st, en, as.character(t$guide_genomic_strand))
    called_all <- present == 23L; exact <- TRUE; pam_called <- TRUE; pam_intact <- TRUE
    called_alleles <- 0L; nonref_alleles <- 0L; missing <- integer()
    for (j in seq_along(pos)) {
      p <- pos[[j]]; g <- vals[[j]]
      if (is.na(g)) {
        called_all <- FALSE; exact <- FALSE; missing <- c(missing, p)
        if (p %in% ppam) { pam_called <- FALSE; pam_intact <- FALSE }
        next
      }
      z <- parse_gt(g)
      called_alleles <- called_alleles + as.integer(z[["ncall"]])
      nonref_alleles <- nonref_alleles + as.integer(z[["nnonref"]])
      if (z[["called"]] == 0L) {
        called_all <- FALSE; exact <- FALSE
        if (p %in% ppam) { pam_called <- FALSE; pam_intact <- FALSE }
      } else if (z[["refonly"]] == 0L) {
        exact <- FALSE
        if (p %in% ppam) pam_intact <- FALSE
      }
    }
    if (!called_all) exact <- FALSE
    if (!pam_called) pam_intact <- FALSE
    rows[[i]] <- data.table(
      sample_id = sid, taxon = taxon, population_site_id = as.character(t$population_site_id), contig = cc,
      records_present_23bp = present, callable_23bp = called_all,
      exact_23bp = if (called_all) exact else NA,
      pam_callable = pam_called, pam_intact = if (pam_called) pam_intact else NA,
      called_alleles_23bp = called_alleles, nonref_alleles_23bp = nonref_alleles,
      missing_positions = paste(missing, collapse = ";")
    )
  }
  rbindlist(rows)
}

states_list <- list(); errors <- list(); ns <- nrow(selected)
max_fail <- as.numeric(ev$max_sample_failure_fraction %||% 0.10)
for (j in seq_len(ns)) {
  sid <- as.character(selected$sample_id[[j]]); tx <- as.character(selected$taxon[[j]]); url <- as.character(selected$vcf_url[[j]])
  log_step(STEP, paste0("[", j, "/", ns, "] Ag3 targeted all-sites VCF: ", sid, " (", tx, ")"))
  if (!identical(digest::digest(file = lock_file, algo = "sha256", serialize = FALSE), lock_hash)) {
    fail_step(STEP, "External target-lock SHA-256 changed during acquisition; aborting.")
  }
  ans <- tryCatch(sample_states(url, sid, tx), error = function(e) e)
  if (inherits(ans, "error")) {
    errors[[length(errors) + 1L]] <- data.table(sample_id = sid, taxon = tx, vcf_url = url, error = conditionMessage(ans))
    log_warning(STEP, paste0("Sample failed: ", sid, ": ", conditionMessage(ans)))
    if (j >= 10L && length(errors) / j > max_fail) {
      err <- if (length(errors)) rbindlist(errors, fill = TRUE) else data.table(sample_id=character(),taxon=character(),vcf_url=character(),error=character())
      atomic_fwrite(err, file.path(out_dir, "sample_errors.tsv"), sep = "\t")
      fail_step(STEP, paste0("Sample failure fraction exceeded ", scales::percent(max_fail), "; stopping."))
    }
  } else {
    states_list[[length(states_list) + 1L]] <- ans
  }
}
states <- if (length(states_list)) rbindlist(states_list, use.names = TRUE, fill = TRUE) else data.table()
err <- if (length(errors)) rbindlist(errors, use.names = TRUE, fill = TRUE) else data.table(sample_id=character(),taxon=character(),vcf_url=character(),error=character())
atomic_fwrite(err, file.path(out_dir, "sample_errors.tsv"), sep = "\t")
if (!nrow(states)) fail_step(STEP, "No Ag3 external sample-site states were generated.")
good <- uniqueN(states$sample_id)
if (good < max(10L, floor(0.8 * ns))) fail_step(STEP, paste0("Only ", good, "/", ns, " selected samples succeeded."))
expected <- good * nrow(targets)
if (nrow(states) != expected) fail_step(STEP, paste0("External state matrix incomplete: ", nrow(states), " rows; expected ", expected, "."))

agg_one <- function(d) {
  called <- d[callable_23bp == TRUE]
  pamc <- d[pam_callable == TRUE]
  ca <- sum(d$called_alleles_23bp, na.rm = TRUE); na <- sum(d$nonref_alleles_23bp, na.rm = TRUE)
  list(
    n_samples_total = nrow(d), n_samples_called_23bp = nrow(called),
    n_samples_exact_23bp = sum(called$exact_23bp == TRUE, na.rm = TRUE),
    n_samples_pam_called = nrow(pamc), n_samples_pam_intact = sum(pamc$pam_intact == TRUE, na.rm = TRUE),
    callable_sample_fraction = if (nrow(d)) nrow(called) / nrow(d) else NA_real_,
    exact_23bp_fraction = if (nrow(called)) mean(called$exact_23bp, na.rm = TRUE) else NA_real_,
    pam_intact_fraction = if (nrow(pamc)) mean(pamc$pam_intact, na.rm = TRUE) else NA_real_,
    called_alleles_23bp = ca, nonreference_alleles_23bp = na,
    nonreference_allele_fraction_23bp = if (ca > 0) na / ca else NA_real_
  )
}
site <- states[, agg_one(.SD), by = population_site_id]
tax <- states[, agg_one(.SD), by = .(population_site_id, taxon)]
ann <- unique(targets[, .(population_site_id, gene_id, selection_reason)])
site <- merge(site, ann, by = "population_site_id", all.x = TRUE, sort = FALSE)
tax <- merge(tax, ann, by = "population_site_id", all.x = TRUE, sort = FALSE)
atomic_fwrite(states, file.path(out_dir, "sample_site_states.tsv"), sep = "\t")
atomic_fwrite(site, file.path(out_dir, "site_summary.tsv"), sep = "\t")
atomic_fwrite(tax, file.path(out_dir, "site_taxon_summary.tsv"), sep = "\t")

manifest <- list(
  schema_version = "ag3-external-windows-rsamtools-v3",
  run_mode = run_mode,
  external_dataset = "Ag3.0 public Sanger per-sample all-sites VCF",
  access_engine = "R/Bioconductor Rsamtools + VariantAnnotation remote TabixFile",
  target_lock_file = lock_file,
  target_lock_sha256 = lock_hash,
  n_locked_targets = nrow(targets), n_target_windows = nrow(windows), target_window_logical_bp = window_bp,
  n_selected_samples = ns, n_successful_samples = good, n_failed_samples = nrow(err),
  whole_vcf_files_downloaded = 0L, discovery_retuning_allowed = FALSE,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(manifest, file.path(out_dir, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
atomic_writeLines("ag3-external-windows-rsamtools-v3", file.path(out_dir, "ACQUISITION_COMPLETE.ok"))
write_session_info(file.path("logs", paste0("05c_ag3_", run_mode, "_acquisition_sessionInfo.txt")))
log_step(STEP, paste0("Ag3 ", run_mode, " acquisition completed: ", good, " samples x ", nrow(targets), " locked targets; whole VCF downloads = 0"))
