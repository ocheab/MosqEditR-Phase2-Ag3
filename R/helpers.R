# MosqEditR Manuscript 2 -- shared robust helper functions
# This file contains no analysis-specific side effects.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

locate_project_root <- function(start = getwd()) {
  candidates <- character()
  add_candidate <- function(x) {
    if (!is.null(x) && length(x) && nzchar(x)) {
      x <- tryCatch(normalizePath(x, winslash = "/", mustWork = FALSE), error = function(e) x)
      candidates <<- unique(c(candidates, x))
    }
  }

  add_candidate(start)

  # When sourced, sys.frames() often contains the current script path in $ofile.
  for (i in rev(seq_along(sys.frames()))) {
    of <- sys.frames()[[i]]$ofile %||% NULL
    if (!is.null(of) && nzchar(of)) {
      of <- tryCatch(normalizePath(of, winslash = "/", mustWork = FALSE), error = function(e) of)
      add_candidate(dirname(of))
      add_candidate(dirname(dirname(of)))
    }
  }

  # Walk upward from every candidate.
  expanded <- character()
  for (x in candidates) {
    cur <- x
    for (j in 0:8) {
      expanded <- unique(c(expanded, cur))
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }

  for (x in expanded) {
    if (file.exists(file.path(x, "analysis_config.yml")) &&
        dir.exists(file.path(x, "R")) &&
        dir.exists(file.path(x, "data_input"))) {
      return(x)
    }
  }

  stop(
    "Could not locate the MosqEditR project root. Run the script from the project folder, ",
    "or set the working directory to the folder containing analysis_config.yml, R/, and data_input/.",
    call. = FALSE
  )
}

set_project_root <- function() {
  root <- locate_project_root()
  if (!identical(normalizePath(getwd(), winslash = "/", mustWork = FALSE), root)) {
    setwd(root)
  }
  ensure_directories()
  invisible(root)
}

ensure_directories <- function() {
  dirs <- c(
    "logs", "data_raw", "data_processed", "figures", "tables",
    "supplement", "release"
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

require_packages <- function(pkgs, bioc = character()) {
  all_pkgs <- unique(c(pkgs, bioc))
  missing <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". Run source(\"R/install_dependencies.R\") first.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

read_config <- function(path = "analysis_config.yml") {
  assert_file(path)
  require_packages("yaml")
  cfg <- yaml::read_yaml(path)
  if (!is.list(cfg)) stop("analysis_config.yml did not parse to a YAML mapping.", call. = FALSE)
  cfg
}

cfg_get <- function(cfg, ..., default = NULL, required = FALSE) {
  keys <- list(...)
  x <- cfg
  for (k in keys) {
    if (!is.list(x) || is.null(x[[k]])) {
      if (required) {
        stop("Missing required configuration key: ", paste(unlist(keys), collapse = " -> "), call. = FALSE)
      }
      return(default)
    }
    x <- x[[k]]
  }
  x
}

log_step <- function(step, message, level = "INFO") {
  ensure_directories()
  stamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  line <- paste(stamp, level, step, message, sep = "\t")
  cat(line, "\n")
  cat(line, "\n", file = "logs/pipeline.log", append = TRUE)
  invisible(line)
}

log_warning <- function(step, message) {
  log_step(step, message, level = "WARNING")
  warning(message, call. = FALSE, immediate. = TRUE)
}

fail_step <- function(step, message) {
  log_step(step, message, level = "ERROR")
  stop(message, call. = FALSE)
}

assert_file <- function(path, nonempty = TRUE) {
  if (!file.exists(path)) stop("Required file is missing: ", path, call. = FALSE)
  if (nonempty && isTRUE(file.info(path)$size <= 0)) stop("Required file is empty: ", path, call. = FALSE)
  invisible(path)
}

assert_files <- function(paths, nonempty = TRUE) {
  invisible(lapply(paths, assert_file, nonempty = nonempty))
}

assert_columns <- function(x, cols, object_name = deparse(substitute(x))) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) {
    stop(object_name, " is missing required column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_unique <- function(x, col, object_name = deparse(substitute(x)), allow_na = FALSE) {
  assert_columns(x, col, object_name)
  v <- x[[col]]
  if (!allow_na && anyNA(v)) stop(object_name, "$", col, " contains missing values.", call. = FALSE)
  dup <- duplicated(v) & !is.na(v)
  if (any(dup)) {
    vals <- unique(v[dup])
    stop(object_name, "$", col, " is not unique. Example duplicate(s): ",
         paste(utils::head(vals, 10L), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_no_missing <- function(x, cols, object_name = deparse(substitute(x))) {
  assert_columns(x, cols, object_name)
  bad <- vapply(cols, function(z) any(is.na(x[[z]]) | (is.character(x[[z]]) & !nzchar(trimws(x[[z]])))), logical(1))
  if (any(bad)) {
    stop(object_name, " contains missing/blank values in required column(s): ",
         paste(cols[bad], collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_numeric_range <- function(x, col, lower = -Inf, upper = Inf, object_name = deparse(substitute(x))) {
  assert_columns(x, col, object_name)
  v <- suppressWarnings(as.numeric(x[[col]]))
  bad <- !is.na(v) & (v < lower | v > upper)
  if (any(bad)) {
    stop(object_name, "$", col, " contains values outside [", lower, ", ", upper, "].", call. = FALSE)
  }
  invisible(TRUE)
}

normalize_flag <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

normalize_taxon <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("^anopheles\\s+", "", y)
  y <- gsub("^a\\.\\s*", "", y)
  y <- gsub("^an\\.\\s*", "", y)
  y <- gsub("[_-]", " ", y)
  y <- trimws(y)
  out <- y
  out[grepl("(^|\\s)gambiae($|\\s)", y)] <- "gambiae"
  out[grepl("(^|\\s)coluzzii($|\\s)", y)] <- "coluzzii"
  out[is.na(x) | !nzchar(trimws(as.character(x)))] <- NA_character_
  out
}

resolve_column <- function(x, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) return(hit[[1]])
  if (required) {
    stop("Could not identify ", label, ". Tried column names: ", paste(candidates, collapse = ", "),
         ". Available columns: ", paste(names(x), collapse = ", "), call. = FALSE)
  }
  NA_character_
}

safe_mean <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  mean(x)
}

safe_sd <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(default)
  stats::sd(x)
}

safe_min <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  min(x)
}

safe_max <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  max(x)
}

safe_median <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  stats::median(x)
}

safe_quantile <- function(x, p, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 8, na.rm = TRUE))
}

safe_cor <- function(x, y, method = "spearman", default = NA_real_) {
  x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L) return(default)
  if (length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(default)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

renorm_weighted <- function(values, weights) {
  values <- suppressWarnings(as.numeric(values))
  weights <- suppressWarnings(as.numeric(weights))
  ok <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(ok)) return(NA_real_)
  sum(values[ok] * weights[ok]) / sum(weights[ok])
}

validate_weights <- function(weights, label = "weights", tol = 1e-8) {
  x <- unlist(weights, use.names = TRUE)
  if (!length(x) || any(!is.finite(x)) || any(x < 0)) stop(label, " must contain finite non-negative numbers.", call. = FALSE)
  if (abs(sum(x) - 1) > tol) stop(label, " must sum to 1. Current sum = ", signif(sum(x), 8), call. = FALSE)
  invisible(x)
}

parse_gt_haplotypes <- function(gt, require_phased = FALSE) {
  gt <- as.character(gt)
  n <- length(gt)
  a1 <- rep(NA_character_, n)
  a2 <- rep(NA_character_, n)
  is_missing <- is.na(gt) | gt %in% c("", ".", "./.", ".|.")
  is_unphased_nonmissing <- !is_missing & grepl("/", gt, fixed = TRUE)
  if (require_phased && any(is_unphased_nonmissing)) {
    ex <- unique(gt[is_unphased_nonmissing])
    stop("Encountered non-missing unphased genotype(s) in a phased analysis: ",
         paste(utils::head(ex, 10L), collapse = ", "), call. = FALSE)
  }

  idx <- which(!is_missing)
  for (i in idx) {
    parts <- strsplit(gt[[i]], "[|/]", perl = TRUE)[[1]]
    if (length(parts) >= 1L && parts[[1]] != "." && nzchar(parts[[1]])) a1[[i]] <- parts[[1]]
    if (length(parts) >= 2L && parts[[2]] != "." && nzchar(parts[[2]])) a2[[i]] <- parts[[2]]
    # Haploid calls deliberately keep a2 = NA rather than duplicating a1.
  }
  list(a1 = a1, a2 = a2, unphased_nonmissing = is_unphased_nonmissing)
}

allele_is_reference <- function(a) {
  a <- as.character(a)
  out <- rep(NA, length(a))
  ok <- !is.na(a) & nzchar(a)
  out[ok] <- a[ok] == "0"
  out
}

interleave_haplotypes <- function(a1, a2) {
  if (length(a1) != length(a2)) stop("Haplotype vectors must have equal length.", call. = FALSE)
  as.vector(rbind(a1, a2))
}

make_haplotype_names <- function(sample_ids) {
  interleave_haplotypes(paste0(sample_ids, "|h1"), paste0(sample_ids, "|h2"))
}

haplotype_to_sample <- function(hap_names) {
  sub("\\|h[12]$", "", hap_names)
}

combine_exact_state <- function(current, new_state) {
  if (length(current) != length(new_state)) stop("State vectors have different lengths.", call. = FALSE)
  out <- as.logical(current)
  new_state <- as.logical(new_state)
  out[!is.na(new_state) & !new_state] <- FALSE
  out[is.na(new_state) & (is.na(out) | out)] <- NA
  out
}

at_least_one_true <- function(x) {
  x <- as.logical(x)
  if (any(x %in% TRUE, na.rm = TRUE)) return(TRUE)
  if (anyNA(x)) return(NA)
  FALSE
}

all_true_state <- function(x) {
  x <- as.logical(x)
  if (any(x %in% FALSE, na.rm = TRUE)) return(FALSE)
  if (anyNA(x)) return(NA)
  TRUE
}

atomic_fwrite <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  data.table::fwrite(x, tmp, ...)
  if (!file.rename(tmp, path)) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    if (!ok) stop("Failed to write output: ", path, call. = FALSE)
    unlink(tmp)
  }
  invisible(path)
}

atomic_saveRDS <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (!file.rename(tmp, path)) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    if (!ok) stop("Failed to write output: ", path, call. = FALSE)
    unlink(tmp)
  }
  invisible(path)
}

atomic_writeLines <- function(text, path, useBytes = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(text, tmp, useBytes = useBytes)
  if (!file.rename(tmp, path)) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    if (!ok) stop("Failed to write output: ", path, call. = FALSE)
    unlink(tmp)
  }
  invisible(path)
}

write_checksum <- function(files, outfile) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) return(invisible(NULL))
  tab <- data.table::data.table(file = files, md5 = unname(tools::md5sum(files)))
  atomic_fwrite(tab, outfile, sep = "\t")
  invisible(tab)
}


sha256_file <- function(path) {
  assert_file(path)
  require_packages("digest")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_sha256 <- function(files, outfile) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) return(invisible(NULL))
  require_packages(c("data.table", "digest"))
  tab <- data.table::data.table(
    file = files,
    sha256 = vapply(files, sha256_file, character(1))
  )
  atomic_fwrite(tab, outfile, sep = "\t")
  invisible(tab)
}

write_session_info <- function(path = "logs/sessionInfo.txt") {
  txt <- capture.output(utils::sessionInfo())
  atomic_writeLines(txt, path)
  invisible(path)
}

empty_table <- function(schema) {
  # schema is a named list of zero-length typed vectors.
  data.table::as.data.table(schema)
}

nonempty_string <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}
