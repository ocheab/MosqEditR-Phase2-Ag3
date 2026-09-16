
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20260916)

infile <- "ag3_phase2_external_validation_comparison_SEQUENCE_AWARE.csv"
outdir <- "outputs/07_ag3_manuscript_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) stop("Missing input: ", infile, call. = FALSE)
d <- fread(infile)

if (nrow(d) != 708L) stop("Expected 708 taxon-level rows.", call. = FALSE)
if (uniqueN(d$population_site_id) != 354L) stop("Expected 354 locked sites.", call. = FALSE)
if (uniqueN(d$gene_id) != 102L) stop("Expected 102 genes.", call. = FALSE)

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

# Site-level table: one row per locked site, Ag3 averaged across taxa.
site <- d[, .(
  gene_id = first(gene_id),
  phase2_global_exact = first(phase2_global_exact),
  phase2_population_q10_exact = first(phase2_population_q10_exact),
  phase2_pam_intact = first(phase2_pam_intact),
  phase2_max_alt_af = first(phase2_max_alt_af),
  ag3_strict_exact_mean = mean(strict_genotype_exact_23bp_fraction, na.rm = TRUE),
  ag3_pam_exact_mean = mean(pam_exact_3bp_fraction, na.rm = TRUE),
  ag3_pam_ngg_mean = mean(pam_ngg_intact_fraction, na.rm = TRUE),
  ag3_maxpos_nonref_af_mean =
    mean(max_position_nonreference_allele_fraction_23bp, na.rm = TRUE)
), by = population_site_id]

defs_site <- list(
  list(
    analysis = "Phase2 global exact vs mean Ag3 strict exact Spearman",
    fun = function(z) safe_spearman(
      z$phase2_global_exact,
      z$ag3_strict_exact_mean
    )
  ),
  list(
    analysis = "Phase2 population q10 exact vs mean Ag3 strict exact Spearman",
    fun = function(z) safe_spearman(
      z$phase2_population_q10_exact,
      z$ag3_strict_exact_mean
    )
  ),
  list(
    analysis = "Phase2 PAM intact vs mean Ag3 exact 3bp PAM Spearman",
    fun = function(z) safe_spearman(
      z$phase2_pam_intact,
      z$ag3_pam_exact_mean
    )
  ),
  list(
    analysis = "Phase2 PAM intact vs mean Ag3 functional NGG PAM Spearman",
    fun = function(z) safe_spearman(
      z$phase2_pam_intact,
      z$ag3_pam_ngg_mean
    )
  ),
  list(
    analysis = "Phase2 max alt AF vs mean Ag3 max-position AF Spearman",
    fun = function(z) safe_spearman(
      z$phase2_max_alt_af,
      z$ag3_maxpos_nonref_af_mean
    )
  )
)

defs_pooled <- list(
  list(
    analysis = "Phase2 global exact vs Ag3 strict exact Spearman",
    fun = function(z) safe_spearman(
      z$phase2_global_exact,
      z$strict_genotype_exact_23bp_fraction
    )
  ),
  list(
    analysis = "Phase2 population q10 exact vs Ag3 strict exact Spearman",
    fun = function(z) safe_spearman(
      z$phase2_population_q10_exact,
      z$strict_genotype_exact_23bp_fraction
    )
  ),
  list(
    analysis = "Phase2 PAM intact vs Ag3 exact 3bp PAM Spearman",
    fun = function(z) safe_spearman(
      z$phase2_pam_intact,
      z$pam_exact_3bp_fraction
    )
  ),
  list(
    analysis = "Phase2 PAM intact vs Ag3 functional NGG PAM Spearman",
    fun = function(z) safe_spearman(
      z$phase2_pam_intact,
      z$pam_ngg_intact_fraction
    )
  ),
  list(
    analysis = "Phase2 max alt AF vs Ag3 max-position AF Spearman",
    fun = function(z) safe_spearman(
      z$phase2_max_alt_af,
      z$max_position_nonreference_allele_fraction_23bp
    )
  )
)

boot_gene <- function(dt, stat_fun, B = 5000L, seed = 1L) {
  set.seed(seed)

  gene_rows <- split(seq_len(nrow(dt)), dt$gene_id, drop = TRUE)
  genes <- names(gene_rows)
  n_genes <- length(genes)

  if (n_genes != 102L) {
    stop("Gene-cluster bootstrap expected 102 genes; found ", n_genes, call. = FALSE)
  }

  vals <- numeric(B)

  for (b in seq_len(B)) {
    sampled <- sample.int(n_genes, n_genes, replace = TRUE)
    idx <- unlist(gene_rows[sampled], use.names = FALSE)
    vals[b] <- stat_fun(dt[idx])

    if (b %% 1000L == 0L) {
      cat(
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        "\tINFO\t07g\tBootstrap progress: ",
        b, "/", B, "\n", sep = ""
      )
    }
  }

  c(
    estimate = stat_fun(dt),
    ci_low = unname(quantile(vals, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(vals, 0.975, na.rm = TRUE))
  )
}

run_defs <- function(dt, defs, scope, seed_base) {
  out <- lapply(seq_along(defs), function(i) {
    cat(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      "\tINFO\t07g\t", scope, " ", i, "/", length(defs),
      ": ", defs[[i]]$analysis, "\n", sep = ""
    )

    ci <- boot_gene(
      dt,
      defs[[i]]$fun,
      B = 5000L,
      seed = seed_base + i
    )

    data.table(
      scope = scope,
      analysis = defs[[i]]$analysis,
      estimate = unname(ci["estimate"]),
      ci_low = unname(ci["ci_low"]),
      ci_high = unname(ci["ci_high"]),
      bootstrap_unit = "gene_id; all locked sites and relevant taxon rows retained",
      n_genes = uniqueN(dt$gene_id),
      B = 5000L
    )
  })
  rbindlist(out)
}

res <- rbindlist(list(
  run_defs(site, defs_site, "primary_site_mean", 2026091600L),
  run_defs(d, defs_pooled, "pooled_taxon_rows", 2026091700L)
))

outfile <- file.path(outdir, "Table_S10_gene_cluster_bootstrap_CI.csv")
fwrite(res, outfile)

cat("\nGene-cluster sensitivity analysis completed:\n")
print(res)
cat("\nWritten: ", outfile, "\n", sep = "")
