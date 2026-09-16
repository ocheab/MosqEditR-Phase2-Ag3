# MosqEditR: Phase-2 Discovery and Ag3 External Validation

## Overview

This repository contains the reproducible analysis workflow used to identify, rank, freeze, and externally validate CRISPR target sites in *Anopheles gambiae* and *Anopheles coluzzii*.

The workflow is intentionally divided into two analytically independent stages:

1. **Discovery and target freezing using Ag1000G Phase 2 AR1 phased haplotypes**
2. **Independent external validation using Ag1000G Ag3 genotype data**

Ag3 data are used only for external validation. They do **not** alter Phase-2 target selection, discovery scores, ranks, classes, thresholds, or optimized multiplex portfolios.

## Frozen discovery panel

The original frozen CRISPR panel contains:

- 106 genes
- 530 candidate sites
- exactly 5 candidate sites per gene
- 515 sites on canonical *Anopheles gambiae* chromosome arms
- 15 unresolved/unplaced sites

The external-validation lock contains:

- 354 sites
- 102 genes
- SHA-256 target lock:
  `eb8dd496c21690750a6a6d1fe6c4f19acbdfb138e2a31de3990687d6e942a2f3`

The target lock must not be regenerated or modified after external validation begins.

## Primary Phase-2 cohort

The primary discovery cohort contains 938 samples with explicit Phase-2 release-population taxon assignments.

Population codes ending in:

- `gam` are treated as *An. gambiae*
- `col` are treated as *An. coluzzii*

Mixed or undetermined groups are excluded from the primary species-assigned discovery analysis.

Contig-specific denominators are preserved. Missing sample/contig combinations, including the reduced X-chromosome sample axis, are treated as missing and are never interpreted as reference.

## Independent Ag3 validation cohort

The final external-validation cohort contains:

- 500 mosquitoes
- 250 *An. gambiae*
- 250 *An. coluzzii*

The validation sampling frame is the Ag3.0 West African collection defined from published Ag3 metadata. Sampling is deterministic, country-stratified, and SHA-256 based.

All 354 locked sites passed the prespecified external sample-coverage requirement.

## Sequence-aware validation metrics

The Ag3 validator computes:

- `strict_genotype_exact_23bp_fraction`
- `protospacer_exact_20bp_fraction`
- `pam_exact_3bp_fraction`
- `pam_ngg_intact_fraction`
- `functional_target_intact_fraction`
- `max_position_nonreference_allele_fraction_23bp`

Important interpretation:

- Phase-2 exactness is a **phased-haplotype** quantity.
- Ag3 strict exactness is an **unphased-genotype** quantity.

These metrics are therefore compared primarily using rank-concordance analyses rather than treated as numerically interchangeable measurements.

For SpCas9, PAM functionality is evaluated as **NGG**. Variation at the PAM N position may preserve a functional PAM even when the exact three-base PAM sequence changes.

## Final external-validation results

Primary 354-site analysis, averaging the two Ag3 taxa per locked site:

| Comparison | Spearman rho | Gene-cluster bootstrap 95% CI |
|---|---:|---:|
| Phase-2 global exact vs Ag3 strict 23-bp exact | 0.823 | 0.762–0.872 |
| Phase-2 population-q10 exact vs Ag3 strict exact | 0.779 | 0.710–0.839 |
| Phase-2 PAM intact vs Ag3 exact 3-bp PAM | 0.613 | 0.526–0.688 |
| Phase-2 PAM intact vs Ag3 functional NGG PAM | 0.509 | 0.389–0.614 |
| Phase-2 max alternate AF vs Ag3 max-position nonreference AF | 0.748 | 0.668–0.817 |

Allele-frequency quantitative agreement:

- Pearson r: 0.706
- Lin's CCC: 0.695
- mean Ag3-minus-Phase-2 difference: -0.043
- MAE: 0.106
- RMSE: 0.191
- calibration slope: 0.673
- calibration intercept: 0.010

Ag3 validation changed:

- discovery ranks: **0**
- discovery thresholds: **0**

## Repository layout

```text
.
├── R/
│   ├── 00_preflight_and_freeze.R
│   ├── 01_prepare_phase2_metadata.R
│   ├── 02_extract_phase2_target_variation.R
│   ├── 03_compute_population_site_metrics.R
│   ├── 04_optimize_multiplex_portfolios.R
│   ├── 05_gene_robustness_and_reranking.R
│   ├── 05b_freeze_external_validation_targets.R
│   ├── 06_sensitivity_and_validation.R
│   ├── 06a_ag3_external_validation.R
│   ├── 06b_build_ag3_phase2_sequence_aware_comparison.R
│   ├── 07_ag3_manuscript_grade_validation_FAST.R
│   └── 07g_ag3_gene_cluster_bootstrap_sensitivity.R
├── scripts/
│   ├── acquire_ag3_external_validation.py
│   └── acquire_ag3_external_validation.ps1
├── metadata/
├── data_raw/                # not committed
├── data_processed/          # selected compact derived tables only
├── outputs/                 # selected final tables/figures only
├── docs/
├── environment/
├── README.md
├── CITATION.cff
└── .gitignore
```

## Recommended execution order

### Phase-2 discovery

```bash
Rscript R/00_preflight_and_freeze.R
Rscript R/01_prepare_phase2_metadata.R
Rscript R/02_extract_phase2_target_variation.R
Rscript R/03_compute_population_site_metrics.R
Rscript R/04_optimize_multiplex_portfolios.R
Rscript R/05_gene_robustness_and_reranking.R
Rscript R/05b_freeze_external_validation_targets.R
Rscript R/06_sensitivity_and_validation.R
```

### Ag3 acquisition

On Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass `
  -File ".\scripts\acquire_ag3_external_validation.ps1" `
  --final `
  --workers 2
```

The acquisition is resumable. Completed sample checkpoints are validated against:

- target-lock SHA-256
- sample ID and taxon
- exact locked-site set
- complete 23-position VCF coverage
- sequence-aware checkpoint schema

Interrupted runs should be resumed by rerunning the same command. Do not delete the final checkpoint directory.

### Ag3 external validation

```bash
Rscript R/06a_ag3_external_validation.R
Rscript R/06b_build_ag3_phase2_sequence_aware_comparison.R
Rscript R/07_ag3_manuscript_grade_validation_FAST.R
Rscript R/07g_ag3_gene_cluster_bootstrap_sensitivity.R
```

## Reproducibility safeguards

The workflow enforces the following rules:

1. Ag3 data cannot modify Phase-2 discovery.
2. Missing external VCF records remain uncallable and are never treated as reference.
3. The external-validation target lock is checked using SHA-256.
4. External validation requires the complete locked sample cohort.
5. Multiple CRISPR sites within genes are handled using gene-cluster bootstrap sensitivity analyses.
6. Phase-2 exactness and Ag3 strict exactness are not treated as numerically identical estimands.
7. PAM functionality is sequence aware and respects the NGG rule.

See `docs/REPRODUCIBILITY.md` for additional details.

## Data availability

Large public genomic resources and raw external VCF data should not be committed to GitHub.

The repository should contain:

- analysis code
- metadata needed to reconstruct the cohort
- target-lock files
- compact processed tables required to reproduce figures/statistics where redistribution is permitted
- checksums
- documentation

Raw public data should instead be reacquired from their authoritative public repositories using the included acquisition scripts.

## Citation

xxxx

## License

MIT. Dataset licenses and upstream Ag1000G/Ag3 terms remain separate from the software license.
