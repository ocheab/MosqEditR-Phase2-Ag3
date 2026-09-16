# Workflow

## Analytical boundary

The central design principle is a strict boundary between **discovery** and **external validation**.

### Discovery phase

Ag1000G Phase 2 AR1 phased haplotypes are used to:

- quantify target-sequence variation
- compute population-level robustness metrics
- classify candidate sites
- rank genes
- optimize multiplex portfolios
- freeze external-validation targets

### Validation phase

Ag3 genotype data are used to:

- evaluate locked sites
- quantify sequence conservation
- characterize PAM and protospacer integrity
- assess cross-release allele-frequency concordance
- identify taxon heterogeneity
- quantify robustness using bootstrap sensitivity analyses

Ag3 is not used to:

- add or remove locked sites
- alter discovery weights
- alter thresholds
- change discovery classes
- rerank Phase-2 results
- re-optimize multiplexes

## Script map

| Script | Purpose |
|---|---|
| `00_preflight_and_freeze.R` | preflight checks and frozen configuration |
| `01_prepare_phase2_metadata.R` | construct primary species-assigned Phase-2 metadata |
| `02_extract_phase2_target_variation.R` | target-level Phase-2 variation extraction |
| `03_compute_population_site_metrics.R` | population/site robustness metrics |
| `04_optimize_multiplex_portfolios.R` | pair/triple multiplex optimization |
| `05_gene_robustness_and_reranking.R` | gene-level robustness and ranking |
| `05b_freeze_external_validation_targets.R` | freeze independent-validation targets and hash |
| `06_sensitivity_and_validation.R` | discovery sensitivity checks |
| `acquire_ag3_external_validation.py` | public Ag3 per-sample VCF target extraction |
| `06a_ag3_external_validation.R` | prespecified external-validation evaluation |
| `06b_build_ag3_phase2_sequence_aware_comparison.R` | build taxon-level Phase-2/Ag3 comparison |
| `07_ag3_manuscript_grade_validation_FAST.R` | final manuscript-grade validation statistics and figures |
| `07g_ag3_gene_cluster_bootstrap_sensitivity.R` | gene-cluster bootstrap sensitivity analysis |

## Frozen target hash

The external target lock must resolve to:

```text
eb8dd496c21690750a6a6d1fe6c4f19acbdfb138e2a31de3990687d6e942a2f3
```

Any mismatch should stop the validation workflow.
