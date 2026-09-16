# Ag3 External Validation Methods

## Cohort

The final validation cohort contains 500 Ag3 mosquitoes:

- 250 *Anopheles gambiae*
- 250 *Anopheles coluzzii*

Samples are drawn from an Ag3.0 West African sampling frame using deterministic country-stratified SHA-256 selection.

## Target set

The validation target lock contains 354 sites from 102 genes.

Each target is represented as an exact 23-bp CRISPR target sequence.

## Acquisition engine

The validation acquisition engine uses:

- public per-sample compressed VCF files
- Tabix indexes
- HTTP byte-range requests
- BGZF block decoding
- exact interval extraction

Acquisition is checkpointed per sample and is safe to resume after transient network failure.

## Missingness

Missing VCF records are **not** interpreted as reference genotypes.

A target is considered fully covered only when all required positions are represented according to the validator's callability rules.

## Sequence-aware metrics

### Strict genotype exactness

`strict_genotype_exact_23bp_fraction`

Fraction of callable samples whose diploid genotype is sequence-consistent with the full frozen 23-bp target.

### Protospacer exactness

`protospacer_exact_20bp_fraction`

Fraction of callable samples matching the frozen 20-bp protospacer.

### Exact PAM

`pam_exact_3bp_fraction`

Fraction matching the exact frozen three-base PAM.

### Functional NGG PAM

`pam_ngg_intact_fraction`

Fraction retaining an SpCas9-compatible NGG PAM. Variation at the first PAM base can remain functionally compatible.

### Functional target integrity

`functional_target_intact_fraction`

Sequence-aware target integrity combining protospacer conservation with functional PAM status.

### Maximum-position nonreference allele frequency

`max_position_nonreference_allele_fraction_23bp`

Maximum observed nonreference/mismatch allele frequency over the 23 target positions.

## Comparison with Phase 2

Phase-2 exactness is derived from phased haplotypes; Ag3 strict exactness is based on unphased genotypes.

Therefore exactness comparisons are interpreted primarily using:

- Spearman rank correlation
- Pearson correlation as a descriptive complement
- site-level mean across taxa
- taxon-specific analyses
- locked-site bootstrap confidence intervals
- gene-cluster bootstrap sensitivity confidence intervals

For allele frequencies, additional quantitative agreement measures include:

- bias
- MAE
- RMSE
- Lin's CCC
- calibration regression

## Bootstrap design

### Locked-site bootstrap

Resampling unit: `population_site_id`

Both taxon rows for a selected site are retained together.

Replicates: 5,000.

### Gene-cluster bootstrap

Resampling unit: `gene_id`

All locked sites within a selected gene and all relevant taxon rows are retained together.

Replicates: 5,000.

The gene-cluster analysis is a sensitivity analysis for within-gene dependence.
