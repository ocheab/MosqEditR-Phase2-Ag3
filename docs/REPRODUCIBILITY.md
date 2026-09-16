# Reproducibility and Integrity

## Immutable discovery

External-validation data must never be used to alter the Phase-2 discovery process.

The following objects are frozen before Ag3 evaluation:

- target set
- target sequences
- discovery weights
- thresholds
- robustness classes
- site rankings
- gene rankings
- optimized pair/triple portfolios

## Hash validation

The frozen validation target set is protected using SHA-256.

Expected hash:

```text
eb8dd496c21690750a6a6d1fe6c4f19acbdfb138e2a31de3990687d6e942a2f3
```

## Checkpoint safety

A per-sample Ag3 checkpoint is reusable only when all of the following match:

- schema version
- target-lock hash
- sample ID
- taxon
- complete locked-site set
- complete positional state representation

A failed or interrupted sample is retried; completed samples are not reacquired.

## Final cohort requirement

The final external-validation run requires all selected samples to complete successfully.

The pipeline must not write a final completion marker for a partial final cohort.

## Missing data

Unobserved external VCF positions are not imputed as reference.

Missingness remains missing/uncallable.

## Numerical reproducibility

Where bootstrap procedures are used, deterministic seeds are specified in the analysis scripts.

## Recommended archival practice

For each tagged release archive:

- source scripts
- metadata
- frozen target lock
- SHA-256 checksums
- compact processed tables
- final figures
- `sessionInfo()`
- Python version and package versions

Do not archive multi-gigabyte public raw VCF files in Git.
