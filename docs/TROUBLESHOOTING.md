# Troubleshooting

## Ag3 public VCF request stalls

The acquisition script is resumable.

Do not delete the final checkpoint directory.

Rerun the same command:

```powershell
powershell -ExecutionPolicy Bypass `
  -File ".\scripts\acquire_ag3_external_validation.ps1" `
  --final `
  --workers 2
```

Valid checkpoints are reused.

## `phase2_max_alt_af.x` / `.y`

Cause: duplicate Phase-2 aliases entering a merge.

The corrected Step 06a removes pre-existing convenience aliases before attaching canonical discovery metrics.

## `log_info` not found

Use the repository's existing Step-06a logger (`log_step`) rather than introducing an undefined logging helper.

## Step 07 appears frozen at bootstrap

Do not use repeated `data.table` filtering inside each bootstrap replicate.

Use the optimized indexed implementation in:

```text
R/07_ag3_manuscript_grade_validation_FAST.R
```

## ggplot2 built under a newer patch version of R

A warning such as:

```text
package 'ggplot2' was built under R version 4.5.3
```

is normally informational when running R 4.5.2 and is not by itself evidence of analytical failure.

## Missing external records

Never replace missing Ag3 records with reference genotype state.
