$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
$projectRoot = (Get-Location).Path

# Resolve the short-path helper Python environment (metadata + HDF5 only).
$venvPython = $null
if ($env:MOSQEDIT_AG3_VENV) {
    $candidate = Join-Path ([System.IO.Path]::GetFullPath($env:MOSQEDIT_AG3_VENV)) "Scripts\python.exe"
    if (Test-Path $candidate) { $venvPython = $candidate }
}
if ($null -eq $venvPython) {
    $pointer = Join-Path $projectRoot "metadata\ag3_python_path.txt"
    if (Test-Path $pointer) {
        $candidate = (Get-Content $pointer -Raw).Trim()
        if ($candidate -and (Test-Path $candidate)) { $venvPython = $candidate }
    }
}
if ($null -eq $venvPython -and $env:LOCALAPPDATA) {
    $candidate = Join-Path $env:LOCALAPPDATA "MosqEditR\ag3_venv\Scripts\python.exe"
    if (Test-Path $candidate) { $venvPython = $candidate }
}
if ($null -eq $venvPython) { throw "Helper Python environment not found. Run scripts/setup_hybrid_python.ps1 first." }

Write-Host "Using helper Python: $venvPython"
& $venvPython scripts\acquire_phase2_discovery.py @args
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

# Phase-2 phased data are acquired directly from the anonymous legacy public
# ag1000g-release GCS Zarr hierarchy. No R/Rsamtools Phase-2 VCF stage follows.
if ($LASTEXITCODE -eq 0) {
    if ($args -contains "--dry-run") { Write-Host "Dry run complete. No network access performed." }
    elseif ($args -contains "--probe-only") { Write-Host "Anonymous Phase-2 Zarr probe complete." }
    else { Write-Host "Phase-2 targeted Zarr acquisition complete." }
}
exit $LASTEXITCODE
