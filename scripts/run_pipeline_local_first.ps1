$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
$projectRoot = (Get-Location).Path
$venvPython = $null

# 1) Explicit override.
if ($env:MOSQEDIT_AG3_VENV) {
    $candidate = Join-Path ([System.IO.Path]::GetFullPath($env:MOSQEDIT_AG3_VENV)) "Scripts\python.exe"
    if (Test-Path $candidate) { $venvPython = $candidate }
}

# 2) Interpreter pointer written by setup_ag3_python.ps1.
if ($null -eq $venvPython) {
    $pointer = Join-Path $projectRoot "metadata\ag3_python_path.txt"
    if (Test-Path $pointer) {
        $candidate = (Get-Content $pointer -Raw).Trim()
        if ($candidate -and (Test-Path $candidate)) { $venvPython = $candidate }
    }
}

# 3) New default short-path environment.
if ($null -eq $venvPython -and $env:LOCALAPPDATA) {
    $candidate = Join-Path $env:LOCALAPPDATA "MosqEditR\ag3_venv\Scripts\python.exe"
    if (Test-Path $candidate) { $venvPython = $candidate }
}

# 4) Legacy project-local environment, retained as a compatibility fallback.
if ($null -eq $venvPython) {
    $candidate = Join-Path $projectRoot ".venv_ag3\Scripts\python.exe"
    if (Test-Path $candidate) { $venvPython = $candidate }
}

if ($null -eq $venvPython) {
    throw "Ag3 Python environment not found. Run: powershell -ExecutionPolicy Bypass -File scripts/setup_ag3_python.ps1"
}

Write-Host "Using Python: $venvPython"
Write-Host "[1/3] Running frozen-panel preflight..."
& Rscript --vanilla R\00_preflight_and_freeze.R
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[2/3] Acquiring ONLY targeted Ag3 metadata/haplotypes/accessibility..."
& $venvPython scripts\acquire_ag3_targeted.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[3/3] Running local-only R analysis from Step 01 onward..."
$env:MOSQEDIT_START_STEP = "01"
& Rscript --vanilla run_all.R
exit $LASTEXITCODE
