$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
$projectRoot = (Get-Location).Path

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
if ($null -eq $venvPython) {
    throw "Hybrid Python environment not found. Run scripts\setup_hybrid_python.ps1 first."
}

Write-Host "Using Python sequence-aware pure-Tabix/BGZF engine: $venvPython"
& $venvPython scripts\acquire_ag3_external_validation.py @args
exit $LASTEXITCODE
