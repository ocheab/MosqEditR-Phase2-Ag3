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
& $venvPython scripts\acquire_ag3_targeted.py @args
exit $LASTEXITCODE
