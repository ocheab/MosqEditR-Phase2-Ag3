$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
$projectRoot = (Get-Location).Path

# IMPORTANT (Windows): keep the virtual environment outside the project tree.
# The project path may be long, and some Python packages (e.g. jedi/typeshed)
# contain deeply nested files that can exceed the legacy Windows MAX_PATH limit.
# Override with MOSQEDIT_AG3_VENV if desired.
if ($env:MOSQEDIT_AG3_VENV) {
    $venvDir = [System.IO.Path]::GetFullPath($env:MOSQEDIT_AG3_VENV)
} else {
    if (-not $env:LOCALAPPDATA) {
        throw "LOCALAPPDATA is unavailable. Set MOSQEDIT_AG3_VENV to a short writable path, e.g. C:\\ag3env."
    }
    $venvDir = Join-Path $env:LOCALAPPDATA "MosqEditR\ag3_venv"
}

$launcher = $null
$versionArg = $null
$pyCmd = Get-Command py -ErrorAction SilentlyContinue
if ($null -ne $pyCmd) {
    foreach ($v in @("3.12", "3.11", "3.10")) {
        $code = 1
        try {
            & py "-$v" -c "import sys; raise SystemExit(0 if (3,10) <= sys.version_info[:2] < (3,13) else 1)" 2>$null
            $code = $LASTEXITCODE
        } catch {
            $code = 1
        }
        if ($code -eq 0) {
            $launcher = "py"
            $versionArg = "-$v"
            break
        }
    }
}

if ($null -eq $launcher) {
    $pythonCmd = Get-Command python -ErrorAction SilentlyContinue
    if ($null -ne $pythonCmd) {
        $code = 1
        try {
            & python -c "import sys; raise SystemExit(0 if (3,10) <= sys.version_info[:2] < (3,13) else 1)"
            $code = $LASTEXITCODE
        } catch {
            $code = 1
        }
        if ($code -eq 0) {
            $launcher = "python"
            $versionArg = $null
        }
    }
}

if ($null -eq $launcher) {
    throw "Python 3.10-3.12 is required (3.12 recommended). Install Python 3.12, then rerun this script."
}

$venvParent = Split-Path -Parent $venvDir
New-Item -ItemType Directory -Force -Path $venvParent | Out-Null

# A failed previous install can leave a half-populated environment. If the
# interpreter is missing, remove the directory before recreating it.
$venvPython = Join-Path $venvDir "Scripts\python.exe"
if ((Test-Path $venvDir) -and -not (Test-Path $venvPython)) {
    Write-Host "Removing incomplete Python environment: $venvDir"
    Remove-Item -Recurse -Force $venvDir
}

if (-not (Test-Path $venvPython)) {
    Write-Host "Creating short-path Python environment: $venvDir"
    if ($launcher -eq "py") {
        & py $versionArg -m venv $venvDir
    } else {
        & python -m venv $venvDir
    }
    if ($LASTEXITCODE -ne 0) { throw "Failed to create Python environment at $venvDir." }
}

$venvPython = Join-Path $venvDir "Scripts\python.exe"
if (-not (Test-Path $venvPython)) { throw "Virtual-environment Python was not created: $venvPython" }

Write-Host "Using Python environment: $venvDir"
& $venvPython -m pip install --upgrade pip setuptools wheel
if ($LASTEXITCODE -ne 0) { throw "Failed to upgrade pip/setuptools/wheel." }

$requirements = Join-Path $projectRoot "scripts\python_requirements_ag3.txt"
& $venvPython -m pip install --upgrade -r $requirements
if ($LASTEXITCODE -ne 0) { throw "Failed to install the Ag3 Python requirements." }

& $venvPython -c "import sys,malariagen_data,pandas,numpy,yaml; print('Python:',sys.version.split()[0]); print('malariagen_data:',getattr(malariagen_data,'__version__','unknown')); print('pandas:',pandas.__version__); print('numpy:',numpy.__version__); print('Targeted Ag3 environment ready.')"
if ($LASTEXITCODE -ne 0) { throw "Python environment import test failed." }

New-Item -ItemType Directory -Force -Path (Join-Path $projectRoot "metadata") | Out-Null
$pointer = Join-Path $projectRoot "metadata\ag3_python_path.txt"
Set-Content -Path $pointer -Value $venvPython -Encoding UTF8
& $venvPython -m pip freeze | Out-File -Encoding utf8 (Join-Path $projectRoot "metadata\python_environment_freeze.txt")

Write-Host ""
Write-Host "Targeted Ag3 Python environment is ready."
Write-Host "Interpreter: $venvPython"
Write-Host "Pointer written to: $pointer"
Write-Host "Next (no network/data download):"
Write-Host "  & `"$venvPython`" scripts\acquire_ag3_targeted.py --dry-run"
