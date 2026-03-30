Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$python = Join-Path $PSScriptRoot "..\.venv\Scripts\python.exe"
$python = [System.IO.Path]::GetFullPath($python)

if (-not (Test-Path $python)) {
    throw "Missing .venv\Scripts\python.exe. Create the virtualenv and install requirements-desktop.txt first."
}

$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))
Push-Location $projectRoot
try {
    & $python -m PyInstaller --noconfirm --clean .\VinaDockStudio.spec
}
finally {
    Pop-Location
}
