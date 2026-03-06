# PowerShell script to build and publish package to PyPI

# Clean dist directory
if (Test-Path dist) {
	Remove-Item -Recurse -Force dist
}

# Set UTF-8 encoding
$env:PYTHONIOENCODING = "utf-8"

# Build package
Write-Host "Building package..." -ForegroundColor Cyan
python -m build

if ($LASTEXITCODE -ne 0) {
	Write-Host "Build failed!" -ForegroundColor Red
	exit 1
}

# Upload to PyPI
Write-Host "Uploading to PyPI..." -ForegroundColor Cyan

# SECURITY: Use environment variable for API token
# Set the token before running: $env:PYPI_TOKEN = "your-token-here"
if (-not $env:PYPI_TOKEN) {
	Write-Host "ERROR: PYPI_TOKEN environment variable not set" -ForegroundColor Red
	Write-Host "Set it using: `$env:PYPI_TOKEN = `"your-token-here`"" -ForegroundColor Yellow
	exit 1
}

python -m twine upload dist/* --username __token__ --password $env:PYPI_TOKEN

if ($LASTEXITCODE -eq 0) {
	Write-Host "Upload successful!" -ForegroundColor Green
} else {
	Write-Host "Upload failed!" -ForegroundColor Red
}

# Wait for key press before exiting
Write-Host ""
Write-Host "Press any key to exit..." -ForegroundColor Cyan
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
