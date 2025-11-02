# Build script for creating PyPI distribution packages
#
# This script builds both source distribution (sdist) and wheel distribution (bdist_wheel)
# for uploading to PyPI.
#
# Prerequisites:
# 1. Run Build.ps1 first to build radia.pyd
# 2. (Optional) Run Build_NGSolve.ps1 to build rad_ngsolve.pyd
# 3. Install build tools: pip install build twine

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Radia PyPI Package Build Script" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Check if radia.pyd exists
if (-not (Test-Path "dist\radia.pyd")) {
	Write-Host "ERROR: dist\radia.pyd not found" -ForegroundColor Red
	Write-Host "Please run Build.ps1 first to build the core module" -ForegroundColor Yellow
	exit 1
}

Write-Host "[OK] Found dist\radia.pyd" -ForegroundColor Green

# Check for rad_ngsolve.pyd (optional)
if (Test-Path "build\Release\rad_ngsolve.pyd") {
	Write-Host "[OK] Found build\Release\rad_ngsolve.pyd" -ForegroundColor Green
} else {
	Write-Host "[INFO] rad_ngsolve.pyd not found (optional)" -ForegroundColor Yellow
}

# Clean previous builds
Write-Host ""
Write-Host "Cleaning previous builds..." -ForegroundColor Cyan

$cleanPaths = @(
	"build\lib",
	"build\bdist.*",
	"src\python\radia.egg-info",
	"radia.egg-info"
)

foreach ($path in $cleanPaths) {
	if (Test-Path $path) {
		Remove-Item $path -Recurse -Force -ErrorAction SilentlyContinue
	}
}

# Remove wheel and tar.gz files
Get-ChildItem "dist" -Filter "*.whl" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem "dist" -Filter "*.tar.gz" -ErrorAction SilentlyContinue | Remove-Item -Force

# Build source distribution
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Building Source Distribution (sdist)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

python -m build --sdist

if ($LASTEXITCODE -ne 0) {
	Write-Host "ERROR: Source distribution build failed" -ForegroundColor Red
	exit 1
}

Write-Host "[OK] Source distribution created" -ForegroundColor Green

# Build wheel distribution
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Building Wheel Distribution (bdist_wheel)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

python -m build --wheel

if ($LASTEXITCODE -ne 0) {
	Write-Host "ERROR: Wheel distribution build failed" -ForegroundColor Red
	exit 1
}

Write-Host "[OK] Wheel distribution created" -ForegroundColor Green

# Show results
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Build Complete" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

Write-Host "Distribution files created in dist/:" -ForegroundColor Cyan
Get-ChildItem "dist" -Filter "*.whl" | Select-Object Name
Get-ChildItem "dist" -Filter "*.tar.gz" | Select-Object Name

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Next Steps" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "To test the package locally:" -ForegroundColor Yellow
Write-Host "  pip install dist\radia-ngsolve-1.0.0-py3-none-any.whl" -ForegroundColor White
Write-Host ""
Write-Host "To upload to Test PyPI:" -ForegroundColor Yellow
Write-Host "  python -m twine upload --repository testpypi dist/*" -ForegroundColor White
Write-Host ""
Write-Host "To upload to PyPI:" -ForegroundColor Yellow
Write-Host "  python -m twine upload dist/*" -ForegroundColor White
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
