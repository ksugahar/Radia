# Tetrahedron Nonlinear Material Benchmark
#
# Runs benchmark_tetra_netgen.py with various mesh sizes
# Generates results in tetrahedron/lu/ and tetrahedron/bicgstab/ directories
#
# Usage:
#   .\run_tetra_benchmark.ps1           # Run both LU and BiCGSTAB
#   .\run_tetra_benchmark.ps1 -LU       # Run LU only
#   .\run_tetra_benchmark.ps1 -BiCGSTAB # Run BiCGSTAB only
#
# Author: Radia Development Team
# Date: 2025-12-06

param(
    [switch]$LU,
    [switch]$BiCGSTAB,
    [string[]]$MaxH = @("0.4", "0.3", "0.25", "0.2")
)

$ErrorActionPreference = "Stop"

# Get script directory
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$BenchmarkScript = Join-Path $ScriptDir "benchmark_tetra_netgen.py"

# Check if benchmark script exists
if (-not (Test-Path $BenchmarkScript)) {
    Write-Error "Benchmark script not found: $BenchmarkScript"
    exit 1
}

# Header
Write-Host "=" * 70
Write-Host "RADIA TETRAHEDRON NONLINEAR BENCHMARK"
Write-Host "=" * 70
Write-Host ""
Write-Host "Benchmark script: $BenchmarkScript"
Write-Host "maxh values: $($MaxH -join ', ')"
Write-Host ""

# If neither -LU nor -BiCGSTAB specified, run both
$RunLU = $LU -or (-not $LU -and -not $BiCGSTAB)
$RunBiCGSTAB = $BiCGSTAB -or (-not $LU -and -not $BiCGSTAB)

# Run benchmarks
if ($RunLU) {
    Write-Host "Running LU solver benchmark..."
    Write-Host "-" * 70
    & python $BenchmarkScript --lu $MaxH
    Write-Host ""
}

if ($RunBiCGSTAB) {
    Write-Host "Running BiCGSTAB solver benchmark..."
    Write-Host "-" * 70
    & python $BenchmarkScript --bicgstab $MaxH
    Write-Host ""
}

# Summary
Write-Host "=" * 70
Write-Host "BENCHMARK COMPLETE"
Write-Host "=" * 70
Write-Host ""
Write-Host "Results saved to:"
if ($RunLU) {
    Write-Host "  - tetrahedron/lu/*.json"
}
if ($RunBiCGSTAB) {
    Write-Host "  - tetrahedron/hacapk/*.json"
}
Write-Host ""
