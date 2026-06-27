# PowerShell script to run all 18 simulation files (Sphere_3D + Cylinder_3D)
# Each model: order=2,3,4 x 3 methods (Refine_all_elements, Refine_with_zz_estimator, metric_based)

$ErrorActionPreference = "Continue"
$startTime = Get-Date

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "Starting All Simulations (Sphere_3D + Cylinder_3D)" -ForegroundColor Cyan
Write-Host "Start Time: $startTime" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

# Base directory
$baseDir = "S:\NGSolve\NGSolve\2025_12_14_Kelvin変換\AdaptiveMesh\Omega_ReducedOmega"

# Define all simulation files with literal paths
$simulations = @(
    # Sphere_3D - order=2
    @{Name="Sphere_3D order=2 Refine_all_elements"; Dir="Sphere_3D\order=2\Refine_all_elements"; File="Sphere_3D_adaptive_with_Kelvin.py"},
    @{Name="Sphere_3D order=2 Refine_with_zz_estimator"; Dir="Sphere_3D\order=2\Refine_with_zz_estimator"; File="Sphere_3D_adaptive_with_Kelvin.py"},
    @{Name="Sphere_3D order=2 metric_based"; Dir="Sphere_3D\order=2\metric_based"; File="Sphere_3D_metric_with_Kelvin.py"},

    # Sphere_3D - order=3
    @{Name="Sphere_3D order=3 Refine_all_elements"; Dir="Sphere_3D\order=3\Refine_all_elements"; File="Sphere_3D_adaptive_with_Kelvin.py"},
    @{Name="Sphere_3D order=3 Refine_with_zz_estimator"; Dir="Sphere_3D\order=3\Refine_with_zz_estimator"; File="Sphere_3D_adaptive_with_Kelvin.py"},
    @{Name="Sphere_3D order=3 metric_based"; Dir="Sphere_3D\order=3\metric_based"; File="Sphere_3D_metric_with_Kelvin.py"},

    # Sphere_3D - order=4
    @{Name="Sphere_3D order=4 Refine_all_elements"; Dir="Sphere_3D\order=4\Refine_all_elements"; File="Sphere_3D_adaptive_with_Kelvin.py"},
    @{Name="Sphere_3D order=4 Refine_with_zz_estimator"; Dir="Sphere_3D\order=4\Refine_with_zz_estimator"; File="Sphere_3D_adaptive_with_Kelvin.py"},
    @{Name="Sphere_3D order=4 metric_based"; Dir="Sphere_3D\order=4\metric_based"; File="Sphere_3D_metric_with_Kelvin.py"},

    # Cylinder_3D - order=2
    @{Name="Cylinder_3D order=2 Refine_all_elements"; Dir="Cylinder_3D\order=2\Refine_all_elements"; File="Cylinder_3D_adaptive_with_Kelvin.py"},
    @{Name="Cylinder_3D order=2 Refine_with_zz_estimator"; Dir="Cylinder_3D\order=2\Refine_with_zz_estimator"; File="Cylinder_3D_adaptive_with_Kelvin.py"},
    @{Name="Cylinder_3D order=2 metric_based"; Dir="Cylinder_3D\order=2\metric_based"; File="Cylinder_3D_metric_with_Kelvin.py"},

    # Cylinder_3D - order=3
    @{Name="Cylinder_3D order=3 Refine_all_elements"; Dir="Cylinder_3D\order=3\Refine_all_elements"; File="Cylinder_3D_adaptive_with_Kelvin.py"},
    @{Name="Cylinder_3D order=3 Refine_with_zz_estimator"; Dir="Cylinder_3D\order=3\Refine_with_zz_estimator"; File="Cylinder_3D_adaptive_with_Kelvin.py"},
    @{Name="Cylinder_3D order=3 metric_based"; Dir="Cylinder_3D\order=3\metric_based"; File="Cylinder_3D_metric_with_Kelvin.py"},

    # Cylinder_3D - order=4
    @{Name="Cylinder_3D order=4 Refine_all_elements"; Dir="Cylinder_3D\order=4\Refine_all_elements"; File="Cylinder_3D_adaptive_with_Kelvin.py"},
    @{Name="Cylinder_3D order=4 Refine_with_zz_estimator"; Dir="Cylinder_3D\order=4\Refine_with_zz_estimator"; File="Cylinder_3D_adaptive_with_Kelvin.py"},
    @{Name="Cylinder_3D order=4 metric_based"; Dir="Cylinder_3D\order=4\metric_based"; File="Cylinder_3D_metric_with_Kelvin.py"}
)

$totalSimulations = $simulations.Count
$completedSimulations = 0
$failedSimulations = @()

Write-Host "`nTotal simulations to run: $totalSimulations" -ForegroundColor Yellow
Write-Host ""

foreach ($sim in $simulations) {
    $simStartTime = Get-Date
    $simIndex = $completedSimulations + 1
    $simDir = Join-Path $baseDir $sim.Dir
    $fileName = $sim.File

    Write-Host "------------------------------------------------------------" -ForegroundColor Gray
    Write-Host "[$simIndex/$totalSimulations] $($sim.Name)" -ForegroundColor Green
    Write-Host "Directory: $simDir" -ForegroundColor Gray
    Write-Host "File: $fileName" -ForegroundColor Gray
    Write-Host "Started: $simStartTime" -ForegroundColor Gray
    Write-Host ""

    # Change to simulation directory and run
    Push-Location -LiteralPath $simDir
    try {
        python $fileName
        $exitCode = $LASTEXITCODE

        if ($exitCode -eq 0) {
            Write-Host ""
            Write-Host "  SUCCESS" -ForegroundColor Green
        } else {
            Write-Host ""
            Write-Host "  FAILED (exit code: $exitCode)" -ForegroundColor Red
            $failedSimulations += $sim.Name
        }
    }
    catch {
        Write-Host ""
        Write-Host "  ERROR: $_" -ForegroundColor Red
        $failedSimulations += $sim.Name
    }
    finally {
        Pop-Location
    }

    $simEndTime = Get-Date
    $simDuration = $simEndTime - $simStartTime
    Write-Host "  Duration: $($simDuration.ToString('hh\:mm\:ss'))" -ForegroundColor Gray

    $completedSimulations++
    $remainingSimulations = $totalSimulations - $completedSimulations
    Write-Host "  Progress: $completedSimulations/$totalSimulations completed, $remainingSimulations remaining" -ForegroundColor Gray
    Write-Host ""
}

$endTime = Get-Date
$totalDuration = $endTime - $startTime

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "All Simulations Completed" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Summary:" -ForegroundColor Yellow
Write-Host "  Total simulations: $totalSimulations"
Write-Host "  Successful: $($totalSimulations - $failedSimulations.Count)" -ForegroundColor Green
Write-Host "  Failed: $($failedSimulations.Count)" -ForegroundColor $(if ($failedSimulations.Count -gt 0) { "Red" } else { "Green" })

if ($failedSimulations.Count -gt 0) {
    Write-Host ""
    Write-Host "Failed simulations:" -ForegroundColor Red
    foreach ($failed in $failedSimulations) {
        Write-Host "  - $failed" -ForegroundColor Red
    }
}

Write-Host ""
Write-Host "Start Time: $startTime" -ForegroundColor Gray
Write-Host "End Time: $endTime" -ForegroundColor Gray
Write-Host "Total Duration: $($totalDuration.ToString('hh\:mm\:ss'))" -ForegroundColor Yellow
Write-Host ""
