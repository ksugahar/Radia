#==============================================================================
# BuildRadiaInternal.ps1 - DEPRECATED
#
# This script is DEPRECATED. Use Build.ps1 or BuildMSVC.ps1 instead.
#
# Usage: powershell.exe -ExecutionPolicy Bypass -File Build.ps1 [-Rebuild]
#==============================================================================

param(
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$Verbose
)

Write-Host ""
Write-Host "========================================" -ForegroundColor Yellow
Write-Host "  DEPRECATED: Use Build.ps1" -ForegroundColor Yellow
Write-Host "========================================" -ForegroundColor Yellow
Write-Host ""
Write-Host "Redirecting to BuildMSVC.ps1..." -ForegroundColor Cyan
Write-Host ""

# Forward to BuildMSVC.ps1
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$args = @()
if ($Rebuild) { $args += "-Rebuild" }
if ($Test) { $args += "-Test" }

& powershell.exe -ExecutionPolicy Bypass -File "$ScriptDir\BuildMSVC.ps1" @args
exit $LASTEXITCODE
