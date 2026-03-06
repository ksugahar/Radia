#==============================================================================
# Build.ps1 - Build Radia with MSVC + Intel MKL
#
# Builds both radia.pyd and radia_ngsolve.pyd using MSVC compiler
# with Intel MKL for BLAS/LAPACK operations.
#
# Usage:
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Rebuild
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Test
#
# Options:
#   -Rebuild  Clean build directory before building
#   -Test     Run import test after build
#   -NoOpenMP Disable OpenMP (for debugging)
#==============================================================================

param(
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$NoOpenMP
)

# Forward to BuildMSVC.ps1
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$args = @()
if ($Rebuild) { $args += "-Rebuild" }
if ($Test) { $args += "-Test" }
if ($NoOpenMP) { $args += "-NoOpenMP" }

& powershell.exe -ExecutionPolicy Bypass -File "$ScriptDir\BuildMSVC.ps1" @args
exit $LASTEXITCODE
