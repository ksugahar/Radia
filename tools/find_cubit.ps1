<#
.SYNOPSIS
    Find the latest Coreform Cubit install. Dot-source from build scripts.

.DESCRIPTION
    Centralizes Cubit install discovery so Build.ps1, cubit_plugin
    builds, and CI no longer hardcode a specific version. Replaces
    the prior pattern:
        $env:CUBIT_DIR = 'C:\Program Files\Coreform Cubit 2025.3\cmake'
    with:
        . tools/find_cubit.ps1
        $env:CUBIT_DIR = $CubitCmakeDir   # always points at latest

    Discovery rules (first hit wins):
      1. $env:CUBIT_INSTALL_DIR if set and exists
      2. Latest 'C:\Program Files\Coreform Cubit *' dir
         (version-aware sort: 2025.12 > 2025.3)
      3. $null + warning if nothing found

    When dot-sourced, sets script-scoped variables:
      $CubitInstallDir   = e.g. 'C:\Program Files\Coreform Cubit 2025.12'
      $CubitBinDir       = "$CubitInstallDir\bin"
      $CubitCmakeDir     = "$CubitInstallDir\cmake"
      $CubitVersion      = e.g. '2025.12'

    When run as a standalone script, prints those values one-per-line
    (KEY=VALUE) so a .bat caller can `FOR /F ... DO` parse them.

.EXAMPLE
    # PowerShell consumer (Build.ps1 style):
    . "$PSScriptRoot\..\tools\find_cubit.ps1"
    if (-not $CubitInstallDir) { Write-Error "Cubit not installed"; exit 1 }
    Write-Host "Using Cubit: $CubitVersion at $CubitBinDir"

.EXAMPLE
    # Batch consumer (build_cubit_only.bat style):
    @echo off
    for /f "tokens=1,2 delims==" %%A in ('powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0..\tools\find_cubit.ps1"') do (
        set "%%A=%%B"
    )
    rem now %CubitInstallDir%, %CubitBinDir%, %CubitCmakeDir%, %CubitVersion% are set

.EXAMPLE
    # Override for an off-default install:
    $env:CUBIT_INSTALL_DIR = 'D:\Tools\Cubit2025.12'
    . tools/find_cubit.ps1
#>

function Find-CubitInstall {
    # Honor explicit env-var override first.
    if ($env:CUBIT_INSTALL_DIR -and (Test-Path $env:CUBIT_INSTALL_DIR)) {
        return $env:CUBIT_INSTALL_DIR
    }
    # Scan C:\Program Files for the highest-version install.
    $base = 'C:\Program Files'
    $cands = Get-ChildItem -Path $base -Filter 'Coreform Cubit*' `
        -Directory -ErrorAction SilentlyContinue
    if (-not $cands) { return $null }
    # Version-aware sort: parse "Coreform Cubit X.Y" -> (X, Y).
    $sorted = $cands | Sort-Object @{
        Expression = {
            $verStr = ($_.Name -replace '^Coreform Cubit\s*', '').Trim()
            $parts  = $verStr -split '\.'
            $major  = 0; $minor = 0
            [int]::TryParse($parts[0], [ref]$major) | Out-Null
            if ($parts.Length -gt 1) {
                [int]::TryParse($parts[1], [ref]$minor) | Out-Null
            }
            # 5-digit minor handles future "Cubit 2026.12345" without overflow.
            $major * 100000 + $minor
        }
        Descending = $true
    }
    return $sorted[0].FullName
}

# ===== Compute discovery results (always runs, dot-sourced or scripted) =====
$script:CubitInstallDir = Find-CubitInstall
if ($script:CubitInstallDir) {
    $script:CubitBinDir   = Join-Path $script:CubitInstallDir 'bin'
    $script:CubitCmakeDir = Join-Path $script:CubitInstallDir 'cmake'
    $script:CubitVersion  = (Split-Path -Leaf $script:CubitInstallDir) `
        -replace '^Coreform Cubit\s*', ''
} else {
    $script:CubitBinDir   = $null
    $script:CubitCmakeDir = $null
    $script:CubitVersion  = $null
}

# ===== Script mode: emit KEY=VALUE for batch / shell consumers =====
if ($MyInvocation.InvocationName -ne '.') {
    if ($script:CubitInstallDir) {
        "CubitInstallDir=$($script:CubitInstallDir)"
        "CubitBinDir=$($script:CubitBinDir)"
        "CubitCmakeDir=$($script:CubitCmakeDir)"
        "CubitVersion=$($script:CubitVersion)"
    } else {
        Write-Error "No Coreform Cubit install found under C:\Program Files. Set CUBIT_INSTALL_DIR to override."
        exit 1
    }
}
