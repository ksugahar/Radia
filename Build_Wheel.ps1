# ============================================================================
# Build_Wheel.ps1 - Build and verify Radia wheel for PyPI
# ============================================================================
#
# Usage:
#   CI (default):     .\Build_Wheel.ps1 -DryRun    # build + verify only
#   Local publish:    $env:PYPI_TOKEN = "pypi-..."
#                     .\Build_Wheel.ps1             # build + verify + upload
#
# What this script does:
#   1. Clean previous build artifacts
#   2. Build wheel (no sdist - source build requires CMake)
#   3. Remove any accidentally bundled MKL DLLs from wheel
#   4. Repack wheel with correct platform tag (cp312-cp312-win_amd64)
#   5. Verify with twine check
#   6. Upload to PyPI (skipped with -DryRun; CI uses Trusted Publishers instead)
#
# DLL policy: radia_motor_rom.dll is the public standalone C ABI and is
# bundled. Third-party runtime DLLs such as MKL are not bundled.
# ============================================================================

param(
    [switch]$DryRun,       # Build and verify only, skip upload
    [switch]$SkipBuild     # Skip build, use existing wheel in dist/
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$env:PYTHONIOENCODING = "utf-8"

# Configuration
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$DistDir = Join-Path $ProjectRoot "dist"
$BuildDir = Join-Path $ProjectRoot "build"
$BuildLibRadia = Join-Path (Join-Path $BuildDir "lib") "radia"
$PlatformTag = "cp312-cp312-win_amd64"

Write-Host "============================================" -ForegroundColor Cyan
Write-Host "  Radia Wheel Builder" -ForegroundColor Cyan
Write-Host "============================================" -ForegroundColor Cyan
Write-Host ""

# ----------------------------------------------------------------------------
# Step 0: Check prerequisites
# ----------------------------------------------------------------------------
if (-not $DryRun -and -not $env:PYPI_TOKEN) {
    Write-Host "ERROR: PYPI_TOKEN environment variable not set" -ForegroundColor Red
    Write-Host 'Set it using: $env:PYPI_TOKEN = "pypi-..."' -ForegroundColor Yellow
    Write-Host "Or use -DryRun to build and verify without uploading." -ForegroundColor Yellow
    exit 1
}

# Read version from pyproject.toml
$PyprojectPath = Join-Path $ProjectRoot "pyproject.toml"
$VersionLine = Select-String -Path $PyprojectPath -Pattern '^version\s*=\s*"(.+)"' | Select-Object -First 1
if ($VersionLine) {
    $Version = $VersionLine.Matches[0].Groups[1].Value
} else {
    Write-Host "ERROR: Could not read version from pyproject.toml" -ForegroundColor Red
    exit 1
}
Write-Host "Version: $Version" -ForegroundColor White
Write-Host "Platform: $PlatformTag" -ForegroundColor White
Write-Host ""

if (-not $SkipBuild) {
    # ------------------------------------------------------------------------
    # Step 1: Clean previous build artifacts
    # ------------------------------------------------------------------------
    Write-Host "[1/6] Cleaning build artifacts..." -ForegroundColor Cyan

    if (Test-Path $DistDir) {
        Remove-Item -Recurse -Force $DistDir
        Write-Host "  Removed dist/" -ForegroundColor Gray
    }

    # Clean MKL DLLs from build/lib/radia/ to prevent accidental bundling
    if (Test-Path $BuildLibRadia) {
        $DLLs = @(Get-ChildItem -Path $BuildLibRadia -Filter "*.dll" -ErrorAction SilentlyContinue)
        if ($DLLs.Count -gt 0) {
            $DLLs | Remove-Item -Force
            Write-Host "  Removed $($DLLs.Count) DLL(s) from build/lib/radia/" -ForegroundColor Yellow
        }
    }

    Write-Host "  Done." -ForegroundColor Green
    Write-Host ""

    # ------------------------------------------------------------------------
    # Step 2: Build wheel (no sdist)
    # ------------------------------------------------------------------------
    Write-Host "[2/6] Building wheel..." -ForegroundColor Cyan

    # Fail on the first build error. CI retries hide reproducible packaging
    # failures and delay diagnosis; the captured log preserves the subprocess
    # output needed to investigate transient infrastructure failures separately.
    Push-Location $ProjectRoot
    try {
        $logFile = Join-Path $env:TEMP "radia_wheel_build.log"
        python -m build --wheel 2>&1 | Tee-Object -FilePath $logFile
        if ($LASTEXITCODE -ne 0) {
            Write-Host "" -ForegroundColor Red
            Write-Host "ERROR: Wheel build failed." -ForegroundColor Red
            Write-Host "Build log (tail):" -ForegroundColor Red
            if (Test-Path $logFile) {
                Get-Content $logFile -Tail 40 | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
            }
            exit 1
        }
    } finally {
        Pop-Location
    }

    Write-Host "  Done." -ForegroundColor Green
    Write-Host ""

    # ------------------------------------------------------------------------
    # Step 3: Remove third-party DLLs from wheel (safety check)
    # ------------------------------------------------------------------------
    Write-Host "[3/6] Checking wheel for bundled DLLs..." -ForegroundColor Cyan

    $OrigWheel = Get-ChildItem -Path $DistDir -Filter "radia-*.whl" | Select-Object -First 1
    if (-not $OrigWheel) {
        Write-Host "ERROR: No wheel found in dist/" -ForegroundColor Red
        exit 1
    }

    # Extract wheel to temp directory for inspection and repack
    $TempDir = Join-Path $env:TEMP "radia_wheel_repack"
    if (Test-Path $TempDir) { Remove-Item -Recurse -Force $TempDir }
    New-Item -ItemType Directory -Path $TempDir | Out-Null

    # Rename .whl to .zip for extraction
    $ZipPath = Join-Path $env:TEMP "radia_wheel_temp.zip"
    Copy-Item $OrigWheel.FullName $ZipPath -Force
    Expand-Archive -Path $ZipPath -DestinationPath $TempDir -Force
    Remove-Item $ZipPath

    # Keep the Radia-owned standalone C ABI, but remove third-party runtime
    # DLLs (especially MKL) that may have leaked in from a local build tree.
    $AllowedDllNames = @("radia_motor_rom.dll")
    $BundledDLLs = @(Get-ChildItem -Path $TempDir -Filter "*.dll" -Recurse -ErrorAction SilentlyContinue)
    $RemovedDLLs = @($BundledDLLs | Where-Object { $_.Name -notin $AllowedDllNames })
    if ($RemovedDLLs.Count -gt 0) {
        Write-Host "  WARNING: Found third-party DLLs - removing:" -ForegroundColor Yellow
        foreach ($dll in $RemovedDLLs) {
            Write-Host "    - $($dll.Name)" -ForegroundColor Yellow
            Remove-Item $dll.FullName -Force
        }
    } else {
        Write-Host "  No third-party DLLs found (good)." -ForegroundColor Green
    }

    Write-Host ""

    # ------------------------------------------------------------------------
    # Step 4: Repack wheel with correct platform tag
    # ------------------------------------------------------------------------
    Write-Host "[4/6] Repacking wheel with platform tag: $PlatformTag ..." -ForegroundColor Cyan

    # Fix WHEEL metadata
    $WheelMetaPath = Get-ChildItem -Path $TempDir -Filter "WHEEL" -Recurse | Select-Object -First 1
    if ($WheelMetaPath) {
        $WheelContent = Get-Content $WheelMetaPath.FullName -Raw
        # Replace Tag line
        $WheelContent = $WheelContent -replace 'Tag: .+', "Tag: $PlatformTag"
        # Ensure Root-Is-Purelib is false
        $WheelContent = $WheelContent -replace 'Root-Is-Purelib: true', 'Root-Is-Purelib: false'
        Set-Content -Path $WheelMetaPath.FullName -Value $WheelContent -NoNewline
    }

    # Also fix RECORD file to remove entries for DLLs that were removed.
    $RecordPath = Get-ChildItem -Path $TempDir -Filter "RECORD" -Recurse | Select-Object -First 1
    if ($RecordPath -and $RemovedDLLs.Count -gt 0) {
        $RemovedNames = @($RemovedDLLs | ForEach-Object { [regex]::Escape($_.Name) })
        $RemovedPattern = '/(' + ($RemovedNames -join '|') + '),'
        $RecordLines = Get-Content $RecordPath.FullName | Where-Object { $_ -notmatch $RemovedPattern }
        Set-Content -Path $RecordPath.FullName -Value ($RecordLines -join "`n") -NoNewline
    }

    # Remove old wheel and create new one with correct name
    Remove-Item $OrigWheel.FullName -Force
    $NewWheelName = "radia-$Version-$PlatformTag.whl"
    $NewWheelPath = Join-Path $DistDir $NewWheelName

    # Repack as zip
    $NewZipPath = Join-Path $env:TEMP "radia_wheel_new.zip"
    Compress-Archive -Path (Join-Path $TempDir "*") -DestinationPath $NewZipPath -Force
    Move-Item $NewZipPath $NewWheelPath -Force

    # Cleanup
    Remove-Item -Recurse -Force $TempDir

    Write-Host "  Created: $NewWheelName" -ForegroundColor Green
    Write-Host ""

} else {
    # SkipBuild: use existing wheel
    $NewWheelPath = Get-ChildItem -Path $DistDir -Filter "radia-*-$PlatformTag.whl" | Select-Object -First 1
    if (-not $NewWheelPath) {
        Write-Host "ERROR: No wheel with tag $PlatformTag found in dist/" -ForegroundColor Red
        exit 1
    }
    $NewWheelPath = $NewWheelPath.FullName
    Write-Host "Using existing wheel: $NewWheelPath" -ForegroundColor Cyan
}

# ----------------------------------------------------------------------------
# Step 5: Verify wheel
# ----------------------------------------------------------------------------
Write-Host "[5/6] Verifying wheel..." -ForegroundColor Cyan

# twine check
python -m twine check $NewWheelPath
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: twine check failed!" -ForegroundColor Red
    exit 1
}

# Content summary
Write-Host ""
Write-Host "  Wheel contents summary:" -ForegroundColor White

$InspectTempDir = Join-Path $env:TEMP "radia_wheel_inspect"
if (Test-Path $InspectTempDir) { Remove-Item -Recurse -Force $InspectTempDir }
$InspectZip = Join-Path $env:TEMP "radia_inspect.zip"
Copy-Item $NewWheelPath $InspectZip -Force
Expand-Archive -Path $InspectZip -DestinationPath $InspectTempDir -Force
Remove-Item $InspectZip

$AllFiles = @(Get-ChildItem -Path $InspectTempDir -Recurse -File)
$PyFiles = @($AllFiles | Where-Object { $_.Extension -eq ".py" })
$PydFiles = @($AllFiles | Where-Object { $_.Extension -eq ".pyd" })
$DllFiles = @($AllFiles | Where-Object { $_.Extension -eq ".dll" })
$WheelSize = (Get-Item $NewWheelPath).Length

Write-Host "    Total files: $($AllFiles.Count)" -ForegroundColor White
Write-Host "    .py files:   $($PyFiles.Count)" -ForegroundColor White
Write-Host "    .pyd files:  $($PydFiles.Count)" -ForegroundColor White
Write-Host "    .dll files:  $($DllFiles.Count)" -ForegroundColor White
Write-Host "    Wheel size:  $([math]::Round($WheelSize / 1MB, 2)) MB" -ForegroundColor White

Remove-Item -Recurse -Force $InspectTempDir

if (-not ($DllFiles | Where-Object { $_.Name -eq "radia_motor_rom.dll" })) {
    Write-Host ""
    Write-Host "  ERROR: Wheel is missing radia_motor_rom.dll." -ForegroundColor Red
    exit 1
}

$UnexpectedDLLs = @($DllFiles | Where-Object { $_.Name -ne "radia_motor_rom.dll" })
if ($UnexpectedDLLs.Count -gt 0) {
    Write-Host ""
    Write-Host "  ERROR: Wheel contains unexpected third-party DLL files." -ForegroundColor Red
    Write-Host "  Unexpected DLLs:" -ForegroundColor Red
    foreach ($dll in $UnexpectedDLLs) { Write-Host "    - $($dll.Name)" -ForegroundColor Red }
    exit 1
}

Write-Host ""
Write-Host "  Verification PASSED" -ForegroundColor Green
Write-Host ""

# ----------------------------------------------------------------------------
# Step 6: Upload to PyPI
# ----------------------------------------------------------------------------
if ($DryRun) {
    Write-Host "[6/6] DRY RUN - Skipping upload." -ForegroundColor Yellow
    Write-Host "  Wheel ready at: $NewWheelPath" -ForegroundColor White
    Write-Host ""
    Write-Host "To upload manually:" -ForegroundColor Cyan
    Write-Host "  python -m twine upload `"$NewWheelPath`" --username __token__ --password `$env:PYPI_TOKEN" -ForegroundColor White
} else {
    Write-Host "[6/6] Uploading to PyPI..." -ForegroundColor Cyan

    python -m twine upload $NewWheelPath --username __token__ --password $env:PYPI_TOKEN

    if ($LASTEXITCODE -eq 0) {
        Write-Host ""
        Write-Host "============================================" -ForegroundColor Green
        Write-Host "  Upload successful!  radia $Version" -ForegroundColor Green
        Write-Host "  https://pypi.org/project/radia/$Version/" -ForegroundColor Green
        Write-Host "============================================" -ForegroundColor Green
    } else {
        Write-Host ""
        Write-Host "Upload failed. To retry manually:" -ForegroundColor Red
        Write-Host "  python -m twine upload `"$NewWheelPath`" --username __token__ --password `$env:PYPI_TOKEN" -ForegroundColor Yellow
    }
}

Write-Host ""
# Skip interactive prompt in non-interactive mode (e.g., CI)
if ([Environment]::UserInteractive -and -not $env:CI) {
    Write-Host "Press any key to exit..." -ForegroundColor Cyan
    $null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
}
