# IEEE Access Paper Compilation Script
# Usage: .\compile.ps1 [-Clean] [-Draft] [-Final] [-View]
#
# Options:
#   -Clean  : Remove all auxiliary files before compilation
#   -Draft  : Quick compilation (single pass, no bibliography)
#   -Final  : Full compilation (pdflatex -> bibtex -> pdflatex x2)
#   -View   : Open PDF after compilation
#
# Requirements:
#   - MiKTeX or TeX Live with pdflatex and bibtex
#   - SumatraPDF or other PDF viewer (optional)

param(
    [switch]$Clean,
    [switch]$Draft,
    [switch]$Final,
    [switch]$View,
    [string]$TexFile = "urn_paper"
)

# Set working directory to script location
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

# Colors for output
function Write-Success { param($msg) Write-Host $msg -ForegroundColor Green }
function Write-Info { param($msg) Write-Host $msg -ForegroundColor Cyan }
function Write-Warn { param($msg) Write-Host $msg -ForegroundColor Yellow }
function Write-Err { param($msg) Write-Host $msg -ForegroundColor Red }

# Check if TeX file exists
if (-not (Test-Path "$TexFile.tex")) {
    # Fall back to access.tex if urn_paper.tex does not exist
    if (Test-Path "access.tex") {
        $TexFile = "access"
        Write-Warn "Using access.tex (urn_paper.tex not found)"
    } else {
        Write-Err "Error: No .tex file found!"
        exit 1
    }
}

Write-Info "=== IEEE Access Paper Compilation ==="
Write-Info "TeX file: $TexFile.tex"
Write-Info "Working directory: $ScriptDir"

# Clean auxiliary files
if ($Clean) {
    Write-Info "Cleaning auxiliary files..."
    $extensions = @("aux", "bbl", "blg", "log", "out", "toc", "lof", "lot",
                    "fls", "fdb_latexmk", "synctex.gz", "nav", "snm", "vrb")
    foreach ($ext in $extensions) {
        Remove-Item -Path "*.$ext" -Force -ErrorAction SilentlyContinue
    }
    Write-Success "Cleaned auxiliary files"
}

# Find pdflatex
$pdflatex = Get-Command pdflatex -ErrorAction SilentlyContinue
if (-not $pdflatex) {
    # Try common installation paths
    $paths = @(
        "C:\texlive\2024\bin\windows\pdflatex.exe",
        "C:\texlive\2023\bin\windows\pdflatex.exe",
        "C:\Program Files\MiKTeX\miktex\bin\x64\pdflatex.exe",
        "$env:LOCALAPPDATA\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
    )
    foreach ($p in $paths) {
        if (Test-Path $p) {
            $pdflatex = $p
            break
        }
    }
}

if (-not $pdflatex) {
    Write-Err "Error: pdflatex not found. Please install MiKTeX or TeX Live."
    exit 1
}

Write-Info "Using pdflatex: $pdflatex"

# Find bibtex
$bibtex = Get-Command bibtex -ErrorAction SilentlyContinue
if (-not $bibtex) {
    $btpaths = @(
        "C:\texlive\2024\bin\windows\bibtex.exe",
        "C:\texlive\2023\bin\windows\bibtex.exe",
        "C:\Program Files\MiKTeX\miktex\bin\x64\bibtex.exe",
        "$env:LOCALAPPDATA\Programs\MiKTeX\miktex\bin\x64\bibtex.exe"
    )
    foreach ($p in $btpaths) {
        if (Test-Path $p) {
            $bibtex = $p
            break
        }
    }
}

# Compilation function
function Invoke-PdfLatex {
    param([string]$Pass)
    Write-Info "Running pdflatex ($Pass)..."
    $proc = Start-Process -FilePath $pdflatex -ArgumentList "-interaction=nonstopmode", "-file-line-error", "$TexFile.tex" -NoNewWindow -Wait -PassThru
    if ($proc.ExitCode -ne 0) {
        Write-Warn "pdflatex returned exit code $($proc.ExitCode)"
        # Check log for errors
        if (Test-Path "$TexFile.log") {
            $errors = Select-String -Path "$TexFile.log" -Pattern "^!" -Context 0,2
            if ($errors) {
                Write-Err "Errors found in log:"
                $errors | ForEach-Object { Write-Host $_.Line }
            }
        }
    }
    return $proc.ExitCode
}

function Invoke-BibTeX {
    if (-not $bibtex) {
        Write-Warn "bibtex not found, skipping bibliography"
        return 0
    }
    Write-Info "Running bibtex..."
    $proc = Start-Process -FilePath $bibtex -ArgumentList "$TexFile" -NoNewWindow -Wait -PassThru
    if ($proc.ExitCode -ne 0) {
        Write-Warn "bibtex returned exit code $($proc.ExitCode)"
        if (Test-Path "$TexFile.blg") {
            $errors = Select-String -Path "$TexFile.blg" -Pattern "Error|Warning"
            if ($errors) {
                $errors | ForEach-Object { Write-Host $_.Line }
            }
        }
    }
    return $proc.ExitCode
}

# Compilation modes
$startTime = Get-Date

if ($Draft) {
    # Quick draft mode: single pdflatex pass
    Write-Info "=== Draft Mode (single pass) ==="
    $exitCode = Invoke-PdfLatex -Pass "draft"
} else {
    # Full compilation: pdflatex -> bibtex -> pdflatex -> pdflatex
    Write-Info "=== Full Compilation ==="

    # First pass
    $exitCode = Invoke-PdfLatex -Pass "1st pass"
    if ($exitCode -ne 0 -and -not $Final) {
        Write-Err "First pass failed. Use -Final to force continue."
        exit $exitCode
    }

    # BibTeX
    if (Test-Path "references.bib") {
        Invoke-BibTeX
    } else {
        Write-Warn "references.bib not found, skipping bibliography"
    }

    # Second pass (resolve references)
    $exitCode = Invoke-PdfLatex -Pass "2nd pass"

    # Third pass (final)
    $exitCode = Invoke-PdfLatex -Pass "3rd pass"
}

$endTime = Get-Date
$elapsed = ($endTime - $startTime).TotalSeconds

# Check result
if (Test-Path "$TexFile.pdf") {
    $pdfSize = (Get-Item "$TexFile.pdf").Length / 1KB
    Write-Success "=== Compilation Complete ==="
    Write-Success "Output: $TexFile.pdf ($([math]::Round($pdfSize, 1)) KB)"
    Write-Success "Time: $([math]::Round($elapsed, 1)) seconds"

    # Check for warnings
    if (Test-Path "$TexFile.log") {
        $warnings = (Select-String -Path "$TexFile.log" -Pattern "Warning" -AllMatches).Matches.Count
        $underfull = (Select-String -Path "$TexFile.log" -Pattern "Underfull" -AllMatches).Matches.Count
        $overfull = (Select-String -Path "$TexFile.log" -Pattern "Overfull" -AllMatches).Matches.Count
        if ($warnings -gt 0 -or $underfull -gt 0 -or $overfull -gt 0) {
            Write-Warn "Warnings: $warnings, Underfull: $underfull, Overfull: $overfull"
        }
    }

    # Open PDF if requested
    if ($View) {
        Write-Info "Opening PDF..."
        Start-Process "$TexFile.pdf"
    }
} else {
    Write-Err "=== Compilation Failed ==="
    Write-Err "Check $TexFile.log for details"
    exit 1
}
