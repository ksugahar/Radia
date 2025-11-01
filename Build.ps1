#==============================================================================
# Radia Python 3.12 Module (.pyd) Build Script for Visual Studio 2022
# This script builds ONLY the Python 3.12 extension module (64-bit)
#==============================================================================

param(
    [ValidateSet("Release", "Debug", "RelWithDebInfo")]
    [string]$BuildType = "Release",

    [switch]$Clean,
    [switch]$Rebuild,
    [switch]$Install,
    [switch]$Verbose
)

#==============================================================================
# Configuration
#==============================================================================

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$BuildDir = Join-Path $ScriptDir "build"
$LogFile = Join-Path $ScriptDir "build.log"
$RequiredPythonVersion = "3.12"

#==============================================================================
# Functions
#==============================================================================

function Write-Header {
    param([string]$Message)
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host "  $Message" -ForegroundColor Cyan
    Write-Host "========================================" -ForegroundColor Cyan
}

function Write-Info {
    param([string]$Message)
    Write-Host "[INFO]  $Message" -ForegroundColor Green
}

function Write-Warning {
    param([string]$Message)
    Write-Host "[WARN]  $Message" -ForegroundColor Yellow
}

function Write-ErrorMsg {
    param([string]$Message)
    Write-Host "[ERROR] $Message" -ForegroundColor Red
}

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host ">>> $Message" -ForegroundColor Cyan
}

function Find-CMake {
    Write-Step "Locating CMake from Visual Studio 2022..."

    $VSBasePath = "C:\Program Files\Microsoft Visual Studio\2022"
    $VSEditions = @("Enterprise", "Professional", "Community")

    foreach ($edition in $VSEditions) {
        $cmakePath = Join-Path $VSBasePath "$edition\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
        if (Test-Path $cmakePath) {
            Write-Info "Found CMake: $cmakePath"
            Write-Info "Edition: Visual Studio 2022 $edition"
            return $cmakePath
        }
    }

    # Try system PATH as fallback
    $systemCMake = Get-Command cmake -ErrorAction SilentlyContinue
    if ($systemCMake) {
        Write-Warning "Using system CMake from PATH: $($systemCMake.Source)"
        return $systemCMake.Source
    }

    throw "CMake not found! Please install CMake component in Visual Studio 2022."
}

function Find-Python312 {
    Write-Step "Locating Python 3.12 (64-bit)..."

    # Try to find Python in PATH
    $pythonCmd = Get-Command python -ErrorAction SilentlyContinue
    if ($pythonCmd) {
        $pythonExe = $pythonCmd.Source

        # Check version
        $pythonVersionOutput = & $pythonExe --version 2>&1
        if ($pythonVersionOutput -match "Python (\d+\.\d+)") {
            $pythonVersion = $matches[1]
        } else {
            throw "Could not determine Python version"
        }

        # Check if it's 64-bit
        $pythonBits = & $pythonExe -c "import struct; print(struct.calcsize('P') * 8)" 2>&1

        # Verify version is 3.12
        if ($pythonVersion -eq $RequiredPythonVersion) {
            if ($pythonBits -eq "64") {
                Write-Info "Found Python: $pythonExe"
                Write-Info "Version: Python $pythonVersion"
                Write-Info "Architecture: 64-bit"
                return $pythonExe
            } else {
                throw "Python 3.12 found but it is $pythonBits-bit. 64-bit Python 3.12 is required."
            }
        } else {
            throw "Python $pythonVersion found, but Python $RequiredPythonVersion is required."
        }
    }

    throw "Python 3.12 (64-bit) not found in PATH. Please install Python 3.12 64-bit."
}

function Invoke-Build {
    param(
        [string]$CMakePath,
        [string]$BuildDir,
        [string]$BuildType
    )

    # Start timing
    $buildStartTime = Get-Date

    # Prepare CMake configuration options
    $cmakeArgs = @(
        "..",
        "-G", "Visual Studio 17 2022",
        "-A", "x64"
    )

    # Configure
    Write-Step "Configuring CMake project (Python 3.12 .pyd - 64-bit)..."
    Write-Info "Configuration: $BuildType"
    Write-Info "Generator: Visual Studio 17 2022"
    Write-Info "Architecture: x64 (64-bit)"
    Write-Info "Target: radia.pyd (Python 3.12 extension module)"

    if ($Verbose) {
        Write-Host "CMake command: $CMakePath $($cmakeArgs -join ' ')" -ForegroundColor Gray
    }

    Push-Location $BuildDir
    try {
        & $CMakePath $cmakeArgs 2>&1 | Tee-Object -FilePath $LogFile

        if ($LASTEXITCODE -ne 0) {
            throw "CMake configuration failed with exit code $LASTEXITCODE"
        }

        Write-Info "Configuration completed successfully"

        # Build
        Write-Step "Building Python 3.12 module ($BuildType)..."

        $buildArgs = @(
            "--build", ".",
            "--config", $BuildType,
            "--parallel"
        )

        if ($Verbose) {
            $buildArgs += "--verbose"
        }

        & $CMakePath $buildArgs 2>&1 | Tee-Object -FilePath $LogFile -Append

        if ($LASTEXITCODE -ne 0) {
            throw "Build failed with exit code $LASTEXITCODE"
        }

        # Calculate build time
        $buildEndTime = Get-Date
        $buildDuration = $buildEndTime - $buildStartTime

        # Success
        Write-Header "Build Completed Successfully"
        Write-Info "Build time: $($buildDuration.ToString('mm\:ss'))"
        Write-Info "Configuration: $BuildType"

        # Find and display output files
        $pydFiles = Get-ChildItem -Path $BuildDir -Recurse -Filter "radia.pyd" -ErrorAction SilentlyContinue

        if ($pydFiles) {
            Write-Info "Output files:"
            foreach ($file in $pydFiles) {
                $sizeMB = [math]::Round($file.Length / 1MB, 2)
                Write-Host "  - $($file.FullName) ($sizeMB MB)" -ForegroundColor White

                # Save the main output path
                $script:OutputPydPath = $file.FullName
            }
        } else {
            Write-Warning "radia.pyd not found in build directory"
        }

        Write-Info "Log file: $LogFile"

        # Install option
        if ($Install) {
            Write-Step "Installing Python module..."

            $installArgs = @(
                "--install", ".",
                "--config", $BuildType
            )

            & $CMakePath $installArgs 2>&1 | Tee-Object -FilePath $LogFile -Append

            if ($LASTEXITCODE -eq 0) {
                Write-Info "Installation completed successfully"
            } else {
                Write-Warning "Installation failed with exit code $LASTEXITCODE"
            }
        }

    } finally {
        Pop-Location
    }
}

function Test-PythonModule {
    Write-Step "Testing Python 3.12 module..."

    if (-not $script:OutputPydPath -or -not (Test-Path $script:OutputPydPath)) {
        Write-Warning "Cannot test: radia.pyd not found"
        return
    }

    # Get the directory containing the .pyd file
    $pydDir = Split-Path -Parent $script:OutputPydPath

    # Try to import the module
    $testScript = @"
import sys
sys.path.insert(0, r'$pydDir')
try:
    import radia
    print('SUCCESS: radia module imported successfully')
    print(f'Module location: {radia.__file__}')
    print(f'Module attributes: {dir(radia)}')
except Exception as e:
    print(f'ERROR: Failed to import radia module: {e}')
    sys.exit(1)
"@

    $pythonExe = (Get-Command python -ErrorAction SilentlyContinue).Source
    if ($pythonExe) {
        Write-Info "Testing module import..."
        $testScript | & $pythonExe -

        if ($LASTEXITCODE -eq 0) {
            Write-Info "Module test passed"
        } else {
            Write-Warning "Module test failed"
        }
    } else {
        Write-Warning "Python not found in PATH, skipping module test"
    }
}

#==============================================================================
# Main Script
#==============================================================================

try {
    # Print header
    Write-Header "Radia Python 3.12 Module (.pyd) Build Script"
    Write-Info "Build Type: $BuildType"
    Write-Info "Script Directory: $ScriptDir"
    Write-Info "Build Directory: $BuildDir"
    Write-Info "Target: Python 3.12 extension module (radia.pyd) - 64-bit only"

    # Find CMake
    $CMakePath = Find-CMake

    # Find Python 3.12 (required)
    $PythonExe = Find-Python312

    # Verify architecture
    if ([Environment]::Is64BitOperatingSystem) {
        Write-Info "Operating System: 64-bit Windows"
    } else {
        throw "This script requires 64-bit Windows"
    }

    if ([Environment]::Is64BitProcess) {
        Write-Info "PowerShell Process: 64-bit"
    } else {
        Write-Warning "PowerShell Process: 32-bit (running on 64-bit OS)"
    }

    # Clean or rebuild
    if ($Clean -or $Rebuild) {
        if (Test-Path $BuildDir) {
            Write-Step "Cleaning build directory..."
            Remove-Item -Recurse -Force $BuildDir
            Write-Info "Build directory cleaned"
        }
    }

    # Create build directory
    if (-not (Test-Path $BuildDir)) {
        Write-Step "Creating build directory..."
        New-Item -ItemType Directory -Path $BuildDir | Out-Null
        Write-Info "Build directory created"
    }

    # Execute build
    Invoke-Build -CMakePath $CMakePath -BuildDir $BuildDir -BuildType $BuildType

    # Test the module
    Test-PythonModule

    # Exit success
    Write-Host ""
    exit 0

} catch {
    # Error handling
    Write-Host ""
    Write-Header "Build Failed"
    Write-ErrorMsg $_.Exception.Message
    Write-ErrorMsg "Stack Trace:"
    Write-Host $_.ScriptStackTrace -ForegroundColor Red

    if (Test-Path $LogFile) {
        Write-Warning "Check log file for details: $LogFile"
    }

    Write-Host ""
    exit 1
}
