# Radia Project - Build Guide

Complete guide for building Radia Python modules on Windows, macOS, and Linux.

---

## Quick Start

### Windows

```powershell
# Build both radia and radia_ngsolve (recommended)
.\Build.ps1

# Build with tests
.\Build.ps1 -Test

# Clean rebuild
.\Build.ps1 -Clean
```

**Outputs**: `build/Release/radia.pyd`, `build/Release/radia_ngsolve.pyd`

### macOS / Linux

```bash
# Install dependencies
brew install cmake fftw python3          # macOS
sudo apt install cmake libfftw3-dev python3-dev  # Linux

# Build
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**Outputs**: `build/radia.cpXX-<platform>.so`

---

## Build Scripts (Windows)

### Build.ps1 - Build Both Modules

Builds `radia.pyd` and `radia_ngsolve.pyd` in one command.

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-BuildType` | `Release`, `Debug` | `Release` | Build configuration |
| `-Clean` | switch | false | Clean before build |
| `-Test` | switch | false | Run import tests |

**Examples**:

```powershell
.\Build.ps1                  # Standard build
.\Build.ps1 -Clean -Test     # Clean rebuild with tests
.\Build.ps1 -BuildType Debug # Debug build
```

### Build.ps1 - Build Radia Module Only

Builds only `radia.pyd` (core Radia module).

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-BuildType` | `Release`, `Debug`, `RelWithDebInfo` | `Release` | Build configuration |
| `-Clean` | switch | false | Clean before build |
| `-Rebuild` | switch | false | Clean + Configure + Build |
| `-Verbose` | switch | false | Detailed output |
| `-Install` | switch | false | Install to site-packages |

**Examples**:

```powershell
.\Build.ps1                      # Standard build
.\Build.ps1 -Rebuild             # Clean rebuild
.\Build.ps1 -BuildType Debug     # Debug build
.\Build.ps1 -Install             # Build and install
```

---

## System Requirements

### Windows

- **OS**: Windows 10/11 (64-bit)
- **Compiler**: Visual Studio 2022 (Community or higher)
  - C++ Desktop Development workload
  - CMake tools component
- **Python**: 3.8+ (64-bit)
- **FFTW**: Included in `src/ext/fftw/fftw64_f.lib`

### macOS

- **OS**: macOS 10.15+
- **Tools**: Xcode Command Line Tools
- **CMake**: `brew install cmake`
- **Python**: 3.8+
- **FFTW**: `brew install fftw`

### Linux (Ubuntu/Debian)

- **OS**: Ubuntu 20.04+ / Debian 11+
- **Compiler**: GCC/G++ 9+
- **CMake**: `sudo apt install cmake`
- **Python**: `sudo apt install python3-dev`
- **FFTW**: `sudo apt install libfftw3-dev`

### For radia_ngsolve (all platforms)

- **NGSolve**: `pip install ngsolve`

---

## Python Version Compatibility

### Important Constraint

**Each Python version requires a separate .pyd/.so file** due to ABI incompatibility:

- Python 3.8, 3.9, 3.10, 3.11, 3.12 each have different binary interfaces
- Internal data structures vary between versions
- Function signatures may change across versions

### Solution

Build for each Python version you need to support:

```powershell
# Build for current Python
.\Build.ps1

# Build for multiple versions
# Switch Python and rebuild for each version
```

### File Naming Convention (PEP 3149)

```
radia.cp<version>-<platform>.<ext>

Examples:
  radia.cp312-win_amd64.pyd       # Python 3.12, Windows 64-bit
  radia.cp311-win_amd64.pyd       # Python 3.11, Windows 64-bit
  radia.cp312-darwin.so           # Python 3.12, macOS
  radia.cp312-linux_x86_64.so     # Python 3.12, Linux
```

---

## Build Targets

### radia - Core Python Module

**Output**: `build/Release/radia.pyd` (Windows) or `radia.cpXX-<platform>.so` (Unix)

**Features**:
- Complete Radia API
- Magnetic field computation
- Python callback for custom fields (`ObjBckgCF`)
- Cross-platform support

**Usage**:
```python
import radia as rad

# Create magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])

# Compute field
field = rad.Fld(magnet, 'b', [10, 10, 10])
print(f"B = {field} T")
```

### radia_ngsolve - NGSolve Integration

**Output**: `build/Release/radia_ngsolve.pyd` (Windows only currently)

**Features**:
- NGSolve CoefficientFunction interface
- Automatic unit conversion (m ↔ mm)
- Support for B, H, A, M fields
- Coordinate transformations

**Usage**:
```python
import ngsolve  # MUST import first
import radia_ngsolve
import radia as rad

# Create Radia magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])

# Create NGSolve CoefficientFunction
B_cf = radia_ngsolve.RadiaField(magnet, 'b')  # Flux density
H_cf = radia_ngsolve.RadiaField(magnet, 'h')  # Magnetic field
A_cf = radia_ngsolve.RadiaField(magnet, 'a')  # Vector potential

# Use in NGSolve mesh
from netgen.occ import *
box = Box((0, 0, 0), (0.1, 0.1, 0.1))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.01))

# Project to GridFunction
from ngsolve import *
fes = HDiv(mesh, order=2)
B_gf = GridFunction(fes)
B_gf.Set(B_cf)
```

**Field Types**:
- `'b'`: Magnetic flux density (Tesla)
- `'h'`: Magnetic field (A/m)
- `'a'`: Vector potential (T·m)
- `'m'`: Magnetization (A/m)

---

## Usage After Build

### Windows

```python
import sys
sys.path.insert(0, r'build\Release')

import radia as rad
print(rad.UtiVer())

# If radia_ngsolve was built
import ngsolve
import radia_ngsolve
```

### Unix

```python
import sys
sys.path.insert(0, 'build')

import radia as rad
print(rad.UtiVer())
```

### Install to Site-Packages

```powershell
# Windows
.\Build.ps1 -Install

# Unix
cd build
sudo make install
```

---

## Troubleshooting

### Windows

#### "CMake not found"
- Open Visual Studio Installer
- Modify Visual Studio 2022
- Check "C++ CMake tools for Windows"
- Install

#### "Python 3.x not found"
- Download from [python.org](https://www.python.org/downloads/)
- **Important**: Install 64-bit version
- Check "Add Python to PATH" during installation

#### "Module import failed"
```bash
# Check Python version
python --version  # Should match build version

# Check architecture
python -c "import struct; print(struct.calcsize('P')*8)"  # Should be 64

# For radia_ngsolve
python -c "import ngsolve; print(ngsolve.__version__)"
```

#### "DLL load failed" (radia_ngsolve)
**Solution**: Import ngsolve before radia_ngsolve:
```python
import ngsolve  # Load NGSolve dependencies first
import radia_ngsolve  # Now this will work
```

#### Build fails with linking errors
```powershell
# Clean rebuild
.\Build.ps1 -Clean

# Or full clean
Remove-Item -Recurse -Force build
.\Build.ps1
```

### macOS

#### "FFTW not found"
```bash
brew install fftw
```

#### Python.h not found
```bash
# macOS comes with Python, but install full version:
brew install python@3.12
```

### Linux

#### Missing dependencies
```bash
sudo apt update
sudo apt install build-essential cmake python3-dev libfftw3-dev
```

#### Permission denied during install
```bash
sudo make install
```

---

## Advanced Topics

### Custom Background Fields (ObjBckgCF)

Define magnetic fields via Python callback:

```python
import radia as rad

def gradient_field(pos):
    """
    pos: [x, y, z] in millimeters
    returns: {'B': [Bx, By, Bz] in Tesla}
    """
    x, y, z = pos
    return {
        'B': [0.01 * x/1000, 0.01 * y/1000, 0.01 * z/1000]
    }

# Create background field
bg_field = rad.ObjBckgCF(gradient_field)

# Evaluate
B = rad.Fld(bg_field, 'b', [10, 20, 30])
print(f"B = {B} T")
```

**Limitations**:
- Binary serialization not supported
- Vector potential (A) computation not implemented
- Infinite integral uses simple trapezoidal rule

### Coordinate Transformations (radia_ngsolve)

Transform between global and local coordinate systems:

```python
B_cf = radia_ngsolve.RadiaField(
    magnet, 'b',
    origin=[0.05, 0.05, 0.05],      # Translation (meters)
    u_axis=[1, 0, 0],                # Local x-axis
    v_axis=[0, 1, 0],                # Local y-axis
    w_axis=[0, 0, 1]                 # Local z-axis
)
```

### Parallel Builds

```powershell
# Build multiple configurations in parallel
Start-Job {.\Build.ps1 -BuildType Release}
Start-Job {.\Build.ps1 -BuildType Debug}
Get-Job | Wait-Job | Receive-Job
```

### CI/CD Integration

```yaml
# GitHub Actions example
- name: Build Radia
  shell: pwsh
  run: |
    .\Build.ps1 -Test

- name: Upload artifacts
  uses: actions/upload-artifact@v3
  with:
    name: radia-modules
    path: build/Release/*.pyd
```

---

## Platform-Specific Notes

### Windows

- FFTW library is **bundled** in `src/ext/fftw/fftw64_f.lib`
- No external FFTW installation needed
- Supports both cmd.exe and PowerShell

### macOS

- Both Intel (x86_64) and Apple Silicon (arm64) supported
- FFTW must be installed via Homebrew
- Universal binaries not yet supported (build for specific architecture)

### Linux

- Tested on Ubuntu 20.04+ and Debian 11+
- FFTW available via package manager
- May need to adjust `LD_LIBRARY_PATH` if FFTW is in non-standard location

---

## Deprecated Scripts

The following old scripts have been replaced:

| Old Script | Replacement |
|------------|-------------|
| `Build_NGSolve.ps1` | `.\Build.ps1` |
| `build_radia_ngsolve.bat` | `.\Build.ps1` |
| `build_radia_ngsolve_full.bat` | `.\Build.ps1 -Clean` |

Old scripts are removed from the repository.

---

## Development Workflow

```powershell
# Initial setup
.\Build.ps1 -Clean -Test

# During development (incremental builds)
.\Build.ps1

# Before committing
.\Build.ps1 -Test

# Clean everything
Remove-Item -Recurse -Force build
```

---

## Getting Help

- **Documentation**: See `docs/` folder
- **Examples**: See `examples/` folder
- **Issues**: GitHub Issues
- **API Reference**: `docs/API.md`

---

**Last Updated**: 2025-11-20
**Build System**: CMake 3.21+ with Visual Studio 2022 / GCC / Clang
**Supported Platforms**: Windows, macOS, Linux
**Supported Python**: 3.8, 3.9, 3.10, 3.11, 3.12
