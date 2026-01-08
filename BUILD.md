# Radia Project - Build Guide

Complete guide for building Radia Python modules on Windows, macOS, and Linux.

---

## Quick Start

### Windows (Recommended: MSVC + Intel MKL)

```powershell
# Build both radia and radia_ngsolve (recommended)
.\BuildMSVC.ps1

# Clean rebuild
.\BuildMSVC.ps1 -Rebuild
```

**Outputs**: `build-msvc/radia.pyd`, `build-msvc/radia_ngsolve.pyd`

**Requirements**:
- Visual Studio 2022 (MSVC compiler)
- Intel oneAPI Base Toolkit (for Intel MKL)

### macOS / Linux

```bash
# Install Intel oneAPI Base Toolkit (recommended)
# Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html

# Source Intel environment
source /opt/intel/oneapi/setvars.sh        # Linux
source /opt/intel/oneapi/setvars.sh        # macOS (Intel)

# Build
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**Outputs**: `build/radia.cpXX-<platform>.so`

---

## Build Scripts (Windows)

### BuildMSVC.ps1 - Primary Build Script (Recommended)

Builds `radia.pyd` and `radia_ngsolve.pyd` using MSVC + Intel MKL.

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-Rebuild` | switch | false | Clean + Configure + Build |
| `-Verbose` | switch | false | Detailed output |

**Examples**:

```powershell
.\BuildMSVC.ps1                  # Standard build
.\BuildMSVC.ps1 -Rebuild         # Clean rebuild
```

**Features**:
- Uses MSVC compiler for NGSolve ABI compatibility
- Links Intel MKL for BLAS/LAPACK operations
- Uses Intel OpenMP (`libiomp5md.dll`) for parallelization
- Automatically copies required DLLs to output directory

---

## System Requirements

### Windows (Primary Platform)

- **OS**: Windows 10/11 (64-bit)
- **Compiler**: Visual Studio 2022 (Community or higher)
  - C++ Desktop Development workload
  - CMake tools component
- **Python**: 3.8+ (64-bit)
- **Intel oneAPI Base Toolkit**: Required for Intel MKL
  - Download from [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
  - Provides: Intel MKL (BLAS/LAPACK/FFT), Intel OpenMP

**Required Intel Components**:
| Component | Purpose |
|-----------|---------|
| Intel MKL | BLAS, LAPACK, DFTI (FFT) operations |
| Intel OpenMP (`libiomp5md.dll`) | Thread parallelization |

### macOS

- **OS**: macOS 10.15+
- **Tools**: Xcode Command Line Tools
- **CMake**: `brew install cmake`
- **Python**: 3.8+
- **Intel oneAPI Base Toolkit**: For Intel MKL and OpenMP
  - Download from [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
  - Note: Apple Silicon (arm64) requires Rosetta 2 for Intel MKL

### Linux (Ubuntu/Debian)

- **OS**: Ubuntu 20.04+ / Debian 11+
- **Compiler**: GCC/G++ 9+
- **CMake**: `sudo apt install cmake`
- **Python**: `sudo apt install python3-dev`
- **Intel oneAPI Base Toolkit**: For Intel MKL and OpenMP
  - Add Intel APT repository and install: `sudo apt install intel-oneapi-mkl-devel`
  - Or download from [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)

### For radia_ngsolve (all platforms)

- **NGSolve**: `pip install ngsolve==6.2.2405`
- **Note**: NGSolve 6.2.2405 is required (later versions have periodic BC regression)

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

**Output**: `build-msvc/radia.pyd` (Windows) or `build/radia.cpXX-<platform>.so` (Unix)

**Features**:
- Complete Radia API for magnetic field computation
- Magnetostatic solver with LU, BiCGSTAB, and HACApK methods
- FastImp Conductor module for eddy current analysis (SIBC)
- ExaFMM-t Fast Multipole Method acceleration
- Python callback for custom background fields
- Cross-platform support (Windows primary)

**Usage**:
```python
import radia as rad

# Create hexahedral magnet (20x20x30 mm, magnetization 1.2 T in z)
# ObjHexahedron auto-generates face topology from 8 vertices
vertices = [[-10,-10,-15], [10,-10,-15], [10,10,-15], [-10,10,-15],
            [-10,-10,15], [10,-10,15], [10,10,15], [-10,10,15]]
magnet = rad.ObjHexahedron(vertices, [0, 0, 1.2])

# Compute field
field = rad.Fld(magnet, 'b', [10, 10, 10])
print(f"B = {field} T")
```

### radia_ngsolve - NGSolve Integration

**Output**: `build-msvc/radia_ngsolve.pyd` (Windows only currently)

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

# Create hexahedral magnet (20x20x30 mm, magnetization 1.2 T in z)
# ObjHexahedron auto-generates face topology from 8 vertices
vertices = [[-10,-10,-15], [10,-10,-15], [10,10,-15], [-10,10,-15],
            [-10,-10,15], [10,-10,15], [10,10,15], [-10,10,15]]
magnet = rad.ObjHexahedron(vertices, [0, 0, 1.2])

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
sys.path.insert(0, r'build-msvc')

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

### Install via pip (Recommended)

```bash
pip install radia
```

Or build and install locally:

```powershell
# Windows - build wheel and install
python -m build
pip install dist/radia-*.whl
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

#### "DLL load failed" (Intel MKL)
**Cause**: Intel MKL DLLs not found
**Solution**: Ensure Intel oneAPI is installed and DLLs are copied:
```powershell
# BuildMSVC.ps1 automatically copies required DLLs
.\BuildMSVC.ps1 -Rebuild
```

#### OpenMP issues (two OpenMP runtimes)
**Symptom**: Random crashes or incorrect parallel execution
**Cause**: Both `libiomp5md.dll` (Intel) and `vcomp140.dll` (MSVC) loaded
**Solution**: Verify only Intel OpenMP is loaded:
```python
import psutil, os
process = psutil.Process(os.getpid())
for dll in process.memory_maps():
    if 'omp' in dll.path.lower():
        print(dll.path)
# Should show: libiomp5md.dll
# Should NOT show: vcomp140.dll
```

#### Build fails with linking errors
```powershell
# Clean rebuild
.\BuildMSVC.ps1 -Rebuild

# Or full clean
Remove-Item -Recurse -Force build-msvc
.\BuildMSVC.ps1
```

### macOS

#### Intel MKL not found
```bash
# Install Intel oneAPI Base Toolkit
# Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html

# Source environment before building
source /opt/intel/oneapi/setvars.sh
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
sudo apt install build-essential cmake python3-dev

# Install Intel MKL
# Option 1: APT repository (Ubuntu/Debian)
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | sudo gpg --dearmor -o /usr/share/keyrings/oneapi-archive-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-mkl-devel

# Option 2: Download installer from Intel website
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

### FastImp Conductor Module (Eddy Currents)

Analyze eddy currents in conductors using Surface Impedance Boundary Condition (SIBC):

```python
import radia as rad

# Create conductor block (center, dimensions, conductivity, panels)
cond = rad.CndRecBlock([0, 0, 0], [0.01, 0.01, 0.001], 5.8e7, 4)

# Set operating frequency
rad.CndSetFrequency(cond, 1e6)  # 1 MHz

# Solve and compute field
rad.CndSolve(cond)
B = rad.CndFld(cond, 'b', [0.02, 0, 0])  # Complex B-field [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]
```

**Available Conductor Types**:
- `CndRecBlock`: Rectangular block conductor
- `CndLoop`: Circular loop conductor
- `CndWire`: Arbitrary wire path
- `CndSpiral`: Spiral coil conductor

### ExaFMM-t Fast Multipole Method

For large conductor problems, enable FMM acceleration:

```python
import radia as rad

# Check FMM status
print(rad.CndFmmGetEnabled())  # Returns True if FMM is enabled

# Configure FMM parameters
rad.CndFmmSetParameters(6, 64, 1.0, True)  # (order, ncrit, theta, verbose)

# FMM is automatically used when enabled
rad.CndSolve(cond)
```

**FMM Parameters**:
| Parameter | Description | Default |
|-----------|-------------|---------|
| `order` | Multipole expansion order | 6 |
| `ncrit` | Max particles per leaf | 64 |
| `theta` | Multipole acceptance criterion | 1.0 |
| `verbose` | Print FMM statistics | False |

### CI/CD Integration

```yaml
# GitHub Actions example
- name: Build Radia
  shell: pwsh
  run: |
    .\BuildMSVC.ps1

- name: Upload artifacts
  uses: actions/upload-artifact@v3
  with:
    name: radia-modules
    path: build-msvc/*.pyd
```

---

## Platform-Specific Notes

### Windows (Primary Platform)

- **Intel MKL** for all BLAS/LAPACK/FFT operations
- **Intel OpenMP** (`libiomp5md.dll`) for parallelization
- **MSVC** compiler for NGSolve ABI compatibility
- BuildMSVC.ps1 automatically copies required DLLs
- Supports both cmd.exe and PowerShell

**Required DLLs** (auto-copied by BuildMSVC.ps1):
| DLL | Purpose |
|-----|---------|
| `mkl_rt.2.dll` | MKL runtime (single dynamic library) |
| `mkl_core.2.dll` | MKL core |
| `mkl_intel_thread.2.dll` | MKL threading |
| `mkl_def.2.dll`, `mkl_avx2.dll` | CPU-specific kernels |
| `libiomp5md.dll` | Intel OpenMP runtime |

### macOS

- Intel (x86_64) fully supported with Intel MKL
- Apple Silicon (arm64) requires Rosetta 2 for Intel MKL
- Uses Intel OpenMP (`libiomp5`) for parallelization
- Source `setvars.sh` before building: `source /opt/intel/oneapi/setvars.sh`

### Linux

- Tested on Ubuntu 20.04+ and Debian 11+
- Uses Intel MKL for BLAS/LAPACK/FFT operations
- Uses Intel OpenMP (`libiomp5`) for parallelization
- Source `setvars.sh` before building: `source /opt/intel/oneapi/setvars.sh`
- May need to set `LD_LIBRARY_PATH` for Intel libraries

---

## Deprecated Scripts

The following old scripts have been replaced:

| Old Script | Replacement |
|------------|-------------|
| `Build.ps1` | `.\BuildMSVC.ps1` |
| `Build_NGSolve.ps1` | `.\BuildMSVC.ps1` |
| `build_radia_ngsolve.bat` | `.\BuildMSVC.ps1` |
| `build_radia_ngsolve_full.bat` | `.\BuildMSVC.ps1 -Rebuild` |

Old scripts are removed from the repository.

---

## Development Workflow

```powershell
# Initial setup
.\BuildMSVC.ps1 -Rebuild

# During development (incremental builds)
.\BuildMSVC.ps1

# Test import
python -c "import sys; sys.path.insert(0, 'build-msvc'); import radia; print(radia.UtiVer())"

# Clean everything
Remove-Item -Recurse -Force build-msvc
.\BuildMSVC.ps1
```

---

## Getting Help

- **Documentation**: See `docs/` folder
- **Examples**: See `examples/` folder
- **Issues**: GitHub Issues
- **API Reference**: `docs/API.md`

---

**Last Updated**: 2026-01-08
**Build System**: CMake 3.21+ with MSVC / GCC / Clang
**Primary Platform**: Windows with Intel MKL + Intel OpenMP
**Supported Platforms**: Windows (primary), macOS, Linux
**Supported Python**: 3.8, 3.9, 3.10, 3.11, 3.12
**C++ Standard**: C++17 (required for ExaFMM-t inline variables)
