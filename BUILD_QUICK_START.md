# Radia - Quick Build Guide

## Prerequisites
- Windows 10/11 (64-bit)
- Visual Studio 2022 (Community or higher)
- Python 3.12 (64-bit)
- For rad_ngsolve: NGSolve (`pip install ngsolve`)

## Build Commands

### Build Everything (Recommended)
```powershell
.\Build.ps1
```
Builds both `radia.pyd` and `rad_ngsolve.pyd`.

### Build with Tests
```powershell
.\Build.ps1 -Test
```

### Clean Rebuild
```powershell
.\Build.ps1 -Clean
```

## Build Scripts Organization

| Script | Purpose | Output |
|--------|---------|--------|
| `Build.ps1` | **Build both modules** | `radia.pyd` + `rad_ngsolve.pyd` |

### Deprecated Scripts (Moved to *.old)
- `build_ngsolve.ps1` → Use `Build.ps1`
- `build_rad_ngsolve.bat` → Use `Build.ps1`
- `build_rad_ngsolve_full.bat` → Use `Build.ps1 -Clean`

## Output Location
```
build/Release/radia.pyd
build/Release/rad_ngsolve.pyd
```

## Usage After Build
```python
import sys
sys.path.insert(0, r'build\Release')

# Use Radia
import radia as rad
print(rad.UtiVer())

# Use rad_ngsolve (requires NGSolve)
import ngsolve
import rad_ngsolve
```

## Troubleshooting

### "Python 3.12 not found"
Install from [python.org](https://www.python.org/downloads/) (64-bit version)

### "CMake not found"
Install CMake component via Visual Studio Installer

### "NGSolve not found" (for rad_ngsolve)
```bash
pip install ngsolve
```

### Build fails
```powershell
.\Build.ps1 -Clean  # Clean rebuild
```

## Development Workflow

```powershell
# Initial build
.\Build.ps1 -Clean -Test

# After code changes (incremental)
.\Build.ps1

# Before commit
.\Build.ps1 -Test
```

---
**For detailed documentation**, see `BUILD.md`
**For full documentation**, see `docs/`
