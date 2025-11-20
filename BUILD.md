# Radia Project Build Guide

## Quick Start

### Build Everything (Recommended)
```powershell
.\build_all.ps1
```

### Build Only Radia Module
```powershell
.\Build.ps1
```

### Build Only rad_ngsolve Module
```powershell
# Requires Visual Studio environment
powershell -Command "& { Import-Module 'C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\Microsoft.VisualStudio.DevShell.dll'; Enter-VsDevShell -VsInstallPath 'C:\Program Files\Microsoft Visual Studio\2022\Community' -SkipAutomaticLocation; cmake --build build --config Release --target rad_ngsolve }"
```

---

## Build Scripts

### build_all.ps1 - Build Both Modules

### Parameters

| Parameter | Values | Default | Description |
|-----------|--------|---------|-------------|
| `-Target` | `all`, `radia`, `rad_ngsolve` | `all` | What to build |
| `-BuildType` | `Release`, `Debug`, `RelWithDebInfo` | `Release` | Build configuration |
| `-Clean` | (switch) | false | Clean build directory before build |
| `-Rebuild` | (switch) | false | Clean + Configure + Build |
| `-Test` | (switch) | false | Run module import tests after build |
| `-Verbose` | (switch) | false | Show detailed build output |
| `-Install` | (switch) | false | Install to Python site-packages |

### Examples

**Full clean rebuild with tests**:
```powershell
.\Build.ps1 -Rebuild -Test
```

**Build only radia in Debug mode**:
```powershell
.\Build.ps1 -Target radia -BuildType Debug
```

**Quick incremental build**:
```powershell
.\Build.ps1
```

**Build and install rad_ngsolve**:
```powershell
.\Build.ps1 -Target rad_ngsolve -Install
```

---

## Build Targets

### `radia` - Radia Python Module
- **Output**: `build/Release/radia.pyd`
- **Usage**: `import radia as rad`
- **Requirements**: Python 3.12 (64-bit)

### `rad_ngsolve` - NGSolve Integration
- **Output**: `build/Release/rad_ngsolve.pyd`
- **Usage**: `import rad_ngsolve`
- **Requirements**:
  - Python 3.12 (64-bit)
  - NGSolve installed

### `all` - Build Both
- Builds both `radia.pyd` and `rad_ngsolve.pyd`
- Default target

---

## System Requirements

### Required
- **Windows 10/11** (64-bit)
- **Visual Studio 2022** (Community, Professional, or Enterprise)
  - C++ Desktop Development workload
  - CMake tools component
- **Python 3.12** (64-bit)

### For rad_ngsolve
- **NGSolve** (with Python 3.12 support)
  ```bash
  pip install ngsolve
  ```

---

## Troubleshooting

### "CMake not found"
- Install CMake component in Visual Studio Installer
- Or: Add CMake to system PATH

### "Python 3.12 not found"
- Download from [python.org](https://www.python.org/downloads/)
- Ensure "Add to PATH" is checked during installation
- **Important**: Install 64-bit version

### "Module import failed"
- Check Python version: `python --version` (should be 3.12.x)
- Check architecture: `python -c "import struct; print(struct.calcsize('P')*8)"` (should be 64)
- For rad_ngsolve: Verify NGSolve is installed: `python -c "import ngsolve; print(ngsolve.__version__)"`

### Build fails with linking errors
- Clean rebuild: `.\Build.ps1 -Rebuild`
- Check Visual Studio 2022 is installed correctly
- Verify Python libraries are 64-bit

---

## Old Build Scripts (Deprecated)

The following scripts have been replaced by `Build.ps1`:
- ~~`build_ngsolve.ps1`~~ → Use `Build.ps1 -Target rad_ngsolve -Rebuild`
- ~~`build_rad_ngsolve.bat`~~ → Use `Build.ps1 -Target rad_ngsolve`
- ~~`build_rad_ngsolve_full.bat`~~ → Use `Build.ps1 -Target rad_ngsolve -Rebuild`

Old scripts are preserved as `*.old` files for reference.

---

## Advanced Usage

### Parallel Builds
```powershell
# Build multiple configurations in parallel
Start-Job {.\Build.ps1 -BuildType Release -Target radia}
Start-Job {.\Build.ps1 -BuildType Debug -Target radia}
Get-Job | Wait-Job | Receive-Job
```

### Custom Build Directory
```powershell
# Not directly supported, edit Build.ps1 line 24:
# $BuildDir = Join-Path $ScriptDir "build"
```

### CI/CD Integration
```yaml
# Example GitHub Actions
- name: Build Radia
  run: powershell -ExecutionPolicy Bypass -File Build.ps1 -Rebuild -Test
```

---

**Last Updated**: 2025-11-20
**Build System**: CMake 3.31 + Visual Studio 2022 + Python 3.12
