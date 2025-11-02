# PyPI Distribution Guide for Radia

This guide explains how to build and distribute the Radia package to PyPI.

## Prerequisites

### 1. Build Tools

Install the required build tools:

```bash
pip install build twine
```

### 2. Build the Extension Modules

Before creating the distribution, you must build the C++ extension modules:

#### Core Module (Required)

```powershell
# Build radia.pyd
.\Build.ps1
```

This creates `dist\radia.pyd`

#### NGSolve Module (Optional)

```powershell
# Build rad_ngsolve.pyd
.\build_ngsolve.ps1
```

This creates `build\Release\rad_ngsolve.pyd`

## Building Distribution Packages

### Using the Build Script

Run the automated build script:

```cmd
build_pypi_package.cmd
```

This script will:
1. Check for required `.pyd` files
2. Clean previous builds
3. Build source distribution (`.tar.gz`)
4. Build wheel distribution (`.whl`)
5. Show created files in `dist/`

### Manual Build

If you prefer to build manually:

```bash
# Build source distribution
python -m build --sdist

# Build wheel distribution
python -m build --wheel
```

## Package Structure

The distribution includes:

```
radia/
├── src/python/
│   ├── radia.pyd                    # Core C++ extension
│   ├── rad_ngsolve.pyd              # NGSolve integration (optional)
│   ├── radia_pyvista_viewer.py      # PyVista viewer utilities
│   ├── radia_vtk_export.py          # VTK export utilities
│   └── nastran_reader.py            # Nastran mesh reader
├── examples/                        # Example scripts and data
├── tests/                           # Test suite
├── docs/                            # Documentation
├── LICENSE                          # LGPL-2.1 license
├── COPYRIGHT.txt                    # Original RADIA copyright
└── README.md                        # Main documentation
```

## Verifying the Package

Before uploading, verify the package:

```bash
# Check package metadata and files
python -m twine check dist/*
```

## Uploading to PyPI

### Test PyPI (Recommended First)

Upload to Test PyPI to verify everything works:

```bash
# Using the upload script
pypi_upload.cmd
# Select option 1 (Test PyPI)

# Or manually
python -m twine upload --repository testpypi dist/*
```

Test the installation:

```bash
pip install --index-url https://test.pypi.org/simple/ radia
```

### Production PyPI

Once verified on Test PyPI, upload to production:

```bash
# Using the upload script
pypi_upload.cmd
# Select option 2 (PyPI)

# Or manually
python -m twine upload dist/*
```

## PyPI Credentials

### Using .pypirc

Create `~/.pypirc` with your credentials:

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-your-api-token-here

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-your-test-api-token-here
```

### Using Environment Variables

Alternatively, set environment variables:

```bash
# Windows
set TWINE_USERNAME=__token__
set TWINE_PASSWORD=pypi-your-api-token-here

# Or use interactive prompt (recommended for security)
python -m twine upload dist/*
```

## Package Metadata

The package metadata is defined in `pyproject.toml`:

- **Name**: radia
- **Version**: 4.32
- **License**: LGPL-2.1
- **Python**: >=3.12
- **Platform**: Windows (currently)

## Dependencies

### Required
- numpy>=1.20

### Optional (extras)
- **viz**: pyvista>=0.40, matplotlib>=3.5
- **test**: pytest>=7.0, pytest-cov>=4.0
- **dev**: All of the above plus build tools

Install with extras:

```bash
pip install radia[viz]        # Visualization support
pip install radia[test]       # Testing tools
pip install radia[dev]        # Development tools
```

## Important Notes

### Binary Distribution

This package includes pre-compiled binary extensions (`.pyd` files):
- Platform-specific: Windows only (currently)
- Python version-specific: Python 3.12
- Cannot be used on other platforms/versions

### Source Distribution

The source distribution (`.tar.gz`) includes:
- All source code
- Build scripts
- Documentation
- Examples

**Note**: Users cannot build from source without:
- Visual Studio 2022
- CMake
- OpenMP support
- (Optional) NGSolve installation

### Version Updates

To release a new version:

1. Update version in `pyproject.toml`
2. Update version in `setup.py`
3. Rebuild extensions
4. Run `build_pypi_package.cmd`
5. Upload to PyPI

## Troubleshooting

### "radia.pyd not found"

Run `Build.ps1` to build the core module first.

### "Package check failed"

Run `python -m twine check dist/*` to see specific errors.

### "Upload failed"

Common causes:
- Version already exists on PyPI (cannot overwrite)
- Invalid credentials
- Network issues

Solutions:
- Increment version number
- Check `.pypirc` or credentials
- Try again or use VPN

### Import Errors After Installation

Make sure you're using:
- Windows platform
- Python 3.12
- Correct architecture (x64)

## Additional Resources

- PyPA Packaging Guide: https://packaging.python.org/
- Twine Documentation: https://twine.readthedocs.io/
- setuptools Documentation: https://setuptools.pypa.io/

## License

This distribution includes:
- LGPL-2.1 licensed code (matching NGSolve)
- Original RADIA copyright from ESRF (1997-2018)

See `LICENSE` and `COPYRIGHT.txt` for full text.
