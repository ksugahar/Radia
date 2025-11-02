# NGSolve Integration - DLL Dependency Issues

## Current Status

`rad_ngsolve.pyd` was built successfully, but the following error occurs when importing:

```
ImportError: DLL load failed while importing rad_ngsolve: The specified module could not be found.
```

## Dependency Check Results

### Direct Dependencies of rad_ngsolve.pyd

```
- python312.dll     [OK] - Python runtime
- libngsolve.dll    [Copied] - NGSolve core library
- KERNEL32.dll      [OK] - Windows API
- VCOMP140.DLL      [Copied] - OpenMP runtime
```

### Dependencies of libngsolve.dll (Estimated)

```
- ngcore.dll        [Copied]
- nglib.dll         [Copied]
- mkl_rt.lib        [Unknown] - Intel MKL
- Other MKL DLLs    [Unknown]
```

## Attempted Solutions

### ✗ Failed Methods

1. **Copying DLLs**: Copied necessary DLLs from `netgen` folder
   - Result: Same error persists

2. **PATH Environment Variable**: Added netgen directory to PATH
   - Result: Same error persists

3. **Copy to site-packages**: Copied rad_ngsolve.pyd directly to site-packages
   - Result: Same error persists

## Root Cause Analysis

NGSolve likely depends on Intel MKL (Math Kernel Library), and these DLLs cannot be found.

Contents of NGSolveConfig.cmake:
```cmake
set(NGSOLVE_MKL_LIBRARIES "C:/gitlabci/tools/builds/3zsqG5ns/0/ngsolve/venv_ngs/Library/lib/mkl_rt.lib")
```

This path is **specific to the build environment** and does not exist in the actual installation environment.

## Recommended Solutions

### Option 1: Use Conda (Most Reliable)

Using a Conda environment automatically resolves all dependencies:

```bash
# Create NGSolve environment
conda create -n ngsolve_env python=3.12
conda activate ngsolve_env
conda install -c ngsolve ngsolve

# Build Radia (within this environment)
cd S:\radia\01_GitHub
.\build_ngsolve.ps1

# Test
python -c "import rad_ngsolve; print('OK')"
```

### Option 2: Standalone Build (MKL-free Version)

Avoid linking to NGSolve, using headers only:

1. Modify `rad_ngsolve.cpp` to remove NGSolve library linking
2. Obtain only necessary definitions from headers
3. Import required functionality from ngsolve module at runtime

### Option 3: Install MKL

Install Intel oneMKL on the system:

```powershell
# Download and install Intel oneMKL (free)
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

# After installation, add to PATH environment variable
$env:PATH += ";C:\Program Files (x86)\Intel\oneAPI\mkl\latest\redist\intel64"
```

### Option 4: Dynamic Loading Approach

Modify implementation at Python level:

```python
# rad_ngsolve_wrapper.py
import ctypes
import os
import radia as rad

# Dynamically set DLL path
ngsolve_path = r"C:\Program Files\Python312\Lib\site-packages\netgen"
os.add_dll_directory(ngsolve_path)

# Load rad_ngsolve.pyd
import rad_ngsolve as _rad_ngsolve

# Wrapper function
def RadBfield(radia_obj):
	return _rad_ngsolve.RadBfield(radia_obj)
```

## Temporary Workaround

The most practical solution at this time is to change the NGSolve integration to a **pure Python implementation**:

```python
# rad_ngsolve_py.py - Pure Python implementation
import radia as rad
from ngsolve import CoefficientFunction

class RadiaBFieldPy(CoefficientFunction):
	def __init__(self, radia_obj):
	    super().__init__(3)
	    self.radia_obj = radia_obj

	def Evaluate(self, mip, result):
	    pnt = mip.GetPoint()
	    B = rad.Fld(self.radia_obj, 'b', [pnt[0], pnt[1], pnt[2]])
	    result[0] = B[0]
	    result[1] = B[1]
	    result[2] = B[2]

# Usage
magnet = rad.ObjRecMag([0,0,0], [20,20,30], [0,0,1000])
B = RadiaBFieldPy(magnet)
```

Advantages of this approach:
- ✓ No DLL dependency issues
- ✓ No installation required
- ✓ Easy to debug

Disadvantages:
- ✗ Potentially slower than C++ implementation
- ✗ NGSolve CoefficientFunction inheritance may not work correctly

## Next Steps

1. **Test in Conda environment** (recommended, most reliable)
2. If that fails, **switch to pure Python implementation**
3. Or **install MKL on the system**

## References

- NGSolve Official Documentation: https://ngsolve.org/
- Intel oneMKL: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html
- Windows DLL Search Order: https://learn.microsoft.com/en-us/windows/win32/dlls/dynamic-link-library-search-order
