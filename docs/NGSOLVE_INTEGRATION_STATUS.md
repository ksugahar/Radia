# NGSolve Integration Status

## Implementation Status

### Completed Items ✓

1. **Source Code Implementation**
   - `src/python/rad_ngsolve.cpp` - Fully implemented
   - `RadBfield`, `RadHfield`, `RadAfield` classes
   - NGSolve CoefficientFunction interface implementation

2. **CMake Configuration**
   - NGSolve auto-detection
   - Using `add_ngsolve_python_module`
   - All dependencies configured

3. **Build**
   - `rad_ngsolve.pyd` built successfully
   - Location: `build/Release/rad_ngsolve.pyd`
   - Size: 1.3 MB

4. **Documentation**
   - User guide created
   - Technical documentation created
   - Test scripts created

### Known Issues ⚠

**DLL Dependency Problem**

When importing `rad_ngsolve.pyd`, the following error occurs:
```
ImportError: DLL load failed while importing rad_ngsolve: The specified module could not be found.
```

**Cause:**
- `rad_ngsolve.pyd` depends on `ngsolve.dll`, `libngsolve.dll`, `ngcore.dll`, etc.
- These DLLs are not included in the Windows DLL search path

**Impact:**
- Module is built but cannot be loaded at runtime
- Cannot find DLLs within NGSolve's Python package

## Solutions

### Method 1: Environment Variable Setup (Recommended)

Add DLL paths to system PATH:

```powershell
# Execute in PowerShell
$env:PATH = "C:\Program Files\Python312\Lib\site-packages\ngsolve;$env:PATH"
$env:PATH = "C:\Program Files\Python312\Lib\site-packages\netgen;$env:PATH"
```

Then launch Python:
```python
import sys
sys.path.insert(0, r'S:\radia\01_GitHub\build\Release')
import rad_ngsolve  # Should work now
```

### Method 2: Dynamically Add DLL Directory

Within Python script:
```python
import os
import sys

# For Windows 10/11
if hasattr(os, 'add_dll_directory'):
	os.add_dll_directory(r'C:\Program Files\Python312\Lib\site-packages\ngsolve')
	os.add_dll_directory(r'C:\Program Files\Python312\Lib\site-packages\netgen')

sys.path.insert(0, r'S:\radia\01_GitHub\build\Release')
import rad_ngsolve
```

### Method 3: Copy Dependent DLLs

Copy required DLLs to `build/Release`:
```powershell
Copy-Item "C:\Program Files\Python312\Lib\site-packages\ngsolve\*.dll" `
	      "S:\radia\01_GitHub\build\Release\"
Copy-Item "C:\Program Files\Python312\Lib\site-packages\netgen\*.dll" `
	      "S:\radia\01_GitHub\build\Release\"
```

### Method 4: Install in Same Location as NGSolve

Most reliable method:
```powershell
# Copy rad_ngsolve.pyd to NGSolve's site-packages
Copy-Item "S:\radia\01_GitHub\build\Release\rad_ngsolve.pyd" `
	      "C:\Program Files\Python312\Lib\site-packages\"
```

Then:
```python
import rad_ngsolve  # Can import from anywhere
```

## Technical Implementation Details

### Architecture

```
Python Script
	↓
rad_ngsolve.pyd (Python module)
	↓
ngfem::RadiaBFieldCF (CoefficientFunction)
	↓
RadFld() (Radia C API)
	↓
Radia Core (Magnetic field calculation)
```

### Function Signatures

```python
# Python-side API
rad_ngsolve.RadBfield(radia_obj: int) -> CoefficientFunction
rad_ngsolve.RadHfield(radia_obj: int) -> CoefficientFunction
rad_ngsolve.RadAfield(radia_obj: int) -> CoefficientFunction
```

### C++ Implementation Key Points

```cpp
namespace ngfem {
	class RadiaBFieldCF : public CoefficientFunction {
	    int radia_obj;  // Radia object index

	    virtual void Evaluate(
	        const BaseMappedIntegrationPoint& mip,
	        FlatVector<> result) const override
	    {
	        auto pnt = mip.GetPoint();
	        double coords[3] = {pnt[0], pnt[1], pnt[2]};
	        double B[3];
	        RadFld(B, &nB, radia_obj, "b", coords, 1);
	        result(0) = B[0];
	        result(1) = B[1];
	        result(2) = B[2];
	    }
	};
}
```

## Testing

### Basic Test (After DLL Issue Resolution)

```python
import radia as rad
import rad_ngsolve

# Radia geometry
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1000])
rad.Solve(magnet, 0.0001, 10000)

# Create CoefficientFunction
B = rad_ngsolve.RadBfield(magnet)
print(type(B))  # <class 'ngsolve.fem.CoefficientFunction'>
```

### NGSolve Integration Test

```python
from ngsolve import *
from netgen.csg import *

# Mesh
geo = CSGeometry()
box = OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50))
geo.Add(box)
mesh = Mesh(geo.GenerateMesh(maxh=10))

# Use Radia field in NGSolve
B_integral = Integrate(B, mesh)
Draw(B, mesh, "B_field")
```

## Build Instructions

```powershell
# 1. Install NGSolve and Netgen
conda install -c ngsolve ngsolve

# 2. Run build script
.\Build_NGSolve.ps1

# 3. Verify output
ls build\Release\rad_ngsolve.pyd
```

## Summary

**Implementation**: ✓ Complete
**Build**: ✓ Success
**Execution**: ⚠ DLL dependency configuration required

Once the DLL issue is resolved, Radia's magnetic field calculations can be directly used in NGSolve FEM simulations.

## Reference Links

- [NGSolve Documentation](https://ngsolve.org/)
- [Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
- [Windows DLL Search Order](https://learn.microsoft.com/en-us/windows/win32/dlls/dynamic-link-library-search-order)
