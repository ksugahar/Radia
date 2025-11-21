# NGSolve Integration Usage Guide

## Overview

There are two implementations for Radia and NGSolve integration:

| Implementation | Status | Recommendation |
|------|------|--------|
| **Pure Python version** (`radia_ngsolve_py.py`) | ✓ Working | ⭐⭐⭐ Recommended |
| **C++ version** (`radia_ngsolve.pyd`) | △ DLL issues | Conda environment only |

## Pure Python Version Usage (Recommended)

### Basic Usage

```python
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\src\python")
sys.path.insert(0, r"S:\radia\01_GitHub\build\lib\Release")

import radia as rad
import radia_ngsolve_py as radia_ngsolve
from ngsolve import *
from netgen.csg import *

# Create Radia geometry
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1])

# Create Radia field object
B_radia = radia_ngsolve.RadBfield(magnet)

# Generate mesh
geo = CSGeometry()
box = OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50))
geo.Add(box)
mesh = Mesh(geo.GenerateMesh(maxh=15))

# Create GridFunction for visualization
fes = H1(mesh, order=1)
gf_Bz = GridFunction(fes)

# Evaluate field at mesh vertices
for i in range(mesh.nv):
	pnt = mesh.vertices[i].point
	B_vec = B_radia.Evaluate(pnt[0], pnt[1], pnt[2])
	gf_Bz.vec[i] = B_vec[2]  # Z component

# Visualization
import netgen.gui
netgen.gui.Draw(gf_Bz, mesh, "Bz")

# Integration
integral = Integrate(gf_Bz, mesh)
```

### Complete Example

See `examples/ngsolve_integration/visualize_field.py`.

## C++ Version Usage (Conda Environment)

The C++ version is implemented with the same pattern as EMPY_Field, but has DLL dependency issues.

### Prerequisites

```bash
# Create Conda environment
conda create -n radia_ngsolve python=3.12
conda activate radia_ngsolve
conda install -c ngsolve ngsolve
```

### Build

```powershell
cd S:\radia\01_GitHub
.\Build_NGSolve.ps1
```

### Usage

```python
# Execute within Conda environment
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")

import radia as rad
import radia_ngsolve  # C++ version

magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1])
B = radia_ngsolve.RadBfield(magnet)

# Can be used directly in NGSolve (in theory)
from ngsolve import *
Bz = B[2]
```

**Note**: The C++ version currently does not work due to DLL issues. Use of the Pure Python version is recommended.

## Implementation Comparison

### Pure Python Version Features

✓ **Advantages**:
- No DLL dependencies
- No installation required
- Easy to debug
- Cross-platform

✗ **Disadvantages**:
- Manual conversion to GridFunction required
- Slightly slower than C++ version (not an issue in practice)

### C++ Version Features

✓ **Advantages**:
- Can be used directly as NGSolve CoefficientFunction
- Slightly faster
- Same pattern as EMPY_Field

✗ **Disadvantages**:
- DLL dependency issues
- Conda environment required
- Complex build process

## Technical Details

### Pure Python Version Mechanism

`radia_ngsolve_py.py` works as follows:

1. `RadiaBFieldCF` class provides `Evaluate(x, y, z)` method
2. Calls `rad.Fld()` to calculate Radia field
3. User creates GridFunction and evaluates at mesh vertices
4. Used with NGSolve's `Draw()` and `Integrate()`

### C++ Version Mechanism

`radia_ngsolve.cpp` is implemented as follows:

```cpp
class RadiaBFieldCF : public CoefficientFunction
{
public:
	int radia_obj;

	RadiaBFieldCF(int obj)
	    : CoefficientFunction(3), radia_obj(obj) {}

	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                     FlatVector<> result) const override
	{
	    auto pnt = mip.GetPoint();
	    double coords[3] = {pnt[0], pnt[1], pnt[2]};
	    double B[3];
	    int nB = 3;

	    RadFld(B, &nB, radia_obj, "b", coords, 1);

	    result(0) = B[0];
	    result(1) = B[1];
	    result(2) = B[2];
	}
};
```

This follows the same pattern as EMPY_Field's `B_MAGNET_CF`.

## Important Limitations

### ⚠️ Do NOT Use Inside Permanent Magnets

**Critical Limitation**: Both Pure Python and C++ versions of `radia_ngsolve` should **NOT** be used for field calculations **inside** permanent magnets.

**Reason**:
- Both versions evaluate Radia field using `rad.Fld()`
- `rad.Fld()` uses boundary element method (BEM) which is designed for **air regions** only
- BEM is **inaccurate inside permanent magnets** - this is a fundamental limitation of Radia, not a bug
- Therefore, using `radia_ngsolve` inside magnets will return incorrect field values

**Valid Use Cases** (✓ OK):
```python
# GOOD - Mesh OUTSIDE magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1])  # Magnet at origin
mesh = Mesh(Box((50, 50, 50), (100, 100, 100)))  # Mesh far from magnet
B = radia_ngsolve.RadBfield(magnet)
# Field evaluation in air region - OK
```

**Invalid Use Cases** (✗ DO NOT DO):
```python
# BAD - Mesh overlaps with or inside magnet
magnet = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 1])  # Large magnet
mesh = Mesh(Box((-50, -50, -50), (50, 50, 50)))  # Mesh inside magnet
B = radia_ngsolve.RadBfield(magnet)
# Field evaluation inside magnet - WRONG!
```

**Alternative for Interior Fields**:
- Use direct NGSolve FEM solutions with proper material properties
- Do NOT rely on `radia_ngsolve` for calculations inside magnetic materials

---

## Troubleshooting

### Pure Python Version

**Issue**: `No module named 'radia_ngsolve_py'`

**Solution**:
```python
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\src\python")
```

**Issue**: `No module named 'radia'`

**Solution**:
```python
sys.path.insert(0, r"S:\radia\01_GitHub\build\lib\Release")
```

### C++ Version

**Issue**: `DLL load failed while importing radia_ngsolve`

**Solution**: Use the Pure Python version or use a Conda environment.

## Sample Code

### Example 1: Simple Field Evaluation

```python
import radia_ngsolve_py as radia_ngsolve

magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
B = radia_ngsolve.RadBfield(magnet)

# Evaluate at a point
Bx, By, Bz = B.Evaluate(0, 0, 20)
print(f"B at (0,0,20): Bz = {Bz:.6f} T")
```

### Example 2: NGSolve Integration

```python
# Generate mesh
geo = CSGeometry()
box = OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50))
geo.Add(box)
mesh = Mesh(geo.GenerateMesh(maxh=15))

# Create GridFunction
fes = H1(mesh, order=1)
gf_Bz = GridFunction(fes)

# Field evaluation
for i in range(mesh.nv):
	pnt = mesh.vertices[i].point
	B_vec = B.Evaluate(pnt[0], pnt[1], pnt[2])
	gf_Bz.vec[i] = B_vec[2]

# Integration
integral = Integrate(gf_Bz, mesh)
print(f"∫Bz dV = {integral:.6f} T·mm³")
```

### Example 3: Visualization

```python
import netgen.gui

# Visualize Bz component
netgen.gui.Draw(gf_Bz, mesh, "Bz_component")
```

## References

- [NGSOLVE_PYTHON_SOLUTION.md](NGSOLVE_PYTHON_SOLUTION.md) - Pure Python version details
- [NGSOLVE_DLL_ISSUE.md](NGSOLVE_DLL_ISSUE.md) - C++ version DLL issues
- [EMPY_Field Implementation](S:/radia/02_EMPY_Field) - C++ version reference implementation
- [NGSolve Documentation](https://ngsolve.org/)
- [Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
