# NGSolve CoefficientFunction Background Field Design

**Date**: 2025-10-31
**Status**: Design Phase
**Authors**: User + Claude Code

## Overview

Extend Radia's background field capability to accept arbitrary spatial distributions defined by NGSolve CoefficientFunctions, enabling:
- Non-uniform background fields
- FEM-calculated external fields
- Analytical field expressions
- Coupling between NGSolve FEM and Radia magnetostatics

## Current Implementation

### Existing: Uniform Background Field

**Class**: `radTBackgroundFieldSource` (src/core/radarccu.h:274-358)

```cpp
class radTBackgroundFieldSource : public radTg3d {
public:
	TVector3d BackgrB;  // Constant field vector

	void B_comp(radTField* FieldPtr) {
	    if(FieldPtr->FieldKey.B_) FieldPtr->B += BackgrB;
	    // ...
	}
};
```

**Python Interface**: `radia.ObjBckg([Bx, By, Bz])`

**Limitation**: Only uniform (constant) fields supported

## Proposed Design

### Architecture: Python Callback Approach

```
┌─────────────────────────────────────────────────────────┐
│ Python Layer                                             │
├─────────────────────────────────────────────────────────┤
│                                                           │
│  NGSolve CF → radia.ObjBckgCF() → Radia Container       │
│      ↓                                 ↓                 │
│  [Registration]                   [Field Calculation]    │
│      ↓                                 ↓                 │
│  Store CF callback              Call Python callback     │
│                                                           │
└─────────────────────────────────────────────────────────┘
	                ↕ PyObject*             ↕ GIL
┌─────────────────────────────────────────────────────────┐
│ C++ Core (Radia)                                         │
├─────────────────────────────────────────────────────────┤
│                                                           │
│  radTCoefficientFunctionFieldSource : public radTg3d     │
│    - PyObject* cf_callback                               │
│    - void B_comp(radTField* FieldPtr)                    │
│        1. Acquire GIL                                    │
│        2. Call cf_callback(x, y, z)                      │
│        3. Get [Bx, By, Bz] result                        │
│        4. Add to FieldPtr->B                             │
│                                                           │
└─────────────────────────────────────────────────────────┘
```

### Design Rationale

**Why Python Callback (not C++ direct)?**

✓ **Pros**:
- Minimal changes to Radia C++ core
- NGSolve and Radia both controlled from Python
- Easy debugging and prototyping
- Flexible: supports CF, Python functions, or any callable
- No build system complications

✗ **Cons**:
- Performance: Python-C++ boundary crossing
- GIL overhead for each evaluation point

**Alternative Rejected**: C++ Direct Integration
- Would require linking NGSolve into Radia core
- Complex build dependencies
- Difficult maintenance
- Overkill for typical use cases (background fields evaluated ~1000s of points)

## Implementation Details

### 1. New C++ Class: `radTCoefficientFunctionFieldSource`

**File**: `src/core/radarccu.h` (add after line 358)

```cpp
//-------------------------------------------------------------------------
// CoefficientFunction-based Background Field Source
//-------------------------------------------------------------------------

class radTCoefficientFunctionFieldSource : public radTg3d {
public:
	PyObject* cf_callback;  // Python callable (CF or function)

	radTCoefficientFunctionFieldSource(PyObject* callback)
	    : cf_callback(callback)
	{
	    if(cf_callback) Py_INCREF(cf_callback);
	}

	virtual ~radTCoefficientFunctionFieldSource()
	{
	    if(cf_callback) {
	        // GIL may not be held during destruction
	        // Should be handled by Python GC
	        Py_DECREF(cf_callback);
	    }
	}

	// Copy constructor for duplication
	radTCoefficientFunctionFieldSource(const radTCoefficientFunctionFieldSource& src)
	    : radTg3d(src), cf_callback(src.cf_callback)
	{
	    if(cf_callback) Py_INCREF(cf_callback);
	}

	void B_comp(radTField* FieldPtr)
	{
	    if(!cf_callback) return;

	    // Need Python GIL to call Python code
	    py::gil_scoped_acquire acquire;

	    try {
	        // Prepare coordinates (Radia uses mm)
	        py::list coords;
	        coords.append(FieldPtr->P.x);
	        coords.append(FieldPtr->P.y);
	        coords.append(FieldPtr->P.z);

	        // Call Python callback: B = callback([x, y, z])
	        py::handle cb_handle = py::handle(cf_callback);
	        py::object result = cb_handle(coords);

	        // Extract [Bx, By, Bz]
	        py::list B_list = result.cast<py::list>();
	        if(py::len(B_list) != 3) {
	            throw std::runtime_error("CF callback must return [Bx, By, Bz]");
	        }

	        TVector3d B_from_cf(
	            B_list[0].cast<double>(),
	            B_list[1].cast<double>(),
	            B_list[2].cast<double>()
	        );

	        // Add to field components
	        if(FieldPtr->FieldKey.B_) FieldPtr->B += B_from_cf;
	        if(FieldPtr->FieldKey.H_) FieldPtr->H += B_from_cf;

	        // Vector potential A (if needed)
	        // For arbitrary B(r), need to compute A from B
	        // This is non-trivial for general fields
	        // For now, skip A computation
	        if(FieldPtr->FieldKey.A_) {
	            // TODO: Implement A computation
	            // May need to store/compute curl^-1(B)
	        }

	    } catch (std::exception &e) {
	        std::cerr << "[CFFieldSource] Evaluation error: " << e.what() << std::endl;
	        std::cerr << "  at point: (" << FieldPtr->P.x << ", "
	                  << FieldPtr->P.y << ", " << FieldPtr->P.z << ")" << std::endl;
	        // Don't add anything on error
	    }
	}

	void B_intComp(radTField* FieldPtr)
	{
	    // Field integral computation
	    if(!cf_callback) return;

	    if(FieldPtr->FieldKey.FinInt_)
	    {
	        // For arbitrary field, need numerical integration
	        // Use trapezoidal or Simpson's rule
	        TVector3d P1 = FieldPtr->P;
	        TVector3d P2 = FieldPtr->NextP;
	        TVector3d D = P2 - P1;
	        double L = sqrt(D.x*D.x + D.y*D.y + D.z*D.z);

	        // Simple trapezoidal: Integral ≈ (B(P1) + B(P2))/2 * L
	        radTField F1 = *FieldPtr;
	        F1.P = P1;
	        F1.FieldKey.FinInt_ = 0;  // Disable infinite integral
	        B_comp(&F1);

	        radTField F2 = *FieldPtr;
	        F2.P = P2;
	        F2.FieldKey.FinInt_ = 0;
	        B_comp(&F2);

	        TVector3d B_avg = (F1.B + F2.B) * 0.5;
	        TVector3d BufIb = L * B_avg;

	        if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
	        if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIb;
	    }

	    // Infinite integral: set to zero (formally infinite for non-localized fields)
	}

	radTg3dGraphPresent* CreateGraphPresent()
	{
	    // No geometry to display for field source
	    return 0;
	}

	void Dump(std::ostream& o, int ShortSign = 0)
	{
	    radTg3d::Dump(o);
	    o << "CoefficientFunction-based background field source";
	    if(ShortSign==1) return;
	    o << endl;
	    o << "   Python callback: " << (cf_callback ? "registered" : "none");
	    o << endl;
	    o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut,
	             map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers,
	             int& gUniqueMapKey, int elemKey)
	{
	    // Binary serialization not supported for Python callbacks
	    // Would need to pickle the callback
	    throw std::runtime_error("Binary dump not supported for CF field source");
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
	    return FinishDuplication(new radTCoefficientFunctionFieldSource(*this), hg);
	}

	int SizeOfThis() {
	    return sizeof(radTCoefficientFunctionFieldSource);
	}
};
```

### 2. Application Layer Function

**File**: `src/core/radapl1.cpp` (add after SetBackgroundFieldSource)

```cpp
//-------------------------------------------------------------------------

int radTApplication::SetCoefficientFunctionFieldSource(PyObject* callback)
{
	if(!callback) return 0;
	if(!PyCallable_Check(callback)) {
	    Send.ErrorMessage("Radia::Error: Callback must be callable");
	    return 0;
	}

	try
	{
	    radThg hg(new radTCoefficientFunctionFieldSource(callback));
	    int ElemKey = AddElementToContainer(hg);
	    if(SendingIsRequired) Send.Int(ElemKey);
	    return ElemKey;
	}
	catch(...)
	{
	    Initialize();
	    return 0;
	}
}
```

**File**: `src/core/radappl.h` (add declaration)

```cpp
int SetCoefficientFunctionFieldSource(PyObject* callback);
```

### 3. C Interface

**File**: `src/core/radinter.cpp` (add function)

```cpp
void CoefficientFunctionFieldSource(PyObject* callback)
{
	rad.SetCoefficientFunctionFieldSource(callback);
}
```

### 4. DLL Entry Point

**File**: `src/lib/radentry.h` (add declaration)

```cpp
int CALL RadObjBckgCF(int* n, PyObject* callback);
```

**File**: `src/lib/radentry.cpp` (add implementation)

```cpp
int CALL RadObjBckgCF(int* n, PyObject* callback)
{
	CoefficientFunctionFieldSource(callback);
	*n = ioBuffer.OutInt();
	return ioBuffer.OutErrorStatus();
}
```

**File**: `src/lib/raddll.def` (add export)

```
EXPORTS
	...
	RadObjBckgCF
```

### 5. Python Binding

**File**: `src/python/radpy.cpp` (add function ~line 864)

```cpp
/************************************************************************//**
* Creates a source of arbitrary background field from callable.
* Callable should accept [x, y, z] in mm and return [Bx, By, Bz] in Tesla.
***************************************************************************/
static PyObject* radia_ObjBckgCF(PyObject* self, PyObject* args)
{
	PyObject *oCallback=0, *oResInd=0;
	try
	{
	    if(!PyArg_ParseTuple(args, "O:ObjBckgCF", &oCallback))
	        throw CombErStr(strEr_BadFuncArg, ": ObjBckgCF");
	    if(oCallback == 0)
	        throw CombErStr(strEr_BadFuncArg, ": ObjBckgCF");

	    if(!PyCallable_Check(oCallback))
	        throw CombErStr(strEr_BadFuncArg,
	            ": ObjBckgCF requires callable (CF or function)");

	    int ind = 0;
	    g_pyParse.ProcRes(RadObjBckgCF(&ind, oCallback));
	    oResInd = Py_BuildValue("i", ind);
	}
	catch(const char* erText)
	{
	    PyErr_SetString(PyExc_RuntimeError, erText);
	}
	return oResInd;
}
```

**Add to method table** (around line 3288):

```cpp
{"ObjBckgCF", radia_ObjBckgCF, METH_VARARGS,
 "ObjBckgCF(callback) creates a source of arbitrary background field. "
 "Callback should accept [x,y,z] in mm and return [Bx,By,Bz] in Tesla."},
```

### 6. Python Wrapper Module

**File**: `src/python/radia_ngsolve_field.py` (new file)

See full implementation in Section 4 above.

Key functions:
- `create_cf_field_source(cf, unit='m')` - From NGSolve CF
- `create_analytical_field(field_func, unit='mm')` - From Python function

## Unit Convention

**Critical**: Coordinate system and unit handling

| Component | Length Unit | Field Unit | Notes |
|-----------|-------------|------------|-------|
| Radia Core | mm | Tesla | All internal calculations |
| NGSolve | m (default) | Tesla | Can be scaled |
| Python Callback | mm (input) | Tesla (output) | Interface standard |

**Wrapper handles conversion**:
```python
# NGSolve CF defined in meters
cf = CF((Bx, By, Bz))  # x, y, z in meters

# Wrapper converts:
def cf_wrapper(coords):  # coords in mm
	x_m = coords[0] / 1000.0  # Convert to meters
	y_m = coords[1] / 1000.0
	z_m = coords[2] / 1000.0
	return cf(x_m, y_m, z_m)  # Evaluate in meters
```

## Usage Examples

### Example 1: Uniform Field (Compatibility Check)

```python
import radia as rad
from radia_ngsolve_field import create_analytical_field

# Old way
b1 = rad.ObjBckg([0, 0, 1.0])

# New way (same result)
def uniform_field(x, y, z):
	return [0, 0, 1.0]

b2 = create_analytical_field(uniform_field)

# Both should give identical results
```

### Example 2: Quadrupole Field

```python
from ngsolve import CF, x, y, z
from radia_ngsolve_field import create_cf_field_source
import radia as rad

# Define quadrupole gradient field
gradient = 10.0  # T/m
Bx = gradient * y
By = gradient * x
Bz = CF(0)

B_quad = CF((Bx, By, Bz))

# Create Radia field source
field_obj = create_cf_field_source(B_quad, unit='m')

# Permanent magnet in quadrupole field
magnet = rad.ObjRecMag([0,0,0], [5,5,5], [0,0,1.2])
container = rad.ObjCnt([magnet, field_obj])

# Solve with external field
rad.MatApl(magnet, rad.MatSatIsoFrm([2000], [1.8], [0.1]))
rad.Solve(container, 0.0001, 10000)

# Check field at origin
B_total = rad.Fld(container, 'b', [0, 0, 0])
print(f"Total field at origin: {B_total}")
```

### Example 3: FEM Coil Field + Permanent Magnet

```python
from ngsolve import *
from radia_ngsolve_field import create_cf_field_source
import radia as rad

# Step 1: Calculate coil field with NGSolve FEM
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.05))
fes = HCurl(mesh, order=2)

# Define current density (e.g., solenoid)
J = CF((0, 0, 1e6))  # A/m^2 in z-direction

# Solve for vector potential A
u = fes.TrialFunction()
v = fes.TestFunction()
a = BilinearForm(fes)
a += curl(u)*curl(v)*dx
a.Assemble()

f = LinearForm(fes)
f += J*v*dx
f.Assemble()

A = GridFunction(fes)
A.vec.data = a.mat.Inverse() * f.vec

# Step 2: Get B = curl(A)
B_coil = CF((A[2].Diff(y) - A[1].Diff(z),
	         A[0].Diff(z) - A[2].Diff(x),
	         A[1].Diff(x) - A[0].Diff(y)))

# Step 3: Use as Radia background field
field_obj = create_cf_field_source(B_coil, unit='m')

# Step 4: Add permanent magnet
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1.0])
rad.MatApl(magnet, rad.MatSatIsoFrm([2000], [1.2], [0.05]))

container = rad.ObjCnt([magnet, field_obj])

# Step 5: Solve with FEM background
rad.Solve(container, 0.0001, 10000)

# Step 6: Compare fields
B_magnet_only = rad.Fld(magnet, 'b', [0, 0, 0])
B_total = rad.Fld(container, 'b', [0, 0, 0])

print(f"Magnet only: {B_magnet_only}")
print(f"Magnet + Coil: {B_total}")
print(f"Coil contribution: {[B_total[i]-B_magnet_only[i] for i in range(3)]}")
```

### Example 4: Spatially Varying Dipole Field

```python
from ngsolve import CF, x, y, z, sqrt
from radia_ngsolve_field import create_cf_field_source
import radia as rad

# Magnetic dipole at origin
# B(r) = (μ₀/4π) * (3(m·r)r/r⁵ - m/r³)

m = [0, 0, 1.0]  # Dipole moment (A·m²)
mu0_4pi = 1e-7  # T·m/A

r = sqrt(x**2 + y**2 + z**2 + 1e-10)  # Avoid singularity
r3 = r**3
r5 = r**5

m_dot_r = m[0]*x + m[1]*y + m[2]*z

Bx = mu0_4pi * (3*m_dot_r*x/r5 - m[0]/r3)
By = mu0_4pi * (3*m_dot_r*y/r5 - m[1]/r3)
Bz = mu0_4pi * (3*m_dot_r*z/r5 - m[2]/r3)

B_dipole = CF((Bx, By, Bz))

field_obj = create_cf_field_source(B_dipole, unit='m')

# Calculate field at various points
for pos in [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]]:
	B = rad.Fld(field_obj, 'b', [p*1000 for p in pos])  # Convert m to mm
	print(f"B at {pos} m: {B} T")
```

## Performance Considerations

### Expected Performance

**Typical Use Case**: Background field evaluated at ~1,000-10,000 points

| Operation | Time per call | For 10k points |
|-----------|---------------|----------------|
| Python callback | ~10 μs | ~0.1 s |
| GIL acquire/release | ~1 μs | ~0.01 s |
| CF evaluation | ~1-10 μs | ~0.01-0.1 s |
| **Total** | ~12-21 μs | **~0.12-0.21 s** |

**Acceptable** for typical magnetostatics calculations where:
- Relaxation solve: 10-1000 iterations @ ~0.1-1 s per iteration
- Background field evaluation: <1% of total time

### Optimization Options (if needed)

1. **Caching**: Store recently evaluated points
2. **Vectorization**: Batch evaluate multiple points
3. **C++ CF Evaluation**: Direct NGSolve C++ (future enhancement)

## Testing Plan

### Unit Tests

1. **Uniform field equivalence**: `ObjBckgCF` ≡ `ObjBckg`
2. **Linear field**: Verify gradient fields
3. **Coordinate system**: mm ↔ m conversion
4. **Error handling**: Invalid callbacks, exceptions

### Integration Tests

1. **With permanent magnets**: External field + magnetization
2. **With saturation**: Nonlinear materials in external field
3. **Field integrals**: Line integrals through varying fields
4. **FEM coupling**: NGSolve solution → Radia background

### Example Scripts

All examples from Section "Usage Examples" should be implemented as test cases.

## Documentation Updates

### Files to Update

1. **README.md**: Add NGSolve CF field feature
2. **NGSOLVE_INTEGRATION.md**: Add CF background field section
3. **Python docstrings**: Complete API documentation
4. **examples/**: Add example scripts

### API Documentation

```python
radia.ObjBckgCF(callback)
	"""
	Create arbitrary background field source from callable.

	Parameters
	----------
	callback : callable
	    Function accepting [x, y, z] in mm, returning [Bx, By, Bz] in Tesla
	    Signature: callback([x, y, z]) -> [Bx, By, Bz]

	Returns
	-------
	int : Radia object key

	See Also
	--------
	radia.ObjBckg : Uniform background field
	radia_ngsolve_field.create_cf_field_source : NGSolve CF wrapper
	"""
```

## Implementation Steps

1. ✓ Design architecture
2. ✓ Create design document (this file)
3. ☐ Implement C++ class in radarccu.h
4. ☐ Add application function in radapl1.cpp
5. ☐ Add C interface in radinter.cpp
6. ☐ Add DLL entry in radentry.cpp
7. ☐ Add Python binding in radpy.cpp
8. ☐ Create wrapper module radia_ngsolve_field.py
9. ☐ Write unit tests
10. ☐ Write integration tests
11. ☐ Create example scripts
12. ☐ Update documentation
13. ☐ Build and test
14. ☐ Commit and push

## Open Questions

1. **Vector potential A**: How to compute for arbitrary B(r)?
   - Option: Ignore for background fields (most use cases only need B)
   - Option: Require user to provide A separately
   - Option: Numerical integration (complex)

2. **Serialization**: How to save/load CF field sources?
   - Binary dump not supported (Python callbacks)
   - Could save CF expression as string
   - Or mark as "not serializable"

3. **Performance**: Is callback overhead acceptable?
   - Benchmark with real use cases
   - Consider caching strategy if needed

4. **Thread safety**: Multiple Radia solves in parallel?
   - GIL protects Python side
   - Each solve should have separate containers

## Related Work

- **rad_ngsolve.cpp**: Opposite direction (Radia → NGSolve CF)
- **radTBackgroundFieldSource**: Existing uniform field implementation
- **Python GIL**: Experience from rad_ngsolve.cpp callback

## References

- Radia Documentation: https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia
- NGSolve: https://ngsolve.org/
- pybind11: https://pybind11.readthedocs.io/

---

**Status**: Ready for implementation
**Next Step**: Implement radTCoefficientFunctionFieldSource in radarccu.h
