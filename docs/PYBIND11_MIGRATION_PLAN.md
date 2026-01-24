# pybind11 Migration Plan

## Current State (2026-01-24)

### Binding Technology

| Aspect | Current Status |
|--------|----------------|
| **Binding type** | Python C API (NOT pybind11) |
| **Source file** | `src/radia/radpy_pyapi.cpp` |
| **Python API functions** | 131 |
| **C API functions** | 172 |
| **Coverage** | ~76% (most core functions exposed) |

### Function Categories and Status

| Category | C API Functions | Python API | Status |
|----------|-----------------|------------|--------|
| **Object Creation** | ObjRecMag, ObjHexahedron, etc. | 15 | ✅ Complete |
| **Transformations** | TrfTrsl, TrfRot, TrfMlt, etc. | 10 | ✅ Complete |
| **Materials** | MatLin, MatSatIsoTab, MatPM, etc. | 10 | ✅ Complete |
| **Solver** | Solve, SolveNonl, Rlx*, SetHACApKParams | 15 | ✅ Complete |
| **Field Computation** | Fld, FldBatch, FldA, FldPhi, FldVTS | 20 | ✅ Complete |
| **PEEC Conductor** | Cnd* functions | 20 | ✅ Complete |
| **Coupled Solver** | CplMag* functions | 13 | ✅ Complete |
| **Utilities** | UtiDel, UtiDelAll, UtiDmp, etc. | 6 | ✅ Complete |
| **RWG/ESIM** | RadRwg* | 0 | ❌ Not implemented in C++ |

## Functions NOT Yet Implemented (Future Work)

### RWG/ESIM Functions (Induction Heating - Future)

These functions are **not implemented in C++** yet. They are planned for future induction heating applications.

| Function | Category | Description | Status |
|----------|----------|-------------|--------|
| `RadRwgSolverCreate` | RWG | Create RWG solver | Not implemented |
| `RadRwgMeshCreate` | RWG | Create mesh | Not implemented |
| `RadRwgMeshLoop` | RWG | Loop mesh | Not implemented |
| `RadRwgMeshSpiral` | RWG | Spiral mesh | Not implemented |
| `RadRwgMeshDisk` | RWG | Disk mesh | Not implemented |
| `RadRwgMeshCylinder` | RWG | Cylinder mesh | Not implemented |
| `RadRwgMeshRect` | RWG | Rectangle mesh | Not implemented |
| `RadRwgSetFrequency` | RWG | Set frequency | Not implemented |
| `RadRwgSolve` | RWG | Solve | Not implemented |
| `RadRwgSolveInductionHeating` | RWG | Induction heating solve | Not implemented |
| `RadRwgComputeB` | RWG | Compute B field | Not implemented |
| `RadRwgGetImpedance` | RWG | Get impedance | Not implemented |
| `RadRwgGetWorkpiecePower` | RWG | Get workpiece power | Not implemented |

**Note**: RWG (Rao-Wilton-Glisson) basis functions with ESIM (Effective Surface Impedance Method) are planned for induction heating workpiece analysis.

## Why Consider pybind11 Migration?

### Current Python C API Challenges

| Aspect | Python C API | pybind11 |
|--------|--------------|----------|
| Code complexity | High (manual ref counting) | Low (automatic) |
| Type safety | Manual | C++ templates |
| NumPy integration | Manual PyArray handling | Built-in `py::array_t<>` |
| Error handling | Manual exception conversion | Automatic |
| Docstrings | Manual strings | Automatic from C++ |
| Maintenance effort | Error-prone | Safer |
| Code size | ~5000 lines | ~2000 lines expected |

### Migration Benefits

1. **Cleaner code**: pybind11 reduces boilerplate significantly
2. **Type safety**: C++ templates catch type errors at compile time
3. **NumPy integration**: Direct `py::array_t<double>` support
4. **Exception handling**: C++ exceptions automatically become Python exceptions
5. **Modern C++**: Works naturally with C++11/14/17 features

## Migration Strategy

### Option A: Gradual Migration (Recommended)

1. **Keep existing Python C API** (all 131 functions work)
2. **New functions in pybind11** module (e.g., future RWG functions)
3. **Python wrapper** to unify APIs
4. **Gradual refactoring** of existing functions to pybind11

### File Structure (Proposed)

```
src/radia/
  radpy_pyapi.cpp      # Existing Python C API (131 functions)
  radia_pybind.cpp     # NEW: pybind11 bindings for new functions
  __init__.py          # Python wrapper to combine both

# Build produces:
  radia.pyd            # Existing C API module
  _radia_ext.pyd       # NEW: pybind11 module (optional)

# Python usage:
import radia as rad    # __init__.py combines both seamlessly
```

### Option B: Full Migration (Future)

Convert all 131 functions to pybind11 for maximum code cleanliness.
Estimated effort: 2-3 weeks of focused work.

## Already Working Functions (No Migration Needed)

### Field Computation (All Working)

```python
# These all work correctly:
B = rad.Fld(obj, 'b', [x, y, z])           # Single point field
B_batch, H_batch = rad.FldBatch(obj, pts)  # Batch computation
A = rad.FldA(obj, pts)                      # Vector potential
phi = rad.FldPhi(obj, pts)                  # Scalar potential
rad.FldVTS(obj, 'output.vts', ...)         # VTS export
classification = rad.ClassifyPoints(obj, pts)  # Point classification
```

### PEEC Conductor (All Working)

```python
# These all work correctly:
loop = rad.CndLoop([0,0,0], 0.05, [0,0,1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
rad.CndSetFrequency(loop, 1e6)
rad.CndSetCurrent(loop, 100, 0)
rad.CndSolve(loop)
Z = rad.CndGetImpedance(loop)
```

### Coupled PEEC+MMM (All Working)

```python
# These all work correctly:
solver = rad.CplMagCreate(conductor, magnet)
rad.CplMagSetFrequency(solver, 1e6)
rad.CplMagSetMu(solver, 1000, 0)
result = rad.CplMagSolve(solver)
power = rad.CplMagPower(solver)
```

### Materials (All Working)

```python
# These all work correctly:
mat_linear = rad.MatLin(1000)                    # Linear mu_r=1000
mat_nonlin = rad.MatSatIsoTab(BH_DATA)           # Nonlinear B-H curve
mat_pm = rad.MatMagFixed([0, 0, 954930])         # Fixed PM
mat_sibc = rad.MatSIBC(5.8e7, 1.0)               # Surface impedance
```

## Immediate Action Items

### 1. ✅ Verify existing Python API coverage (DONE)

All Priority 1-4 functions from the original plan are already exposed and working.

### 2. Document current API comprehensively

Create complete API reference with examples for all 131 functions.

### 3. Implement RWG functions in C++ (Future)

When RWG/ESIM is needed for induction heating, implement in C++ first, then add Python bindings.

### 4. Consider pybind11 for new development

Any new C++ functions can use pybind11 from the start.

---

**Last Updated**: 2026-01-24
**Status**: Analysis complete - Most functions already working
**Author**: Claude Opus 4.5
