# Material API Implementation Plan

**Date**: 2025-11-11
**Status**: In Progress

---

## Overview

Implementation plan for new material model APIs:
1. `MatLin(ksi)` - Isotropic linear material
2. `MatLin([ksi_par, ksi_perp], [ex, ey, ez])` - Anisotropic linear material
3. `MatPM(Br, Hc, [mx, my, mz])` - Permanent magnet with demagnetization
4. `MatPMLinear(Br, mu_rec, [mx, my, mz])` - PM with recoil permeability

---

## Implementation Strategy

### Phase 1: Update radpy.cpp (Python Bindings)

**File**: `src/python/radpy.cpp`

**Function**: `radia_MatLin()` - Update to handle variable arguments

Current implementation:
```cpp
static PyObject* radia_MatLin(PyObject* self, PyObject* args)
{
	PyObject *oKsi=0, *oMr=0, *oResInd=0;
	// Current: expects 2 arguments (Ksi array, Mr scalar/vector)
}
```

New implementation:
```cpp
static PyObject* radia_MatLin(PyObject* self, PyObject* args)
{
	// Handle 3 cases:
	// Case 1: MatLin(ksi) - single float
	// Case 2: MatLin([ksi_par, ksi_perp], [ex,ey,ez]) - 2 arrays
	// Case 3: MatLin([ksi_par, ksi_perp], [mrx,mry,mrz]) - DEPRECATED but supported
}
```

**New Functions**:
```cpp
static PyObject* radia_MatPM(PyObject* self, PyObject* args)
{
	// MatPM(Br, Hc, [mx,my,mz])
}

static PyObject* radia_MatPMLinear(PyObject* self, PyObject* args)
{
	// MatPMLinear(Br, mu_rec, [mx,my,mz])
}
```

---

### Phase 2: Update rad_c_interface.cpp (C API)

**File**: `src/core/rad_c_interface.cpp`

**New Functions**:
```cpp
// Isotropic linear material
EXP int CALL RadMatLinIso(int* pIndRes, double ksi);

// Anisotropic linear material with easy axis
EXP int CALL RadMatLinAniso(int* pIndRes, double* arKsi, double* arEasyAxis);

// Permanent magnet (Br, Hc)
EXP int CALL RadMatPM(int* pIndRes, double Br, double Hc, double* arMagAxis);

// Permanent magnet (Br, mu_rec)
EXP int CALL RadMatPMLinear(int* pIndRes, double Br, double mu_rec, double* arMagAxis);
```

---

### Phase 3: Material Class Implementation

**Files**:
- `src/core/rad_material_def.h` - Class definitions
- `src/core/rad_material_impl.cpp` - Class implementations

**New Classes**:

#### 1. radTLinearIsotropicMaterial
```cpp
class radTLinearIsotropicMaterial : public radTMaterial {
	double Ksi;  // Magnetic susceptibility

public:
	radTLinearIsotropicMaterial(double InKsi);

	int Type_Material() { return 101; }  // New type ID
	TVector3d M(const TVector3d& H);  // M = ksi * H
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
};
```

#### 2. radTPermanentMagnet (Br/Hc)
```cpp
class radTPermanentMagnet : public radTMaterial {
	double Br;         // Residual flux density [T]
	double Hc;         // Coercivity [A/m]
	TVector3d MagAxis; // Easy magnetization axis (normalized)
	double mu_rec;     // Recoil permeability (calculated)

public:
	radTPermanentMagnet(double InBr, double InHc, const TVector3d& InMagAxis);

	int Type_Material() { return 102; }  // New type ID
	TVector3d M(const TVector3d& H);  // Linear demagnetization curve
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
};
```

#### 3. radTPermanentMagnetLinear (Br/μ_rec)
```cpp
class radTPermanentMagnetLinear : public radTMaterial {
	double Br;         // Residual flux density [T]
	double mu_rec;     // Recoil permeability
	TVector3d MagAxis; // Easy magnetization axis (normalized)

public:
	radTPermanentMagnetLinear(double InBr, double Inmu_rec, const TVector3d& InMagAxis);

	int Type_Material() { return 103; }  // New type ID
	TVector3d M(const TVector3d& H);  // Linear demagnetization curve
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
};
```

---

## Implementation Details

### Isotropic Linear Material

**Physics**:
- Magnetization: M = χ · H
- Permeability: μ_r = 1 + χ
- Isotropic: Same response in all directions

**Implementation**:
```cpp
TVector3d radTLinearIsotropicMaterial::M(const TVector3d& H)
{
	return Ksi * H;
}

void radTLinearIsotropicMaterial::DefineInstantKsiTensor(
	const TVector3d& H, TMatrix3d& InstKsiTensor, TVector3d& InstMr)
{
	InstKsiTensor.Str0.x = Ksi; InstKsiTensor.Str0.y = 0; InstKsiTensor.Str0.z = 0;
	InstKsiTensor.Str1.x = 0; InstKsiTensor.Str1.y = Ksi; InstKsiTensor.Str1.z = 0;
	InstKsiTensor.Str2.x = 0; InstKsiTensor.Str2.y = 0; InstKsiTensor.Str2.z = Ksi;
	InstMr = TVector3d(0, 0, 0);
}
```

---

### Permanent Magnet with Demagnetization

**Physics**:
- B-H curve: B = B_r + μ₀ · μ_rec · H
- M = B/μ₀ - H
- M = (B_r/μ₀) + (μ_rec - 1) · H
- Recoil permeability: μ_rec = B_r / (μ₀ · H_c)

**Constants**:
- μ₀ = 1.25663706212e-6 T/(A/m)

**Implementation**:
```cpp
radTPermanentMagnet::radTPermanentMagnet(double InBr, double InHc, const TVector3d& InMagAxis)
{
	Br = InBr;
	Hc = InHc;
	MagAxis = InMagAxis;
	MagAxis.Normalize();

	const double mu_0 = 1.25663706212e-6;
	mu_rec = Br / (mu_0 * Hc);

	// Set RemMagn = Mr = Br/μ₀
	RemMagn = (Br / mu_0) * MagAxis;
	EasyAxisDefined = 1;
}

TVector3d radTPermanentMagnet::M(const TVector3d& H)
{
	// M = Mr + (μ_rec - 1) · H_parallel
	double H_parallel = H * MagAxis;  // Dot product
	TVector3d M_induced = (mu_rec - 1.0) * H_parallel * MagAxis;
	return RemMagn + M_induced;
}
```

---

## Testing Strategy

### Test 1: Isotropic Linear Material
```python
import radia as rad

# Soft iron (μ_r = 1000)
mat = rad.MatLin(999)
block = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0])
rad.MatApl(block, mat)

# Apply external field
bg = rad.ObjBckg([0, 0, 1e4])  # 10 kA/m
container = rad.ObjCnt([block, bg])
rad.Solve(container, 0.0001, 1000)

# Check magnetization
M = rad.Fld(container, 'm', [0,0,0])
print(f"M = {M} A/m")
# Expected: M ≈ 999 * 10000 = 9.99e6 A/m
```

### Test 2: Permanent Magnet
```python
# NdFeB N52
Br = 1.43  # Tesla
Hc = 876000  # A/m
mat = rad.MatPM(Br, Hc, [0, 0, 1])

magnet = rad.ObjRecMag([0,0,0], [20,20,10], [0,0,0])
rad.MatApl(magnet, mat)
rad.Solve(magnet, 0.0001, 1000)

B = rad.Fld(magnet, 'b', [0,0,15])
print(f"B = {B} T")
```

---

## Migration Guide

### From Old MatLin to New API

**Old (Deprecated)**:
```python
# Permanent magnet with old MatLin
mat = rad.MatLin([0.06, 0.17], [0, 0, 1e6])
```

**New**:
```python
# Use MatPM instead
Br = 1.25  # Tesla
Hc = 950000  # A/m
mat = rad.MatPM(Br, Hc, [0, 0, 1])
```

---

## File Changes Summary

| File | Changes | Lines |
|------|---------|-------|
| src/python/radpy.cpp | Update radia_MatLin(), add radia_MatPM(), radia_MatPMLinear() | ~150 |
| src/core/rad_c_interface.cpp | Add 4 new functions | ~100 |
| src/core/rad_c_interface.h | Add 4 function declarations | ~20 |
| src/core/rad_material_def.h | Add 3 new classes | ~80 |
| src/core/rad_material_impl.cpp | Implement 3 new classes | ~200 |
| **Total** | | **~550** |

---

## Status Tracking

- [ ] Phase 1: Update radpy.cpp
- [ ] Phase 2: Update rad_c_interface.cpp
- [ ] Phase 3: Implement material classes
- [ ] Phase 4: Build and test
- [ ] Phase 5: Create example scripts
- [ ] Phase 6: Update documentation

---

**Author**: Claude Code
**Date**: 2025-11-11
