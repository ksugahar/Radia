# Tetrahedral Mesh Support Implementation

## Overview

Radia now supports tetrahedral meshes with two computation methods:

1. **Standard Polygon Method** (default): Uses existing Radia polygon field computation
2. **Centroid Charge Method**: Uses analytical formulas with centroid point charge cancellation

## Method Selection

Set environment variable before importing Radia:

```python
import os
os.environ['RADIA_TETRA_METHOD'] = 'CENTROID'  # Use centroid charge method
# or
os.environ['RADIA_TETRA_METHOD'] = 'STANDARD'  # Use standard polygon method (default)

import radia as rad
```

## Implementation Details

### File Changes

#### `src/core/rad_polyhedron.h`

```cpp
// Line ~277
bool IsTetrahedron() const { return AmOfFaces == 4; }
void B_comp_tetrahedron_centroid(radTField*);
```

#### `src/core/rad_polyhedron.cpp`

**Headers** (add after line 20):
```cpp
#include <cstdlib>  // std::getenv
#include <string>   // std::string
```

**B_comp_frM() modification** (line ~674):
```cpp
void radTPolyhedron::B_comp_frM(radTField* FieldPtr)
{
	// Tetrahedral method selection via environment variable
	if(IsTetrahedron())
	{
		const char* method_env = std::getenv("RADIA_TETRA_METHOD");
		if(method_env != nullptr && std::string(method_env) == "CENTROID")
		{
			B_comp_tetrahedron_centroid(FieldPtr);
			return;
		}
		// Fall through to standard polygon method
	}

	// Standard implementation continues...
	TVector3d Zero(0.,0.,0.);
	short PointIsInside = 1;
	// ... rest unchanged
}
```

**B_comp_tetrahedron_centroid()**:
- Line ~437: Function implementation with centroid charge cancellation
- Line ~439-441: Update comments
- Line ~462: Add `double total_magnetic_charge = 0.0;`
- Line ~529 (in face loop): Add `total_magnetic_charge += W * Area;`
- Lines ~663-669: Replace with centroid charge cancellation code

### Centroid Charge Cancellation Code

Replace lines 663-669 in `B_comp_tetrahedron_centroid()`:

```cpp
	// Centroid point charge cancellation
	// For each tetrahedron, place point charge at centroid to cancel face contributions
	TVector3d R = ObsPt - EEC;  // Vector from centroid to observation point
	double R_mag = sqrt(R.x*R.x + R.y*R.y + R.z*R.z);

	if(R_mag > EPS)
	{
		// Point charge formula: H = Q / (4π r³) × r
		// Q = sum of face charges (with area)
		double factor = total_magnetic_charge / (4.0 * PI * R_mag * R_mag * R_mag);
		TVector3d H_centroid;
		H_centroid.x = factor * R.x;
		H_centroid.y = factor * R.y;
		H_centroid.z = factor * R.z;

		// Subtract centroid contribution (cancels interior faces)
		H_total -= H_centroid;
	}

	// Add to field pointer
	if(LocFieldKey.H_) FieldPtr->H += H_total;
	if(LocFieldKey.B_) FieldPtr->B += H_total;
}
```

## Testing

### Test Scripts

**test_both_methods.py**:
```bash
cd examples/ngsolve_integration/mesh_magnetization_import
python test_both_methods.py
```

### Expected Results

- **Reference (ObjDivMag)**: |H| ≈ 0.034
- **Standard method**: May overflow for dense meshes (interior face double-counting)
- **Centroid method**: Should match reference within 10-20% (depends on mesh quality)

## Build Instructions

```bash
# Clean build
cd S:\Radia\01_GitHub\build
cmake --build . --config Release --target radia --clean-first

# Test
cd ../examples/ngsolve_integration/mesh_magnetization_import
python test_both_methods.py
```

## Known Limitations

1. **Standard method**: Interior faces double-counted → error scales with mesh density
2. **Centroid method**: Requires uniform magnetization within each tetrahedron
3. **Performance**: Centroid method adds point charge computation overhead (~10%)

## Future Work

- Implement per-element magnetization variation
- Optimize centroid charge computation
- Add adaptive method selection based on mesh properties

---

**Last Updated**: 2025-11-24
**Author**: Radia Development Team
