# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and refactoring policies for the Radia project when working with Claude Code.

## Radia Solver Method: MMM (Magnetic Moment Method)

**IMPORTANT**: Radia uses the **MMM (Magnetic Moment Method)**, NOT BEM (Boundary Element Method).

**Terminology**:
- Correct: "Radia MMM solver", "MMM computation", "MMM field evaluation"
- Incorrect: "Radia BEM", "BEM solver", "boundary element method"

**Context**:
- MMM represents magnetic objects as distributions of magnetic moments
- This differs from BEM which uses surface integral equations
- All references to "BEM" in Radia-related documentation should be corrected to "MMM"

---

## Polyhedral Element Field Computation (MSC)

### Magnetic Surface Charge (MSC) Method

Both tetrahedral and hexahedral polyhedral elements use the **Magnetic Surface Charge (MSC)** method for field computation. This method computes the magnetic field from each triangular face using closed-form formulas based on the surface magnetic charge density (sigma = M dot n).

**Supported Element Types**:
- **Tetrahedra (4 faces)**: 4 triangular faces
- **Hexahedra (6 faces)**: 6 quadrilateral faces, each split into 2 triangles

**Key Features**:
- Uses **global coordinates** directly (no local coordinate transformations)
- Computes field using **solid angle integration** formula
- Handles **outward normal orientation** following ELF_MAGIC convention
- **ELF_MAGIC compatible**: Uses same 6-face hexahedral format (not 12-face triangular)

**Implementation** (`src/core/rad_polyhedron.cpp`):
- `B_comp_tetrahedron_MSC()`: 4-face tetrahedral field computation
- `B_comp_hexahedron_MSC()`: 6-face hexahedral field computation (quad -> 2 triangles)
- Uses `RadFieldFromTriangleFaceGlobal()` from `rad_poly_analytical.cpp`

**Usage**:

```python
import radia as rad
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits("m")

# Import tetrahedral mesh (automatically uses MSC method)
mag_obj = netgen_mesh_to_radia(mesh,
                                material={"magnetization": [0, 0, 0]},
                                units="m")

# Apply material and solve
mat = rad.MatLin(999.0)  # mu_r = 1000
rad.MatApl(mag_obj, mat)
rad.Solve(mag_obj, 0.0001, 1000)
```

**Notes**:
- Both tetrahedral and hexahedral meshes work with linear and nonlinear materials
- LU solver (Method 0) and BiCGSTAB (Method 1) both work with polyhedral elements
- Hexahedral MSC matches ELF_MAGIC results for single-element cases

---

**Last Updated**: 2025-12-11
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
