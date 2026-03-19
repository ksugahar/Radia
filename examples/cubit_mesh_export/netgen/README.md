# Netgen/NGSolve Examples

This folder contains examples for exporting Cubit meshes to Netgen/NGSolve with high-order curving support.

## Recommended Workflow: Name-based + SetGeomInfo

The best approach for accurate high-order curving combines **name-based face mapping** with **SetGeomInfo API**:

```
1. OCC: Create geometry
2. OCC: name_occ_faces() - assign unique names to faces
3. OCC: Export STEP (names preserved in STEP)
4. Cubit: Import STEP (names visible as surface names)
5. Cubit: Mesh
6. export_netgen_with_names() - correct face indices
7. SetGeomInfo: Set UV parameters for curved surfaces
8. mesh.Curve(order) - high-order curving works correctly!
```

## Test Results (2026-01-23)

### Name-based + SetGeomInfo (RECOMMENDED)

| Example | Shape | Curve(2) Error | Curve(3) Error |
|---------|-------|----------------|----------------|
| Complex | Brick with hole | **0.0021%** | **0.0004%** |
| Cylinder | R=0.5, H=2 | 0.0027% | 0.0006% |
| Sphere | R=0.5 | 0.0027% | 0.0004% |
| Torus | R_major=1, R_minor=0.3 | 0.0010% | 0.0003% |

**Key achievement**: Complex geometries now achieve **Netgen-native accuracy** (<0.001%)!

## Example Files

### Name-based + SetGeomInfo

| File | Description | Status |
|------|-------------|--------|
| `netgen_complex_named.py` | **Complex geometry (Boolean ops)** | PASS |
| `netgen_cylinder_setgeominfo.py` | Cylinder (tet mesh) | PASS |
| `netgen_hex_cylinder_setgeominfo.py` | Cylinder (hex mesh) | PASS |
| `netgen_sphere_setgeominfo.py` | Sphere | PASS |
| `netgen_torus_setgeominfo.py` | Torus | PASS |
| `netgen_cone_setgeominfo.py` | Cone | PASS |

### Test Files (in ../tests/)

| File | Description | Status |
|------|-------------|--------|
| `test_setgeominfo.py` | SetGeomInfo API unit test | PASS |
| `test_curve_workflow.py` | Compare curving methods | PASS |
| `test_setgeominfo_uv.py` | SetGeomInfo with UV from geometry | PASS |

## Code Example

```python
import cubit_mesh_export
from netgen.occ import OCCGeometry, Box, Cylinder, gp_Pnt, gp_Ax2, gp_Dir
from ngsolve import Mesh

# 1. Create geometry in OCC
brick = Box(gp_Pnt(-1,-1,-1), gp_Pnt(1,1,1))
cyl = Cylinder(gp_Ax2(gp_Pnt(0,0,-2), gp_Dir(0,0,1)), 0.3, 4)
shape = brick - cyl

# 2. Name faces (critical for correct mapping!)
cubit_mesh_export.name_occ_faces(shape)

# 3. Export STEP
shape.WriteStep("geometry.step")

# 4. Load geometry reference
geo = OCCGeometry("geometry.step")

# 5. Import into Cubit and mesh
cubit.cmd('import step "geometry.step" noheal')
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# 6. Export with name-based mapping
ngmesh = cubit_mesh_export.export_netgen_with_names(cubit, geo)

# 7. Set UV for curved surfaces
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=0.3, height=2.0, axis='z')

# 8. High-order curving
mesh = Mesh(ngmesh)
mesh.Curve(2)  # Now works correctly!
```

## Available Functions

### Name-based Workflow (RECOMMENDED)

```python
# Name OCC faces before STEP export
name_occ_faces(shape, prefix="occ_face_")

# Export with correct face mapping
export_netgen_with_names(cubit, geometry) -> ngmesh
```

### SetGeomInfo Functions

```python
# Set UV for analytic surfaces
set_cylinder_geominfo(ngmesh, radius, height, center=(0,0,0), axis='z')
set_sphere_geominfo(ngmesh, radius, center=(0,0,0))
set_torus_geominfo(ngmesh, major_radius, minor_radius, center=(0,0,0), axis='z')
set_cone_geominfo(ngmesh, base_radius, height, center=(0,0,0), axis='z')
```

## Why This Works

### The Problem

- Cubit uses ACIS kernel, Netgen uses OpenCASCADE (OCC)
- STEP exchange doesn't preserve face indices
- `mesh.Curve()` requires correct face indices + UV parameters
- Without both, curving fails or produces wrong results

### The Solution

1. **Name-based mapping**: OCC face names survive STEP round-trip
2. **SetGeomInfo**: Provides UV parameters for curved surfaces
3. **Together**: Achieves Netgen-native accuracy for any geometry

## Requirements

- NGSolve: Local build from `S:/NGSolve/01_GitHub/install_ksugahar/`
  - Includes PR #232 (SetGeomInfo API)
  - Includes PR #231 (Periodic BC fix)
- Coreform Cubit 2025.3+

## See Also

- [CLAUDE.md](../../../CLAUDE.md) - Development policy
- [NGSolve Forum: SetGeomInfo API](https://forum.ngsolve.org/t/feature-request-python-api-for-high-order-curving-of-externally-imported-meshes/3810)
