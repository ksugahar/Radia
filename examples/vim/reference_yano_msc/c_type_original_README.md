# Electromagnet Simulation - C-Type Yoke

Complete 3D magnetostatic simulation of beam steering electromagnet with racetrack coil and C-type magnetic yoke.

Reference: private CEFC-2020 C-Type model archive.

## Overview

This example demonstrates the Cubit -> Netgen -> Radia workflow with **separated mesh generation and simulation**.

### Workflow

```
generate_quarter_mesh.py      run_simulation.py
========================      =================
Cubit geometry                Load yoke_quarter.vol
(1/4 model)                        |
    |                         Convert to Radia
Hex mesh                           |
(6x6x6 intervals)             Apply Image symmetry
    |                              |
Export Netgen (.vol)          Create coil
                                   |
                              Solve with image='+x-z'
                                   |
                              Export VTS
```

## Image Symmetry API (2026-01-31)

The new unified Image symmetry API is implemented and verified against ELF_MAGIC:

```python
import radia as rad

# Create quarter model geometry (no TrfMlt needed)
yoke = rad.ObjCnt(hex_elements)

# Solve with Image symmetry - quarter model -> full model
rad.Solve(yoke, 0.0001, 100, 0, image='+x-z')

# Or build matrix first for inspection
handle = rad.BuildMatrix(yoke, image='+x-z')
matrix, dof = rad.GetInteractMatrix(handle)
```

### Image Parameter Format

| Parameter | Meaning |
|-----------|---------|
| `+x` | Symmetric mirror across X=0 plane |
| `-x` | Antisymmetric mirror across X=0 plane |
| `+z` | Symmetric mirror across Z=0 plane |
| `-z` | Antisymmetric mirror across Z=0 plane |
| `+x-z` | Both mirrors (quarter model) |
| `+x+y-z` | Three mirrors (eighth model) |

### Verification Results (mu=1000, 13 elements)

**Matrix Comparison:**
- Radia vs ELF matrix relative error: ~7.2% (after convention fix)
- Matrix symmetry verified

**Field Results:**
- Image API correctly changes magnetization solution
- Field at origin with `image='+x-z'`: Bz = 246.8 mT
- Field without Image: Bz = -66.3 mT
- The difference confirms Image symmetry is working

**Expected Difference:**
- ~0.82% field difference expected due to coil modeling:
  - ELF: Discretized coil elements (MCL8T)
  - Radia: Analytical coil (ObjArcCur)

### Key Test Files

| File | Description |
|------|-------------|
| `mu=1000/quarter/test_all_solvers_ima.py` | All solvers with IMA (LU/BiCGSTAB/HACApK) |
| `mu=1000/quarter/test_hacapk_ima_transition.py` | HACApK state management verification |
| `nonlinear/full/LU/test_all_solvers.py` | Nonlinear solver comparison |
| `nonlinear/full/LU/verify_full_nonlinear.py` | Full nonlinear model verification |

## Quick Start

### Step 1: Generate Mesh (requires Cubit)

```bash
python generate_quarter_mesh.py
```

Output:
- `yoke_quarter.vol` : Netgen mesh file (1/4 model)
- `yoke_quarter.vtu` : VTK mesh for ParaView

### Step 2: Run Simulation

```bash
python run_simulation.py
```

Output:

### Step 3: Visualize

```bash
```

## Validation Results

### Reference: ELF_MAGIC Comparison (CEFC 2020)

Gap center Bz field [T] for different mesh densities and permeabilities:

**Linear Material (mu_r constant):**

| Solver | Mesh | Elements | mu_r=100 | mu_r=1000 | mu_r=10000 |
|--------|------|----------|----------|-----------|------------|
| ELF    | 1x1x1 | 52 | 0.0403 T | 0.0513 T | 0.0527 T |
| ELF    | R288 | 288 | 0.1233 T | 0.2318 T | 0.2490 T |
| ELF    | V4672 | 4672 | 0.1222 T | 0.2308 T | 0.2516 T |
| Radia  | R288 | 288 | 0.1232 T | 0.2320 T | 0.2492 T |
| Radia  | V4672 | 4672 | 0.1220 T | 0.2310 T | 0.2519 T |

**Nonlinear Material (1000 AT excitation):**

| Solver | Mesh | Elements | Bz [T] |
|--------|------|----------|--------|
| ELF    | 1x1x1 | 52 | 0.0264 T |
| ELF    | 6x6x6 | 4800 | 0.1019 T |
| ELF    | R288 | 288 | 0.1246 T |
| ELF    | R3856 | 3856 | 0.1250 T |
| ELF    | V4864 | 4864 | 0.1261 T |

### Notes

- R-mesh: Structured rectangular mesh
- V-mesh: Tetrahedral mesh
- ELF uses MMM (Magnetic Moment Method) with 8-node hexahedra
- Radia uses MSC (Magnetic Surface Charge) with 6-DOF hexahedra

## Files

### Main Workflow

| File | Description |
|------|-------------|
| `generate_quarter_mesh.py` | 1/4 model mesh generation (Cubit) |
| `run_simulation.py` | Radia simulation with symmetry |
| `README.md` | This file |

### Generated Files

| File | Description |
|------|-------------|
| `yoke_quarter.vol` | Netgen mesh (1/4 model) |
| `yoke_quarter.vtu` | VTK mesh for visualization |

### Subdirectories

| Directory | Description |
|-----------|-------------|
| `mu=1000/` | Linear material validation (mu_r=1000) |
| `mu=1000/full/` | Full model tests |
| `mu=1000/quarter/` | Quarter model IMA symmetry tests |
| `mu=1000/single/` | Single element verification |
| `nonlinear/` | Nonlinear material (B-H curve) tests |
| `nonlinear/full/` | Full model with LU/BiCGSTAB/HACApK solvers |
| `nonlinear/quarter/` | Quarter model with all solvers |

## Geometry

### C-Type Yoke

```
          +-----+
          |pole |
    +-----+-----+-----+
    |   yoke back     |
    +--+----------+---+
       |          |
       |   leg    |  (magnetic gap)
       |          |
    +--+----------+---+
    |   yoke back     |
    +-----+-----+-----+
          |pole |
          +-----+
```

**Dimensions** (mm):
- Quarter geometry (before reflection):
  - Main leg: 62.5 x 105 x 25
  - Yoke back: 80 x 50 x 25
  - Pole piece: 40 x 100 x 25
- Full model: ~262.5 x 105 x 50 (after X and Z reflections)

### Racetrack Coil

- Position: (-12.5, 131.25, 78.75) mm (after coordinate transformation)
- Inner dimensions: 60 x 72.5 mm
- Corner radii: 5 mm (inner), 40 mm (outer)
- Height: 105 mm
- Current: -2000 A

## Requirements

- **Coreform Cubit 2025.8+** (for mesh generation only)
- **cubit_mesh_export**: `S:\CoreformCubit\01_GitHub`
- **NGSolve / Netgen**
- **Radia**

## Changelog

### 2026-02-05
- Cleaned up debug/temporary scripts (114 -> 67 files)
- Reorganized folder structure: linear in mu=1000/, nonlinear in nonlinear/
- Fixed HACApK thread-local cache invalidation for IMA transitions
- Verified all solvers work with quarter model (+x-z)

### 2026-01-31
- Implemented new Image symmetry API: `rad.Solve(..., image='+x-z')`
- Removed old TrfMlt-based symmetry approach
- Verified against ELF_MAGIC quarter model (7.2% matrix difference)

### 2026-01-29
- Fixed MSC hexahedron mirror symmetry using reciprocity approach

---

**Last Updated**: 2026-02-05
**Reference**: private CEFC-2020 C-Type model archive.
