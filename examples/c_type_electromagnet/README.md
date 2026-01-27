# Electromagnet Simulation - C-Type Yoke

Complete 3D magnetostatic simulation of beam steering electromagnet with racetrack coil and C-type magnetic yoke.

Reference: `S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type`

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
(6x6x6 intervals)             Apply symmetry
    |                              |
Export Netgen (.vol)          Create coil
                                   |
                              Solve
                                   |
                              Export VTS
```

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
- `field_distribution.vts` : Magnetic field data

### Step 3: Visualize

```bash
paraview yoke_quarter.vtu field_distribution.vts
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

**Current Radia Result (1/4 model, 84 elements, mu_r=1000):**
- Gap center Bz: ~0.007 T (requires mesh refinement)

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
| `field_distribution.vts` | Radia field data |

### Utilities (for development)

| File | Description |
|------|-------------|
| `generate_mesh.py` | Full model mesh generation |
| `generate_hex_mesh.py` | Alternative hex mesh generator |
| `verify_geometry.py` | Geometry verification |
| `compare_coil_ngsolve.py` | NGSolve comparison |

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

- **Coreform Cubit 2025.3+** (for mesh generation only)
- **cubit_mesh_export**: `S:\CoreformCubit\01_GitHub`
- **NGSolve / Netgen**
- **Radia**

## TODO

- [ ] Increase mesh density to match ELF reference
- [ ] Add nonlinear B-H curve material
- [ ] Compare with NGSolve FEM solution

---

**Last Updated**: 2026-01-25
**Reference**: S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type
