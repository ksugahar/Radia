# Electromagnet Simulation - C-Type Yoke

Complete 3D magnetostatic simulation of beam steering electromagnet with racetrack coil and C-type magnetic yoke.

Reference: `S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type`

## Overview

This example demonstrates the Cubit -> Netgen -> Radia workflow with **separated mesh generation and simulation**.

### Workflow

```
generate_mesh.py          run_simulation.py
================          =================
Cubit geometry            Load yoke.vol
(from Trelis.jou)              |
    |                     Convert to Radia
Hex mesh                       |
(6x6x6 intervals)         Create coil
    |                          |
Reflect x2                Solve
    |                          |
Transform                 Export VTS
    |
Export Netgen (.vol)
Export STEP (NGSolve)
```

**Reference**: S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\Cubit\6x6x6\Trelis.jou

## Quick Start

### Step 1: Generate Mesh

```bash
python generate_mesh.py
```

Output:
- `yoke.vol` : Netgen mesh file
- `yoke_mesh.vtk` : VTK mesh for ParaView
- `yoke.step` : STEP geometry

### Step 2: Run Simulation

```bash
python run_simulation.py
```

Output:
- `field_distribution.vts` : Magnetic field data

### Step 3: Visualize

```bash
paraview yoke_mesh.vtk field_distribution.vts
```

## Files

| File | Description |
|------|-------------|
| `generate_mesh.py` | Mesh generation (Cubit hex mesh from Trelis.jou) |
| `run_simulation.py` | Radia simulation |
| `yoke.cub` | Cubit database (generated) |
| `yoke.vol` | Netgen hex mesh (generated) |
| `yoke.step` | STEP geometry for NGSolve (generated) |
| `yoke_mesh.vtk` | VTK mesh for ParaView (generated) |
| `field_distribution.vts` | Radia field data (generated) |

## Geometry

### C-Type Yoke

```
          ┌─────┐
          │pole │
    ┌─────┴─────┴─────┐
    │   yoke back     │
    └──┬──────────┬───┘
       │          │
       │   leg    │  (magnetic gap)
       │          │
    ┌──┴──────────┴───┐
    │   yoke back     │
    └─────┬─────┬─────┘
          │pole │
          └─────┘
```

**Dimensions** (mm, from Cubit model):
- Quarter geometry (before reflection):
  - Main leg: 62.5 x 105 x 25
  - Yoke back: 80 x 50 x 25 at (71.25, -27.5, 0)
  - Pole piece: 40 x 100 x 25 at (131.25, -2.5, 0)
- Full model: ~262.5 x 105 x 50 (after reflections)

### Racetrack Coil

- Position: (0, 131.25, 26.25) mm (original ELF coordinates)
- SQRING: 60 x 72.5 mm inner, 35mm height
- Current: -2000 A

## Requirements

- **Coreform Cubit 2025.3+** (for generate_mesh.py only)
- **cubit_mesh_export**: S:\CoreformCubit\01_GitHub
- **NGSolve / Netgen**
- **Radia**

## Advantages of Separation

1. **Re-run simulation without remeshing**: Change coil current, material, etc.
2. **Share mesh files**: Send `yoke.vol` to collaborators
3. **Faster iteration**: Mesh generation is the slow part
4. **Independent debugging**: Test mesh and simulation separately

---

**Last Updated**: 2026-01-22
**Workflow**: Cubit -> Netgen -> Radia (separated)
