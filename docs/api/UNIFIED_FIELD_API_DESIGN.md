# Radia Unified Field API Design

**Date**: 2026-01-08 (updated 2026-02-22)
**Status**: Phase 2 completed — `rad.Fld()` dispatches to static and PEEC solvers

## Overview

All field computations use a **single unified API**: `rad.Fld(obj, field_type, point)`.
The function automatically detects the object type and dispatches to the appropriate solver.

For surface BEM (low-frequency EFIE, eddy currents), ngbem is used externally via Python.
Radia does **not** implement its own RWG-EFIE; ngbem's product space provides this.

## Design Philosophy

```python
import numpy as np

# Unified API — rad.Fld() handles everything in Radia core
# Single point: shape (3,)
B = rad.Fld(obj, 'b', np.array([0, 0, 0.1]))

# Batch: shape (N, 3) — auto-detected, parallelized
points = np.array([[0, 0, 0.1], [0, 0, 0.2], [0, 0, 0.3]])
B_batch = rad.Fld(obj, 'b', points)   # returns (N, 3) array
```

## Architecture

```
                         User API
    ┌─────────────────────────────────────────────────────────┐
    │                                                         │
    │     rad.Fld(obj, field_type, point)                    │
    │                                                         │
    │     - Detects object type automatically                │
    │     - Returns real for DC, complex for AC              │
    │     - Handles containers with mixed object types       │
    │                                                         │
    └────────────────────────┬────────────────────────────────┘
                             │
                             ▼
    ┌─────────────────────────────────────────────────────────┐
    │              Field Source Dispatcher (C++)              │
    │                                                         │
    │   radTg3d?          -> Static field (Biot-Savart/MSC)  │
    │   radTConductor?    -> AC field (PEEC Loop-Star)       │
    │   ObjCnt?           -> Sum all sources                 │
    │                                                         │
    └────────────────────────┬────────────────────────────────┘
                             │
              ┌──────────────┴──────────────┐
              │                             │
              ▼                             ▼
    ┌───────────────────┐         ┌───────────────────┐
    │ Static Magneto    │         │ PEEC Conductor    │
    │ (Biot-Savart/MSC) │         │ (Loop-Star)       │
    │                   │         │                   │
    │ - ObjRecMag       │         │ - PEECBuilder     │
    │ - ObjHexahedron   │         │ - PEECCircuit     │
    │ - ObjTetrahedron  │         │ - SIBC skin effect│
    │ - MatLin/Nonlin   │         │                   │
    └───────────────────┘         └───────────────────┘


    External (Python, via ngbem):
    ┌───────────────────────────────────────────────────┐
    │ ngbem Product Space (HDivSurface × SurfaceL2)    │
    │                                                   │
    │ - Low-frequency BEM (Weggler EFIE)               │
    │ - Coupled core models (FEM-BEM, Radia MMM, etc.) │
    │ - Port extraction via Schur complement            │
    │                                                   │
    │ See: docs/peec_integration/demos/ngsbem_peec_demo/ │
    └───────────────────────────────────────────────────┘
```

## Implementation Status

### Completed

- **CndFld() removed**: All conductor field computation uses `rad.Fld()`
- **Unified RadFld()**: Detects conductor handles (>= 10000) and dispatches
- **IsConductorHandle()**: Helper function for object type detection
- **ComputeConductorField()**: PEEC conductor field with complex return values

### Return Value Convention

| Object Type | Return Format |
|------------|---------------|
| Magnetic (PM, iron) | `[Bx, By, Bz]` (3 real values) |
| Conductor (AC PEEC) | `[Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]` (6 values) |

### Handle Ranges

| Handle Range | Object Type |
|-------------|-------------|
| 1 - 9999 | Magnetic objects (radTg3d) |
| 10000+ | Conductor objects (radTConductor) |
| 20000+ | SIBC materials |

## BEM Strategy: ngbem (not Radia C++)

**Decision (2026-01)**: Surface BEM is handled by **ngbem** (NGSolve add-on), not by custom C++ in Radia.

**Rationale**:
1. ngbem provides high-order H(div) elements, H-matrix acceleration, singular quadrature
2. Product space (HDivSurface × SurfaceL2) naturally gives Loop-Star decomposition
3. No low-frequency breakdown — condition number O(1) from DC to RF
4. Avoids reimplementing BEM infrastructure (quadrature, H-matrix, etc.)

**Integration point**: Python level. ngbem matrices are assembled in Python, coupled with
Radia MMM via `ngbem_coupled.py`. Field computation uses ngbem's `GridFunction.Evaluate()`.

| Domain | Solver | Language | API |
|--------|--------|----------|-----|
| Static magnets | Radia MMM/MSC | C++ | `rad.Fld()` |
| PEEC conductors | Radia PEEC | C++ | `rad.Fld()` |
| Surface BEM | ngbem | Python | `ngbem_peec.py` |
| Coupled (all) | Python glue | Python | `ngbem_coupled.py` |

## API Reference

### Field Computation (Unified)

| Function | Description |
|----------|-------------|
| `Fld(obj, 'b', point)` | Single point — shape (3,) |
| `Fld(obj, 'b', points)` | Batch — shape (N, 3), auto-detected |
| `Fld(obj, 'h', ...)` | H field (A/m) |
| `Fld(obj, 'a', ...)` | Vector potential A (T*m) |
| `Fld(obj, 'phi', ...)` | Scalar potential (A) |

### PEEC Conductors (Python API)

PEEC conductors use Python-based `PEECBuilder` and `PEECCircuitSolver` (see `peec_topology.py`).
The old C++ conductor APIs (`CndLoop`, `CndSpiral`, `CndWire`) are removed.

## Usage Examples

```python
import radia as rad
import numpy as np

# 1. Static magnetostatics — single point
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1e6])
B_static = rad.Fld(magnet, 'b', np.array([0, 0, 0.15]))

# 2. Batch evaluation — shape (N, 3) auto-detected
points = np.array([[0, 0, 0.1], [0, 0, 0.2], [0, 0, 0.3]])
B_batch = np.asarray(rad.Fld(magnet, 'b', points))  # (3, 3) array

# 3. Surface BEM (via ngbem, not rad.Fld)
# See docs/peec_integration/demos/ngsbem_peec_demo/
```

## References

1. A.E. Ruehli, "Equivalent Circuit Models for 3D Multiconductor Systems," IEEE Trans. MTT, 1974
2. G. Vecchi, "Loop-Star Decomposition of Basis Functions in EFIE," IEEE TAP, 1999
3. S. Weggler, "Stabilized EFIE using product spaces," NGSolve ngbem
