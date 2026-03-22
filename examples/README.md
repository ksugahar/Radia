# Radia Examples

Comprehensive collection of Radia examples demonstrating magnetic field computation, material properties, solver methods, PEEC conductor modeling, and integration with NGSolve FEM.

**Total:** 356 Python scripts across 20 directories

---

## Quick Start

```bash
cd examples/simple_problems
python arc_current_with_magnet.py
paraview arc_current_with_magnet.vts
```

---

## Directory Overview

### Magnetostatics (Beginner)

| Directory | Scripts | Description |
|-----------|---------|-------------|
| [simple_problems/](simple_problems/) | 5 | Basic Radia API: magnets, coils, materials, field computation |
| [smco_magnet_array/](smco_magnet_array/) | 1 | SmCo permanent magnet array design |
| [vtk_export/](vtk_export/) | 1 | VTS field export for ParaView visualization |

### Magnetostatics (Intermediate)

| Directory | Scripts | Description |
|-----------|---------|-------------|
| [background_fields/](background_fields/) | 4 | External background fields with magnetizable materials |
| [complex_coil_geometry/](complex_coil_geometry/) | 3 | Multi-segment coils using CoilBuilder API |
| [c_type_electromagnet/](c_type_electromagnet/) | 43 | C-type electromagnet: racetrack coil + yoke with nonlinear materials |
| [cube_uniform_field/](cube_uniform_field/) | 8 | Cube benchmark: hex/tetra mesh, solver comparison |
| [visualization/](visualization/) | 9 | PyVista, ParaView, Netgen GUI, GMSH workflow demos |

### Solver & Performance

| Directory | Scripts | Description |
|-----------|---------|-------------|
| [solver_benchmarks/](solver_benchmarks/) | 15 | LU vs BiCGSTAB performance, scaling studies |
| [tetra_field_accuracy_evaluation/](tetra_field_accuracy_evaluation/) | 5 | Tetrahedron field accuracy vs analytical/NGSolve reference |

### NGSolve Integration

| Directory | Scripts | Description |
|-----------|---------|-------------|
| [ngsolve_integration/](ngsolve_integration/) | 14 | RadiaField CoefficientFunction, field types, mesh convergence |
| [KelvinTransformation/](KelvinTransformation/) | 96 | Kelvin transformation for unbounded domains (H/A formulation, adaptive mesh) |
| [beam_tracking/](beam_tracking/) | 2 | Particle trajectory in magnetic fields |

### PEEC Conductor Modeling

| Directory | Scripts | Description |
|-----------|---------|-------------|
| [peec_integration/](peec_integration/) | 98 | PEEC Loop-Star solver: coils, SPICE export, ngbem coupling, WPT |
| [Effective_Surface_Impedance/](Effective_Surface_Impedance/) | 5 | ESIM conductor model (Dowell + nonlinear homogenization) |
| [induction_heating/](induction_heating/) | 6 | ESIM induction heating, RWG-EFIE 3D, WPT coupling |
| [ngbem_diagnostics/](ngbem_diagnostics/) | 14 | Eddy current solver validation and diagnostics |

### Machine Learning

| Directory | Scripts | Description |
|-----------|---------|-------------|
| [Universal_Relaxation_Network/](Universal_Relaxation_Network/) | 24 | KAN-inspired URN for BH curves, SPICE time-domain |

---

## Example Selection Guide

| Use Case | Recommended |
|----------|-------------|
| Learn Radia basics | `simple_problems/` |
| Permanent magnets | `simple_problems/`, `smco_magnet_array/` |
| Electromagnets | `c_type_electromagnet/`, `complex_coil_geometry/` |
| External/background fields | `background_fields/` |
| FEM coupling | `ngsolve_integration/`, `KelvinTransformation/` |
| PEEC conductors | `peec_integration/`, `induction_heating/` |
| Eddy currents / shielding | `ngbem_diagnostics/`, `Effective_Surface_Impedance/` |
| Solver performance | `solver_benchmarks/` |
| Visualization | `visualization/`, `vtk_export/` |

---

## Common Patterns

### Unit Convention

- **Radia**: meters (m) — Radia always uses meters
- **NGSolve**: meters (m) — automatic conversion via `rad.RadiaField()`

### VTS Export

```python
```

### Material API

```python
rad.MatLin(1000)                                    # Linear (mu_r)
rad.MatPM(1.2, 900000, [0, 0, 1])                  # Permanent magnet (Br, Hc, dir)
rad.MatSatIsoTab([[0, 0], [100, 0.1], [1000, 1.2]])  # Nonlinear BH curve
```

---

## Prerequisites

**Required:** Python 3.12, Radia (`Build.ps1`), NumPy

**Optional:** NGSolve (Kelvin, ngbem, ESIM examples), ParaView (VTS viewing), SciPy, Matplotlib

---

## Troubleshooting

| Error | Solution |
|-------|----------|
| `ModuleNotFoundError: radia` | Build Radia: `powershell Build.ps1` |
| `ImportError: DLL load failed` | Install Visual C++ 2022 Redistributable |
| `ModuleNotFoundError: radia` (RadiaField) | Build Radia: `powershell Build.ps1` (RadiaField is now part of the main radia module since v2.5.0) |

---

**Last Updated:** 2026-02-22
