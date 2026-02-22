# Radia Documentation

**Version:** 2.3.0

## Core Documentation

| Document | Description |
|----------|-------------|
| [API_REFERENCE.md](API_REFERENCE.md) | **Main documentation** - Quick start, API, elements, solvers, NGSolve |
| [SOLVER_ARCHITECTURE.md](SOLVER_ARCHITECTURE.md) | Solver design philosophy and architecture |
| [APPLICATION_ROADMAP.md](APPLICATION_ROADMAP.md) | Research applications and collaboration roadmap |

## PEEC Solver

| Document | Description |
|----------|-------------|
| [PEEC_CONDUCTOR_MODELING_GUIDE.md](PEEC_CONDUCTOR_MODELING_GUIDE.md) | Conductor shape interface and mesh import workflow |
| [PEEC_PANEL_IMPLEMENTATION.md](PEEC_PANEL_IMPLEMENTATION.md) | 2D panel integration (Wilton formula, Gauss quadrature) |
| [PEEC_SURFACE_IMPEDANCE.md](PEEC_SURFACE_IMPEDANCE.md) | Surface impedance (SIBC, ESIM, Dowell, PyKAN, SPICE export) |
| [PEEC_SHIELD_CONDUCTOR.md](PEEC_SHIELD_CONDUCTOR.md) | Shield conductor modeling |
| [PEEC_MSC_COUPLING.md](PEEC_MSC_COUPLING.md) | PEEC-MSC coupled solver theory |
| [PEEC_FASTIMPINTEGRATION.md](PEEC_FASTIMPINTEGRATION.md) | FastImp Loop-Star integration |
| [PEEC_VALIDATION_PLAN.md](PEEC_VALIDATION_PLAN.md) | Systematic validation phases |
| [PEEC_FORMULA_VALIDATION.md](PEEC_FORMULA_VALIDATION.md) | Formula validation (Neumann integral) |

## Mesh & Visualization

| Document | Description |
|----------|-------------|
| [GMSH_WORKFLOW.md](GMSH_WORKFLOW.md) | CAD -> GMSH -> NGSolve -> Radia workflow |
| [SETGEOMINFO_API.md](SETGEOMINFO_API.md) | SetGeomInfo API for high-order curving (ksugahar/ngsolve fork) |
| [VISUALIZATION_WORKFLOW.md](VISUALIZATION_WORKFLOW.md) | Radia-NGSolve visualization workflow |
| [VIEWER_COMPARISON.md](VIEWER_COMPARISON.md) | Viewer comparison (PyVista, ParaView, GMSH, Netgen) |
| [MESH_AND_SURFACE_ELEMENTS_GUIDE.md](MESH_AND_SURFACE_ELEMENTS_GUIDE.md) | Surface elements FAQ and Netgen GUI guide |
| [GEOMETRY_ACCURACY_COMPARISON.md](GEOMETRY_ACCURACY_COMPARISON.md) | VTS export geometry accuracy |

## ELF Compatibility & Verification

| Document | Description |
|----------|-------------|
| [ELF_COMPATIBILITY_GUIDE.md](ELF_COMPATIBILITY_GUIDE.md) | Unit systems, mesh conventions, matrix formulation |

## Advanced Solvers & Methods

| Document | Description |
|----------|-------------|
| [EDDY_CURRENT_METHODS.md](EDDY_CURRENT_METHODS.md) | Eddy current method comparison (ngbem FEM-BEM, BEM+SIBC, FEM+Kelvin) |
| [GROVER_FORMULA_IMPLEMENTATION.md](GROVER_FORMULA_IMPLEMENTATION.md) | Grover's formula for self-inductance |
| [VECTOR_FEMBEM_ANALYSIS.md](VECTOR_FEMBEM_ANALYSIS.md) | VectorEddyCurrentFEMBEM analysis |
| [NGBEM_INTEGRATION_DESIGN.md](NGBEM_INTEGRATION_DESIGN.md) | NGSolve BEM integration design (low-freq verified) |
| [NGSolve_FEM_Verification.md](NGSolve_FEM_Verification.md) | FEM verification of PEEC+BEM inductance |
| [UNIFIED_FIELD_API_DESIGN.md](UNIFIED_FIELD_API_DESIGN.md) | Unified field computation API |
| [CLN_I_VERIFICATION.md](CLN_I_VERIFICATION.md) | CLN I type coordinate transform verification |

## Specialized Features

| Document | Description |
|----------|-------------|
| [BEAM_TRACKING_DESIGN.md](BEAM_TRACKING_DESIGN.md) | Beam tracking in magnetic fields |
| [CYLINDRICAL_MAGNET.md](CYLINDRICAL_MAGNET.md) | Cylindrical magnet analytical field |
| [IMA_SYMMETRY_DESIGN.md](IMA_SYMMETRY_DESIGN.md) | Image symmetry for MSC |
| [KAN_INSPIRED_URN.md](KAN_INSPIRED_URN.md) | KAN-inspired Universal Relaxation Network |
| [NPORT_BLOCK_LANCZOS_SPICE.md](NPORT_BLOCK_LANCZOS_SPICE.md) | N-port Block Lanczos SPICE generation |

## Quick Links

- [Main README](../README.md) - Project overview and installation
- [CHANGELOG](../CHANGELOG.md) - Version history
- [Examples](../examples/) - Working code examples

## References

### MSC (Magnetic Surface Charge) Method

- **Newell, A. J., Williams, W., & Dunlop, D. J.** (1993).
  "A generalization of the demagnetizing tensor for nonuniform magnetization."
  *Journal of Geophysical Research: Solid Earth*, 98(B6), 9551-9555.
  DOI: [10.1029/93JB00694](https://doi.org/10.1029/93JB00694)

- **Yano, T. & Sugahara, K.** (2023).
  "MMM with the Idea of Magnetic Surface Charge Method."
  *Journal of the Magnetics Society of Japan*, 47(4), 89-94.

### ESIM (Effective Surface Impedance Method)

- **Hollaus, K., Kaltenbacher, M., & Schoberl, J.** (2025).
  "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation."
  *IEEE Transactions on Magnetics*.
  DOI: [10.1109/TMAG.2025.3613932](https://doi.org/10.1109/TMAG.2025.3613932)

---

**Last Updated:** 2026-02-22
