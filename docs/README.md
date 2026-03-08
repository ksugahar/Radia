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
| [PEEC_PANEL_IMPLEMENTATION.md](PEEC_PANEL_IMPLEMENTATION.md) | Panel integration (Wilton, Gauss) + FastImp Loop-Star integration |
| [PEEC_SURFACE_IMPEDANCE.md](PEEC_SURFACE_IMPEDANCE.md) | Surface impedance (SIBC, ESIM, Dowell, PyKAN, SPICE export) |
| [PEEC_SHIELD_CONDUCTOR.md](PEEC_SHIELD_CONDUCTOR.md) | Shield conductor modeling |
| [PEEC_MSC_COUPLING.md](PEEC_MSC_COUPLING.md) | PEEC-MSC coupled solver theory |
| [PEEC_VALIDATION_PLAN.md](PEEC_VALIDATION_PLAN.md) | Validation phases + formula validation + Grover formula |

## Mesh & Visualization

| Document | Description |
|----------|-------------|
| [MESH_GUIDE.md](MESH_GUIDE.md) | Mesh types, GMSH workflows, SetGeomInfo API, surface elements |
| [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md) | Viewer comparison, workflows, geometry accuracy |

## ELF Compatibility & Verification

| Document | Description |
|----------|-------------|
| [ELF_COMPATIBILITY_GUIDE.md](ELF_COMPATIBILITY_GUIDE.md) | Unit systems, mesh conventions, matrix formulation |

## Advanced Solvers & Methods

| Document | Description |
|----------|-------------|
| [EDDY_CURRENT_METHODS.md](EDDY_CURRENT_METHODS.md) | Eddy current methods + VectorFEMBEM cross-validation |
| [NGBEM_INTEGRATION_DESIGN.md](NGBEM_INTEGRATION_DESIGN.md) | NGSolve BEM integration + FEM independent verification |
| [UNIFIED_FIELD_API_DESIGN.md](UNIFIED_FIELD_API_DESIGN.md) | Unified field computation API |

## Hysteresis & Nonlinear Materials

| Document | Description |
|----------|-------------|
| [B_INPUT_PLAY_MODEL.md](B_INPUT_PLAY_MODEL.md) | B-input Play hysteresis model (direct Play + energy-based approximation) |

## Specialized Features

| Document | Description |
|----------|-------------|
| [BEAM_TRACKING_DESIGN.md](BEAM_TRACKING_DESIGN.md) | Beam tracking in magnetic fields |
| [CYLINDRICAL_MAGNET.md](CYLINDRICAL_MAGNET.md) | Cylindrical magnet analytical field |
| [IMA_SYMMETRY_DESIGN.md](IMA_SYMMETRY_DESIGN.md) | Image symmetry for MSC |
| [KAN_INSPIRED_URN.md](KAN_INSPIRED_URN.md) | KAN-inspired Universal Relaxation Network |
| [NPORT_BLOCK_LANCZOS_SPICE.md](NPORT_BLOCK_LANCZOS_SPICE.md) | N-port Block Lanczos SPICE + CLN I verification |

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

**Last Updated:** 2026-03-05
