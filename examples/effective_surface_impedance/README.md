# Effective Surface Impedance (ESIM)

Examples and utilities for computing effective surface impedance boundary conditions (SIBC) for electromagnetic field simulations, supporting both linear and nonlinear (BH-curve) materials with applications to eddy-current loss and reactive power evaluation in conducting regions.

## Scripts

| File | Description |
|------|-------------|
| `quad_mesh_functions3D.py` | Utility functions for generating structured 3D hexahedral (cuboid) meshes with configurable macro-elements, domain numbering, and boundary conditions using Netgen |
| `surface_impedence_ref_Kengo.py` | Reference (volume-resolved) FEM solver that computes eddy-current losses and reactive power in a conducting cube by fully meshing the iron region with an HCurl/H1 mixed formulation, then compares results against effective P' and Q' curves |
| `surface_impedence_eff_Kengo.py` | Effective surface impedance FEM solver that replaces the volume iron mesh with an impedance boundary condition on the conductor surface, solving only the exterior H1 scalar potential problem to compute losses and reactive power |
| `esim_conductor_model.py` | ESIM-based conductor model for PEEC that combines DC bulk resistance with frequency-dependent skin-effect impedance via Dowell's analytical formula (linear) or ESIM cell-problem solver (nonlinear) |
| `esim_correct_implementation.py` | Correct ESIM implementation using a homogenization approach: solves the 1D cell problem for H(z), computes an effective permeability via |H|^2-weighted averaging, then applies Dowell's formula with the effective skin-depth parameter |

## Dependencies

- `ngsolve` and `netgen` -- finite element library and meshing (used by the reference and effective surface impedance solvers)
- `mylibcem` -- custom CEM library providing `BiotSavartCylinder`, `MuNonLinBiro3D`, and related utilities
- `Compumag2025.SIBC.stray_meshes` -- mesh generation routines for the SIBC study (provides `SICube` and related functions)
- `numpy` and `scipy` -- numerical computation and sparse linear algebra (used by the conductor model and homogenization solver)
- `esim_cell_problem` (from `radia`) -- provides `ESIMFiniteSlabSolver` and `BHCurveInterpolator` for the 1D cell-problem solve

## References

- K. Muramatsu, T. Nakata, N. Takahashi, K. Fujiwara, "Comparison of effective surface impedance boundary conditions", IEEE Transactions on Magnetics, 1994
- H. Igarashi, "Effective surface impedance for homogenization of nonlinear conducting regions", applicable to ESIM cell-problem approach
- Dowell, P.L., "Effects of eddy currents in transformer windings", Proceedings of the IEE, 1966 -- analytical R_ac/R_dc formula used for linear skin-effect computation
- Preis, K. and Biro, O., "Nonlinear surface impedance boundary conditions", used for the Preis-Biro iteration scheme in the reference solver
