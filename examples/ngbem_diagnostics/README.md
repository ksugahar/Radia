# ngbem Diagnostics

Validation and diagnostic suite for the eddy current solvers (ShieldBEMSIBC, VectorEddyCurrentFEMBEM, EddyCurrentFEMBEM, EddyCurrentBEMSIBC) used in the Radia project, covering loop-basis construction, BEM operator assembly, loss computation, and cross-method consistency checks.

## Scripts

| File | Description |
|------|-------------|
| `validate_shield_bem.py` | Validates ShieldBEMSIBC solver: loop basis sign sums, V_LL positive-definiteness, loss positivity, and frequency-sweep monotonicity for a wire above an aluminum plate. |
| `validate_vector_fembem.py` | Validates VectorEddyCurrentFEMBEM solver: FEM/BEM matrix assembly, system conditioning (Weggler stabilization), loss computation, mu_r support, frequency sweep, and comparison with the scalar Hz solver. |
| `test_shield_bem_mu_r.py` | Tests ShieldBEMSIBC with mu_r > 1 (aluminum, steel, soft iron), verifying that skin depth, surface impedance, and loss scale correctly with relative permeability. |
| `validate_shield_vs_vector.py` | Cross-validates ShieldBEMSIBC against VectorEddyCurrentFEMBEM on identical meshes for aluminum and steel blocks, including a mesh-convergence study. |
| `diagnose_vector_fembem.py` | Diagnoses frequency-independent loss in VectorEddyCurrentFEMBEM by dumping internal variables, checking for a missing curl-curl RHS term, and analyzing shielding cancellation. |
| `test_sibc_loss.py` | Compares three loss computation methods (volume integral, analytical SIBC on boundary faces, ShieldBEMSIBC) against the half-space analytical reference. |
| `compare_eddy_methods.py` | Runs VectorFEMBEM, ShieldBEMSIBC, and ScalarFEM solvers in isolated subprocesses and produces a cross-method comparison table with analytical references. |
| `debug_fembem_singular.py` | Debugs singularity in the VectorEddyCurrentFEMBEM system matrix by inspecting SVDs and ranks of FEM, BEM, and coupling sub-blocks, and checking FreeDofs and null-space structure. |
| `diagnose_basis_eval.py` | Evaluates RT0 basis functions directly from NGSolve GridFunctions, compares with manual formulas, verifies mass-matrix diagonals, and builds a corrected divergence matrix with proper orientation. |
| `diagnose_eddy_methods.py` | Runs all four eddy current solvers at a single frequency with detailed intermediate output in subprocesses to pinpoint root causes of inter-method inconsistency. |
| `diagnose_mfull.py` | Checks eigenvalues and positive-definiteness of the full mass matrix (M_full) and Maxwell SLP operator (V_full), both in the edge space and projected into the loop subspace. |
| `diagnose_vvec_normalization.py` | Determines the normalization convention of NGSolve's MaxwellSingleLayerPotentialOperator by manually computing V_vec and V_div double integrals and extracting components via kappa variation. |
| `systematic_eddy_validation.py` | Sweeps geometry (thin plate vs. block) and frequency to map out validity regions for each solver method, with a detailed analysis of why scalar Hz formulations underpredict 3D losses. |
| `validate_shielded_peec.py` | Validates the ShieldedPEECSolver (PEEC + ShieldBEMSIBC coupling): no-shield baseline, wire-above-plate impedance modification, FastHenry `.shield` block parsing, and A_scattered field evaluation. |

## Dependencies

- NumPy
- SciPy (used in `diagnose_basis_eval.py` for null-space computation)
- NGSolve (finite element library: `ngsolve`, `netgen.occ`)
- ngsolve.bem (boundary element module: `MaxwellSingleLayerPotentialOperator`, `HelmholtzSL`)
- `ngbem_eddy` module from `src/radia` (ShieldBEMSIBC, VectorEddyCurrentFEMBEM, EddyCurrentFEMBEM, EddyCurrentBEMSIBC, LoopBasisBuilder)
- `peec_matrices`, `peec_topology`, `peec_shielded`, `fasthenry_parser` modules from `src/radia` (used by `validate_shielded_peec.py`)
