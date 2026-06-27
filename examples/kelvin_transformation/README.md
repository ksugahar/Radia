# Kelvin Transformation

Finite element examples for unbounded magnetostatic problems using the **Kelvin transformation**. The Kelvin transformation maps an infinite exterior domain (r > R) to a finite computational domain by the inversion r' = R^2/r, enabling standard FEM discretisation on a bounded mesh while exactly representing the far-field decay. All examples are implemented with [NGSolve](https://ngsolve.org/) / Netgen.

## Canonical convention

All examples and the centralised `radia.kelvin_source` API use the single canonical convention derived in:

> **H. Nagamine, T. Yamaguchi, K. Sugahara**, "A Pullback-Based Formulation of Kelvin Transformation in Electromagnetic Field Analysis," CEFC 2026 (Thessaloniki), id 350.

For 3D spherical (conformal) Kelvin:

```
nu_ext = (rho'/R)^2 * nu_0          [HCurl A-formulation]
mu_ext = (R/rho')^2 * mu_0          [H1 Omega / H-formulation]
```

These are pointwise reciprocals (mu * nu = 1), consistent with Kelvin as a *physical* coordinate transformation. See [CONVENTION.md](CONVENTION.md) for the one-page canonical convention declaration.

**Comprehensive theory** (consolidated 2026-05-04):
[`docs/kelvin/KELVIN_TRANSFORMATION.md`](../../docs/kelvin/KELVIN_TRANSFORMATION.md) — single master doc covering:
- §2: 1-form / 2-form pullback derivations (Convention A)
- §3-5: Geometry (Sugahara two-sphere) and FEM workflow
- §7: Reduced potential formulations + Kelvin (H- and A-formulations, including the **(ν - ν₀) form pitfall** — critical for PEEC + Kelvin)
- §8: API reference, §9: Cubit workflow, §10: usage examples
- §11: Known limitations
- §12: References

This README only describes example problems; formulation theory is not duplicated here.

**Quickstart (student-friendly API)**: Import factors from `radia.kelvin_source` rather than re-deriving inline:

```python
from radia.kelvin_source import (
    kelvin_nu_factor_axisym_cf,   # (rho'/R)^2 for A-form, axisym with Z-offset
    kelvin_mu_factor_3d_cf,       # (R/rho')^2 for Omega-form, 3D sphere
    build_material_cf,            # {material: value} CF builder
)

# A-formulation, axisym with Z-offset:
nu_cf = build_material_cf(
    mesh, nu0,
    kelvin_nu_factor_axisym_cf(z_offset, R_kelvin),
    overrides={"magnetic": nu0 / mu_r},
)
```

If a particular Kelvin setup shows unexpected inductance / field / energy errors, **do not change the convention**. Debug the FEM setup separately: GND placement (Dirichlet at rho'=0), gauge regularisation, mesh refinement near rho'=0 (material vanishes or diverges), integration order (bonus_intorder for rational coefficients), and periodic identification.

## Verifying Kelvin Periodic identification (FES-only, sub-second)

Before running any expensive physics solve, you can verify in <1 s
that the Kelvin Periodic BC is actually being enforced.  Two
complementary checks (no `BilinearForm.Assemble`, no solver call):

1. **Slaved-DOF count** -- Periodic FES should eliminate slave DOFs:

   ```python
   from ngsolve import H1, Periodic
   fes_plain    = H1(mesh, order=p, dirichlet="GND")
   fes_periodic = Periodic(H1(mesh, order=p, dirichlet="GND"))
   slaved = sum(fes_plain.FreeDofs())  - sum(fes_periodic.FreeDofs())
   # slaved > 0 means the Periodic constraint is in place.
   # At p=2, slaved should be roughly 4x the p=1 value because
   # NGSolve adds high-order edge / face DOF coupling automatically.
   ```

2. **Boundary-integral ratio** -- set 1.0 on the slave boundary, integrate
   the same field on the master boundary; if Periodic is wired correctly
   the ratio is 1.0:

   ```python
   from ngsolve import GridFunction, Integrate
   gfu = GridFunction(fes_periodic);  gfu.vec[:] = 0
   gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
   a_int = float(Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_int")))
   a_ext = float(Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_ext")))
   assert abs(a_ext / a_int - 1.0) < 1e-3,  "Periodic Kelvin BC not enforced"
   ```

If either check fails, the issue is in `Identify()` placement (must be
*after* `Glue` for OCC) or the Periodic-FES setup -- not in the
formulation.  Do NOT run a multi-minute solve to discover this.

## Subdirectories

### H-formulation

Scalar-potential (H-field) perturbation formulation for magnetostatics.
Dipole and quadrupole problems in 2D, 3D, and axisymmetric geometries, with and without Kelvin transformation.

The standalone H-formulation demo scripts were promoted to
[`docs/kelvin/kelvin_classic_demos.ipynb`](../../docs/kelvin/kelvin_classic_demos.ipynb).
Their full source text and SHA-256 hashes are preserved in
[`docs/kelvin/kelvin_classic_demos_results.json`](../../docs/kelvin/kelvin_classic_demos_results.json).

### A-formulation

Vector-potential (A-field) formulation for magnetostatics, including axisymmetric coil problems and the z-offset Kelvin transformation variant.

See [A-formulation/README.md](A-formulation/README.md) for details.

The classic standalone A-formulation demos were promoted to
[`docs/kelvin/kelvin_classic_demos.ipynb`](../../docs/kelvin/kelvin_classic_demos.ipynb)
with full source preserved in
[`docs/kelvin/kelvin_classic_demos_results.json`](../../docs/kelvin/kelvin_classic_demos_results.json).
Additional validation-named A-formulation scripts that were not part of
the classic-demo archive remain here until they move to `validation_test`
or a source API.

### Cubit_1_4_p_convergence (**VERIFIED 2026-04-25**)

End-to-end p-convergence demonstration for the Cubit-meshed
`export netgen` -> NGSolve Omega-Reduced Omega + Kelvin pipeline.
Cubit-builds a 1/4 sector (x>=0, y>=0, full z) mu_r=100 sphere in
uniform Hz, exports `.vol` at p=1, 2, 3, solves and probes Hz at
origin.  At p=2 matches analytical to **+0.71%**.  See
[Cubit_1_4_p_convergence/README.md](Cubit_1_4_p_convergence/README.md)
for full results and the two non-obvious Cubit fixes documented there.

**Promoted to all three layers**:
- **tests/** : `tests/cubit/test_kelvin_1_4_p_convergence.py` (mesh
  + solver regression, 2 cases) + `tests/panels/test_kelvin_benchmark_golden.py`
  (full panel chain golden, 4 cases) -- 17 tests pass total when
  combined with `test_panel_qa.py` (9) and `test_kelvin_periodic_fes.py` (4).
- **examples/** : this directory (canonical sample with full README).
- **panels/** : `radia_em.py` ships a "Kelvin Benchmark" formulation
  (4th option in the Formulation combo, alongside Omega / A-Phi /
  MSC) that runs `calc_kelvin_benchmark.py` on the bundled
  `panels/samples/kelvin_benchmark_sphere_1_4.vol`.  `radia-mcp`
  documents this via `kelvin_knowledge.get_kelvin_documentation("benchmark_panel")`.

### Omega_ReducedOmega

Omega-Reduced Omega method (total/reduced scalar potential) with Kelvin transformation for 3D and axisymmetric magnetostatics. Includes sphere and cylinder benchmark geometries.

The classic standalone Omega-ReducedOmega demos were promoted to
[`docs/kelvin/kelvin_classic_demos.ipynb`](../../docs/kelvin/kelvin_classic_demos.ipynb),
including the former sphere and cylinder OCC reference implementations.
Full source text and source hashes are preserved in
[`docs/kelvin/kelvin_classic_demos_results.json`](../../docs/kelvin/kelvin_classic_demos_results.json).
Validation and API-candidate scripts under this directory remain until
their behavior is locked in `validation_test` or lifted into `src`.

### DtN_spectrum

The standalone DtN/open-boundary act scripts were promoted to
[`docs/kelvin/kelvin_dtn_spectrum_archive.ipynb`](../../docs/kelvin/kelvin_dtn_spectrum_archive.ipynb)
and pruned from `examples/`.  Full source text and SHA-256 hashes for
all 122 scripts are preserved in
[`docs/kelvin/kelvin_dtn_spectrum_archive_results.json`](../../docs/kelvin/kelvin_dtn_spectrum_archive_results.json).
Productionized behavior lives in `src/radia/open_boundary` and
`validation_test/open_boundary`.

### AdaptiveMesh

Adaptive mesh refinement studies combining Kelvin transformation with ZZ error estimation, Doerfler marking, metric-based remeshing, and Laplacian smoothing.

#### AdaptiveMesh/Omega_ReducedOmega

Convergence studies for the Omega method on sphere, cylinder, cube, and axisymmetric sphere geometries. Each geometry contains subdirectories for polynomial orders (p=2,3,4) and refinement strategies (uniform, ZZ-adaptive, metric-based, maxh). See [AdaptiveMesh/Omega_ReducedOmega/README.md](AdaptiveMesh/Omega_ReducedOmega/README.md) for full results.

| File | Description |
|------|-------------|
| `Sphere_3D/compare_convergence.py` | Convergence comparison plot for 3D sphere across orders and methods |
| `Cylinder_3D/compare_convergence.py` | Convergence comparison plot for 3D cylinder |
| `Cube_3D/compare_convergence.py` | Convergence comparison plot for 3D cube |
| `Cube_3D/CubeMesh.py` | Cube mesh generator class for the cube benchmark |
| `Sphere_Axisymmetric/compare_convergence.py` | Convergence comparison plot for axisymmetric sphere |

#### AdaptiveMesh/TEAM7

TEAM Problem 7 benchmark: 3D eddy current analysis of an asymmetric aluminum plate with a hole, excited by a racetrack coil. See [AdaptiveMesh/TEAM7/README.md](AdaptiveMesh/TEAM7/README.md).

| File | Description |
|------|-------------|
| `team7_geometry.py` | Geometry creation and material properties for TEAM 7 |
| `team7_solver.py` | A-method solver for time-harmonic eddy current analysis |
| `team7_coil_current.py` | Racetrack coil current density definition as CoefficientFunction |
| `team7_A_method.py` | Main driver: 3D eddy current analysis using A-phi formulation |
| `test_weighted_average.py` | Weighted-average convergence study for A-method and Omega-method |

#### AdaptiveMesh/A-formulation

| File | Description |
|------|-------------|
| `CircularCoil_A_formulation_with_Kelvin.py` | 3D A-formulation for a circular coil with Kelvin transformation (1/8 model) |

#### Other AdaptiveMesh examples

| Directory / File | Description |
|------------------|-------------|
| `adaptive_mesh_with_smoothing.py` | 2D adaptive mesh refinement with ZZ estimator and Laplacian smoothing on an L-shaped domain |
| `Coil_A_formulation_adaptive.py` | Axisymmetric coil A-formulation with adaptive refinement and z-offset Kelvin transformation |
| `compare_convergence.py` (parallel wires) | Convergence comparison for 2D parallel-wire adaptive refinement across orders and methods |

### docs

Reference documentation on the mathematical foundations and implementation details.

| File | Description |
|------|-------------|
| `Kelvin_2D.md` | Mathematical derivation of the 2D Kelvin transformation for H-formulation |
| `Kelvin_3D.md` | Mathematical derivation of the 3D Kelvin transformation for H-formulation |
| `kelvin_transform_cylinder.md` | Kelvin transformation in cylindrical coordinates (r and z directions) |
| `Supplement/ErrorEstimator.md` | Equilibrated error estimator theory for edge-element FEM |
| `Supplement/CG-smoother.md` | CG-smoother acceleration for equilibrated error estimation |
| `Supplement/test_cg_smoother_equilibration.py` | Test script: CG-smoother vs direct solve for equilibrated error estimation |

## Dependencies

- **[NGSolve](https://ngsolve.org/)** / **Netgen** -- finite element library and mesh generator (required)
- **NumPy** -- array operations
- **SciPy** -- special functions (`ellipk`, `ellipe`), `.mat` file I/O, sparse linear algebra
- **Matplotlib** -- convergence plots and field visualisation
