# Kelvin Transformation

Finite element examples for unbounded magnetostatic problems using the **Kelvin transformation**. The Kelvin transformation maps an infinite exterior domain (r > R) to a finite computational domain by the inversion r' = R^2/r, enabling standard FEM discretisation on a bounded mesh while exactly representing the far-field decay. All examples are implemented with [NGSolve](https://ngsolve.org/) / Netgen.

## Subdirectories

### H-formulation

Scalar-potential (H-field) perturbation formulation for magnetostatics.
Dipole and quadrupole problems in 2D, 3D, and axisymmetric geometries, with and without Kelvin transformation.

| File | Description |
|------|-------------|
| `2D_dipole.py` | 2D dipole H-formulation on a finite circular domain (no Kelvin) |
| `2D_dipole_with_Kelvin.py` | 2D dipole H-formulation with Kelvin transformation and periodic BC |
| `2D_dipole_half_with_Kelvin.py` | Half-circle variant to test periodic BC with segmented edges |
| `2D_quadrupole.py` | 2D quadrupole H-formulation on a finite domain (no Kelvin) |
| `2D_quadrupole_with_Kelvin.py` | 2D quadrupole H-formulation with Kelvin transformation |
| `3D_dipole.py` | 3D dipole H-formulation for a magnetic sphere in uniform field |
| `3D_dipole_with_Kelvin.py` | 3D dipole H-formulation with Kelvin transformation and periodic BC |
| `3D_quadrupole.py` | 3D quadrupole H-formulation with background field H_s = (-z, 0, -x) |
| `3D_quadrupole_with_Kelvin.py` | 3D quadrupole with Kelvin transformation and periodic BC |
| `Axisymmetric_dipole.py` | Axisymmetric H-formulation for a magnetic sphere (no Kelvin) |
| `Axisymmetric_dipole_with_Kelvin.py` | Axisymmetric H-formulation with Kelvin transformation |
| `Laplace3D_dipole_with_Kelvin.py` | 3D Laplace equation with Kelvin transformation (conventional, non-perturbation) |
| `Fig_1.py` | 2x2 panel contour/streamline plot for 3D dipole results |
| `Fig_2.py` | 2x2 panel contour/streamline plot for 3D quadrupole results |
| `Fig_3.py` | Line plot comparing Hz along the x-axis for dipole and quadrupole |

### A-formulation

Vector-potential (A-field) formulation for magnetostatics, including axisymmetric coil problems and the z-offset Kelvin transformation variant.

See [A-formulation/README.md](A-formulation/README.md) for details.

| File | Description |
|------|-------------|
| `Coil_A_formulation_simple.py` | Axisymmetric A-formulation for a single coil (no Kelvin, baseline) |
| `Coil_A_formulation_with_Kelvin.py` | Axisymmetric A-formulation with z-offset Kelvin transformation |
| `A_formulation_sphere_simple.py` | Axisymmetric A-formulation for a magnetic sphere in uniform field (no Kelvin) |
| `A_formulation_sphere_with_Kelvin.py` | Axisymmetric A-formulation for a magnetic sphere with z-offset Kelvin transformation |
| `ParallelWires_2D_A_formulation_with_Kelvin.py` | 2D A-formulation with Kelvin transformation for parallel wires |
| `sphere_in_uniform_field.py` | Magnetic sphere in uniform field: compares A-method and Omega-method with equilibrated error estimation |

### Omega_ReducedOmega

Omega-Reduced Omega method (total/reduced scalar potential) with Kelvin transformation for 3D and axisymmetric magnetostatics. Includes sphere and cylinder benchmark geometries.

| File | Description |
|------|-------------|
| `Omega_ReducedOmega.py` | Base class implementing the Omega-Reduced Omega solver |
| `Sphere/3D_sphere_with_Kelvin.py` | 3D magnetic sphere (mu_r=100) solved with Omega-Reduced Omega and Kelvin transformation |
| `Sphere/Axisymmetric_sphere_with_Kelvin.py` | Axisymmetric sphere variant of the Omega-Reduced Omega solver |
| `Cylinder/3D_cylinder_with_Kelvin.py` | 3D magnetic cylinder with Kelvin transformation (full model) |
| `Cylinder/3D_cylinder_with_Kelvin_1_8.py` | 3D magnetic cylinder using 1/8 symmetry model |
| `Cylinder/Axisymmetric_cylinder_with_Kelvin.py` | Axisymmetric cylinder variant of the Omega-Reduced Omega solver |

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
