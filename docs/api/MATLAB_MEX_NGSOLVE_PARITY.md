# MATLAB MEX and NGSolve Boundary

This document defines the MATLAB boundary for the Radia and NGSolve-based
solver stack. The MEX gateway does not embed or launch Python. The current
pip-provided Windows `libngsolve.dll` nevertheless has a transitive
`python312.dll` dependency, so this is Python-process-free but not yet
Python-DLL-free. A fully standalone deployment requires an NGSolve/Netgen
build without Python support.

## Radia parity

The public top-level functions, underscore-prefixed numerical kernels, and
pybind11 class members are compared against the command inventory in
`src/matlab/radia_mex.cpp` during the MCP contract test. Mesh-plumbing helpers
and Python-only regression helpers are classified explicitly instead of being
silently dropped from the count.

| Surface | Current contract |
|---|---:|
| Public top-level pybind11 names | 100 |
| Covered mapped MEX names | 100 |
| Underscore numerical kernels | 21 / 21 covered |
| Stateful pybind11 class surface | 126 / 126 covered |
| Radia MEX gateway commands | 365 |
| Optuna MEX gateway commands | 21 |
| MATLAB Optuna classes | 95 |
| MATLAB Optuna functions | 90 |

The gateways are intentionally independent. `radia_mex` contains no
`optuna.*` commands and owns the Radia/NGSolve/oneMKL dependency graph.
`optuna_mex` contains two API commands and 19 optimizer commands, has no
NGSolve, oneMKL, Radia-core, or Python dependency, and is loaded only when
optimization is used. A missing or non-21-command `optuna_mex` fails loudly; it
is never replaced by a call through `radia_mex` or a silent MATLAB substitute.

Acoustic sphere and convolution-quadrature reference formulas deliberately do
not count as native parity. They remain readable, independent SciPy and MATLAB
implementations under `src/radia/acoustics` and `matlab/+radia/+acoustic` so
they can validate the NGSolve/ngsolve.bem CQ-BEM and FSI routes without sharing
a C++, pybind11, or MEX kernel with the solver under test.

`radia.ObjTetrahedronCurrent` is part of the standalone MEX surface and is
tested through the same object/field API as its pybind11 counterpart. The four
`KelvinRadia*` classes are intentionally different: they are NGSolve
`CoefficientFunction` adapters rather than independent numerical kernels.
MATLAB constructs that batch-only NGSolve composition through
`radia.python.openBoundary`; their bindings and constructors are explicitly
classified in the parity inventory instead of being reported as missing MEX
commands.

The `axifem.q1_magnetic_element_matrices` command shares the Q1 Henrotte
stiffness and sigma-mass implementation with `radia.axifem` and the production
NGSolve bilinear-form integrators. MATLAB receives numeric 4-by-4 element
matrices; NGSolve retains ownership of finite-element spaces and assembly.

The 22 HLU names use the MATLAB MEX snake-case namespace. Three legacy object
constructors have been removed from both public Python and MATLAB APIs because
the underlying implementation is process-unsafe or no longer available:

- `radia.ObjMltExtPgn`
- `radia.ObjMltExtRtg`
- `radia.ObjMltExtTri`

Their old C ABI declarations, definitions, and exports are deleted as well; no
compatibility shims remain. Use Netgen or Cubit geometry and the NGSolve-native
mesh path instead.

## NGSolve boundary

NGSolve is an external finite-element system and its complete Python API is not
duplicated as MATLAB MEX functions. NGSolve remains the owner of:

- mesh loading and curved geometry
- HCurl/HDiv element orientation and Piola mappings
- finite-element spaces and local-to-global DoF transforms
- arbitrary Python-side object identity and callback objects
- TaskManager-owned finite-element assembly

MATLAB receives numeric or value-type contracts at the boundary:

- `ngsolve.space_info` for H1, HCurl, and HDiv space dimensions
- `ngsolve.matrix_dump` for assembled sparse matrices and metadata
- `ngsolve.coefficient_function.*` for native coefficient-expression handles
  and mapped physical-point evaluation
- `ngsolve.grid_function.*` for native GridFunction handles, DoF vectors,
  interpolation, and GridFunction-to-CoefficientFunction views
- `ngsolve.grid_function.vector_handle` and `ngsolve.vector.*` for live
  GridFunction component views, independent work-vector copies, native
  `set_zero`, `scale`, `axpy`, `dot`, `norm`, and explicit value transfer
- `ngsolve.mesh.*`, `ngsolve.fespace.*`, and
  `ngsolve.bilinear_form.*` for persistent native mesh, finite-element-space,
  and assembled-form handles
- `ngsolve.mesh.set_deformation`, `ngsolve.mesh.unset_deformation`, and
  `ngsolve.mesh.trafo_quality` for native VectorH1 deformation and per-element
  GetTrafo Jacobian-ratio / spectral-condition sampling
- `ngsolve.matrix.*` for persistent native BaseMatrix handles, sparse triplet
  export, native matrix-vector products, free-DoF inverse operators, projected
  matrices, reduced block sums, and native diagonal preconditioners
- `ngsolve.radia_field.*` and `radia.RadiaField` for direct Radia B/H/A/M/Phi
  fields as NGSolve CoefficientFunctions, coordinate transforms, point caches,
  and native trilinear voxel CoefficientFunctions
- `radia.HACApKChargeGram.configureMassMatrixNGSolve` and
  `configureGeometryMassMatrixNGSolve` for direct registration of an assembled
  scalar NGSolve sparse matrix without a MATLAB sparse-array round trip
- `ngsolve.bilinear_form.create_from_coefficient` for scalar
  CoefficientFunction-weighted real/complex `mass`, `stiffness`, `curlcurl`,
  and `divdiv` integrators
- `ngsolve.bilinear_form.create_boundary_from_coefficient` for scalar
  CoefficientFunction-weighted real/complex trace bilinear forms
- `ngsolve.grid_function.from_fespace` for GridFunction construction that
  shares an existing persistent FESpace
- `ngsolve.linear_form.*` for persistent real/complex volume right-hand sides
  from constant sources or native CoefficientFunction expressions in H1,
  HCurl, and HDiv spaces, plus their native BaseVector views
- `ngsolve.linear_form.create_boundary_from_coefficient` for real/complex
  boundary right-hand sides using NGSolve trace mappings
- `hcurl.eddy_cln.native_basis` for direct C++ HCurl response-basis generation
- `radia.ngsolve.hcurl_eddy_cln_model` for a Python-free native HCurl
  response projection and local diffusion CLN state-space model
- `simulink.state_space.*` for persistent native discrete state-space handles
  used by Simulink without Python or per-step state-vector copies
- `hdiv.field_evaluator.*` for persistent tet, curved-tet, and source-cloud
  field evaluation
- `hdiv.planar_evaluator.*` for planar field and vector-potential evaluation
- `hdiv.field_evaluator.as_coefficient` and
  `hdiv.planar_evaluator.as_coefficient` for evaluator-backed native
  CoefficientFunctions with algorithm/rotation metadata
- `hacapk.charge_gram.*_directional_derivative` for complete TET, HEX, and
  WEDGE ChargeGram shape derivatives, including analytic self-panel blocks and
  TET charge-map row rates
- `hacapk.charge_gram.directional_derivative_contractions` for batched analytic
  `left.'*dG(k)*right` evaluation directly over parent H-matrix leaves
- `radia.HACApKChargeGramDerivative` for a persistent HACApK directional-
  derivative operator with entry, symmetric matvec, and compression statistics
- HCurl reduced Gram and Radia-owned mixed Galerkin/HDiv operators

The canonical MATLAB names mirror the snake_case MEX/Python names:

| MATLAB | MEX/Python contract |
|---|---|
| `radia.ngsolve.space_info` | `ngsolve.space_info` |
| `radia.ngsolve.matrix_dump` | `ngsolve.matrix_dump` |
| `radia.ngsolve.CoefficientFunction` | `ngsolve.coefficient_function.*` |
| `radia.ngsolve.GridFunction` | `ngsolve.grid_function.*` |
| `radia.ngsolve.Vector` | `ngsolve.grid_function.vector_handle` + `ngsolve.vector.*` |
| `radia.ngsolve.Mesh` | `ngsolve.mesh.*` |
| `radia.ngsolve.FESpace` | `ngsolve.fespace.*` |
| `radia.ngsolve.BilinearForm` | `ngsolve.bilinear_form.*` |
| `radia.ngsolve.Matrix` | `ngsolve.matrix.*` |
| `radia.ngsolve.LinearForm` | `ngsolve.linear_form.*` |
| `radia.RadiaField` | `ngsolve.radia_field.*` |
| `radia.HDivFieldEvaluator.coefficientFunction` | `hdiv.field_evaluator.as_coefficient` |
| `radia.PlanarFieldEvaluator.coefficientFunction` | `hdiv.planar_evaluator.as_coefficient` |
| `radia.HACApKChargeGramDerivative` | `hacapk.charge_gram_derivative.*` |
| `radia.ngsolve.hcurl_eddy_cln_native_basis` | `hcurl.eddy_cln.native_basis` |
| `radia.ngsolve.hcurl_eddy_cln_model` | `hcurl.eddy_cln.native_basis` + native CLN projection |
| `radia.hcurl.tet_reduced_gram` | `hcurl.tet_reduced_gram` |

Older names such as `radia.spaceInfo`, `radia.ngsolveMatrix`, and
`radia.tetHCurlReducedGram` remain supported as compatibility aliases.

MATLAB therefore uses NGSolve through Radia's numerically stable operations,
without pretending that a MATLAB array is an NGSolve object. This preserves
the NGSolve-native ownership policy while making the FEM workflow callable from
MATLAB and Simulink.

The native MEX gateway is linked against the NGSolve C++ libraries and oneMKL,
but not the Python runtime. The `radentry.cpp` C entry layer is
compiled with Python callback support disabled in this target; numeric Radia
operations remain available, while Python callback objects are rejected with a
checked boundary error. This keeps Python and MATLAB as two hosts of the same
C++ kernels rather than making MATLAB embed Python.

## Python numerical test parity

`tests/matlab/test_radia_ngsolve_parity.m` mirrors the important numerical
contracts from the Python Radia-NGSolve tests. The MATLAB tests call the C++
kernels directly through MEX; they do not launch Python or compare against a
Python process at runtime.

| Python numerical contract | MATLAB equivalent |
|---|---|
| `test_ngsolve_integration.py`: RadiaField API, field types, and metadata | `testRadiaFieldApiAndAllFieldTypes` |
| `test_ngsolve_integration.py`: HDiv projection and field accuracy | `testHDivGridFunctionProjectionAndFieldAccuracy` |
| `test_ngsolve_integration.py`: voxel CoefficientFunction | `testRadiaFieldVoxelCoefficient` |
| `test_ngsolve_integration.py`: H1, HCurl, and HDiv GridFunction creation | `testNGSolveSpaceGridFunctionCreation` |
| `test_rad_ngsolve.py`: RadiaField integration | `testRadiaFieldApiAndAllFieldTypes` |
| `test_radia_ngsolve_fields.py`: coil and magnet field relations | `testCurrentAndMagnetFieldRelations` |
| `test_radiafield_transformed_container.py`: transformed batch field | `testTransformedContainerBatchField` |
| `test_radiafield_transformed_container.py`: transformed field assembly | `testTransformedContainerLinearFormAssembly` |
| `test_radiafield_transformed_container.py`: transformed polyhedra and containers | `testTransformedPolyhedronAndContainerFields` |

`VectorH1` is available as `radia.ngsolve.FESpace.create(mesh,"vectorh1",p)`
because the Radia-owned sheet-topology method uses it for NGSolve-native mesh
deformation. It is assembled from NGSolve H1 component spaces; element maps,
DoF transforms, deformation ownership, and GetTrafo evaluation remain in
NGSolve. Python callback behavior and Python object identity remain
Python-only tests.

## Python/MEX NGSolve oracle

`tests/matlab/test_ngsolve_mex_pybind_parity.m` is the fast direct
cross-language gate. It generates a temporary MAT fixture through the public
NGSolve Python API on the same `.vol` file, then runs the MATLAB `radia_mex`
boundary against it. Expected numerical values never come from handwritten
MATLAB data.

The result-bearing campaign under
`validation_test/ngsolve_matlab_parity/` extends that gate without making
routine CI solver-heavy:

- 100 bounded linear cases split evenly between 2D and 3D;
- 500 TRI, QUAD, mixed TRI/QUAD, TET, HEX, and WEDGE cases;
- 20 matrix-free cases from 21,515 to 134,130 DoFs; and
- 15 manufactured-solution cases over five element families.

The checked HIBINO run passed all 535 extended cases on MATLAB R2026a Update
4 and NGSolve 6.2.2606, with zero live native handles at termination. MATLAB
may request a MATLAB `sparse` matrix at an explicit observation or optimization
boundary; `radia.ngsolve.Matrix`, `Vector`, and `Solver` retain NGSolve storage
for native repeated work between those boundaries.

The production sheet-topology driver is
`radia.topopt.optimizeHexSheetTopology`. Each outer batch performs 5--20
continuous-activation inner iterations through VectorH1/GetTrafo. Remove and
restore thresholds are distinct, and pending topology changes are committed
only when their fraction or age reaches the Cubit batch policy. A successful
Cubit rebuild is followed by exactly one H-matrix reconstruction; the inner
loop never rebuilds the H-matrix.

For Simulink, `simulink.state_space.create` stores `A`, `B`, `C`, `D`, and the
initial state in a checked native handle. `output` computes `y=C*x+D*u`
without mutation, while `update` advances `x=A*x+B*u` exactly once per sample.
`snapshot` and `restore` back the Level-2 S-Function `CustomSimState`; `step`
retains the atomic output-plus-update behavior for standalone diagnosis.
`reset`, `info`, and `destroy` complete the handle lifecycle.
`radia.simulink.buildHCurlEddyCLNModel(..., Block="radia-mex")` uses this path.
Induction heating does not: its distributed Eddy and Thermal fields are owned
by `radia_ih_eddy_sfun` and `radia_ih_thermal_sfun`, respectively. The former
IH LUT and lumped state-space builders have been removed.

The handle layer is explicit: MATLAB does not alias a Python object. A
`Mesh` owns the native `MeshAccess`, an `FESpace` shares that mesh, and a
`BilinearForm` retains its assembled native matrix. `Matrix` can export sparse
triplets, create native vectors, apply matvecs, and construct a free-DoF
inverse without rebuilding the mesh or form. A GridFunction vector handle is
a live native view, while `Vector.copy()` creates an independent native work
vector. Numeric DoF data is copied only when `values()` or `setValues()` is
requested. Expression arithmetic, matrix products, and vector iteration stay
in NGSolve, so large field updates do not round-trip through MATLAB at every
operation.

ChargeGram mass registration follows the same rule. Passing a
`radia.ngsolve.Matrix` extracts its scalar sparse rows inside C++, copies them
once into the persistent ChargeGram operator, and never materializes MATLAB
triplets. Affine HEX self-block shape derivatives also stay in the C++ kernel:
they use the same complete-cube outer rule as the affine-exact self block,
differentiate the Laplace kernel analytically, and retain the radial Duffy path
for non-affine Q2 geometry.

The persistent matrix slice is intentionally finite and explicit: it currently
assembles one built-in volume integrator per form (`mass`, `stiffness`,
`curlcurl`, or `divdiv`) and supports real and complex sparse matrices for
native matvec and inverse operations. Scalar CoefficientFunction weighting is
available for those integrators. Arbitrary Python callbacks, tensor-valued
bilinear forms, general preconditioner families, and arbitrary NGSolve solver
objects remain outside this MEX contract; callers must use the existing numeric
kernels or the Python/NGSolve host for those features. Native Krylov solver
handles (`cg`, `gmres`, and `bicgstab` for real matrices; `gmres` for complex
matrices) are supported. LinearForm coefficient expressions,
including complex volume and boundary forms, are supported through
`create_from_coefficient` and `create_boundary_from_coefficient`, retaining
NGSolve's native test-function and trace mappings.

`radia.ngsolve.hcurl_eddy_cln_model` is the first high-level native assembly
path. It calls the C++ HCurl response reduction and returns the exact reduced
projections `M_r = V' M V`, `K_r = V' K V`, and `P_r = V' ports`, then forms
`R = reluctivity*K_r` and `L = conductivity*M_r` for the local HCurl diffusion
convention. This gives MATLAB and Simulink a Python-free high-order FE/CLN
workflow. It is deliberately not advertised as the full HCurl-VIM path: the
external Laplace/BEM inductance, topology-aware eddy-bubble surface/bridge
classification, and frequency-dependent SIBC DtN rationalization remain
separate production stages.

## Optimization boundary

Continuous electromagnetic designs use a separate analytic-adjoint contract:

```matlab
evaluation = @(x) struct( ...
    "objective", objective(x), ...
    "gradient", objectiveAdjoint(x), ...
    "constraints", constraints(x), ...
    "constraint_jacobian", constraintAdjoints(x));

radia.topopt.checkAdjointGradient(x0, evaluation);
result = radia.topopt.optimizeAdjoint(x0, evaluation, ...
    Solver="mma", LowerBounds=xmin, UpperBounds=xmax);
```

Constraint values follow `c <= 0`; Jacobians have one row per design variable
and one column per constraint, matching the native HCurl MEX `case_gradient`
orientation and MATLAB's analytic-gradient solver interface. MMA requires
finite bounds and inequality constraints. SQP additionally accepts equality
constraints. Missing, transposed, nonfinite, or size-changing gradients fail
before an optimizer step. `checkAdjointGradient` uses centered directional
differences only as an explicit QA diagnostic, never as a production fallback.

`radia.topopt.optimizeHCurlActivationAdjoint` connects this contract directly
to `hcurl.topopt.activation_multifrequency_joule`: the MEX kernel performs the
complex state and adjoint solves, MATLAB adds the normalized material-volume
constraint, and MMA or SQP updates the bounded cell activation vector. TPE and
CMA-ES remain useful outer/global optimizers; MMA/SQP own the continuous inner
refinement when an analytic field sensitivity is available.

`radia.optuna` follows the familiar Optuna workflow:

```matlab
study = radia.optuna.create_study( ...
    direction="minimize", ...
    sampler=radia.optuna.TPESampler(Seed=42));
study.optimize(@(trial) objective(trial), 30);
best = study.best_trial();
bestParams = study.bestParams();
bestSnapshot = study.bestSolution();
```

The independently released `radia-optuna` package uses pinned Optuna 5.0.0
as its behavioral oracle. See the [package contract](../../packages/radia-optuna/README.md)
for the current native, Python fallback, persistence and Simulink boundaries.
The Radia gateway and its solver dependencies remain separate from the
21-command `optuna_mex` gateway.

Unified `TPESampler` supports single- and multi-objective studies. Its
`Multivariate=[]` setting represents the upstream automatic default and
`ConstantLiar=true` is the current default. MATLAB table/MAT persistence,
parallel scheduling and Simulink telemetry are extension contracts; their
presence does not establish upstream algorithmic parity.

Compatibility evidence is generated from the pinned version and checked in:

- `tests/matlab/fixtures/optuna50_oracle.json`: differential behavior fixtures;
- `tests/matlab/fixtures/optuna50_public_api.json`: upstream API inventory;
- `matlab/optuna50_api_coverage.json`: mapped, verified, asserted and unresolved coverage;
- `tests/matlab/fixtures/optuna_test_manifest.json`: test classification and oracle ownership.

Read the coverage record's `full_compatibility_complete` and its constructor
default audit before making a completeness claim. Earlier 4.9 timing and
proposal-sequence results do not certify the 5.0 implementation. Performance
claims require fresh equal-budget measurements for the current implementation.
