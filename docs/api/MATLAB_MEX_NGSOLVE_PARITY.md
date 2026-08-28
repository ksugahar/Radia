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
| Public top-level pybind11 names | 99 |
| Covered mapped MEX names | 99 |
| Underscore numerical kernels | 28 / 28 covered |
| Stateful pybind11 class surface | 126 / 126 covered |
| Radia MEX gateway commands | 364 |
| Optuna MEX gateway commands | 20 |
| MATLAB Optuna classes | 38 |
| MATLAB Optuna functions | 30 |

The gateways are intentionally independent. `radia_mex` contains no
`optuna.*` commands and owns the Radia/NGSolve/oneMKL dependency graph.
`optuna_mex` contains two API commands and 18 optimizer commands, has no
NGSolve, oneMKL, Radia-core, or Python dependency, and is loaded only when
optimization is used. A missing or non-20-command `optuna_mex` fails loudly; it
is never replaced by a call through `radia_mex` or a silent MATLAB substitute.

Seven acoustic commands share their implementation with the pybind11
scattering/CQ bindings: four real-wavenumber sphere models, complex-wavenumber
soft-sphere scattering, BDF delta, and CQ grid construction. Full acoustic
CQ-BEM and FSI remain NGSolve/Python workflows and use the declared MATLAB
Python fallback rather than a second numerical implementation.

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

The MATLAB implementation is a differentially verified Optuna 4.9.0 subset,
not a drop-in replacement for the full Python package. It provides `ask`,
`tell` by Trial or trial number, FIFO `enqueueTrial`,
define-by-run suggestions, intermediate reporting, table-backed persistence,
and Simulink execution. `Study.optimize` supports timeout, callbacks,
selective exception catching, `TrialPruned`, and callback/objective `stop`;
`FrozenTrial`, `createTrial`/`create_trial`, `addTrial`, study/trial attributes,
and metric names cover the listed Optuna lifecycle surfaces. `bestValue`,
`bestParams`, and `bestSolution` expose the
persisted single-objective best for warm starts and downstream Simulink setup;
multi-objective studies use `paretoFront` or FrozenTrial-returning
`best_trials`. The MAT table schema and Simulink integration are MATLAB-only
contracts and are not evidence of upstream compatibility;
explicitly seeded Random, univariate and multi-objective TPE, Grid, scrambled
QMC, BruteForce, CMA-ES, NSGA-II, and NSGA-III proposal sequences are locked
to generated Optuna 4.9.0 differential fixtures. GP acquisition after its
startup phase and persisted-history replay use the pinned upstream-Python
backend by default and are covered by that proposal guarantee; the optional
MATLAB-native GP backend remains integration-only evidence.

The separate `optuna49_mcp_oracle.json` fixture is generated by launching the
official `optuna/optuna-mcp` 0.2.0 entry point and calling its tools through an
MCP stdio client. It checks the shared Study/Trial surface (`create_study`,
fixed-distribution `ask`, `tell`, attributes, `add_trial(s)`, metric names,
`best_trial`, and `best_trials`) without importing the MCP server in-process.
Because the official MCP `set_sampler` tool exposes only the sampler name and
not a seed, seeded proposal parity remains owned by the direct Optuna 4.9.0
fixture rather than being approximated through an unseeded MCP run.

Unseeded MATLAB sampler constructors match the upstream `seed=None` boundary
by taking fresh private entropy per instance without modifying MATLAB's global
RNG. Since both sides are intentionally nondeterministic in that mode, exact
proposal-sequence comparisons still use the same explicit seed.

The complete-package denominator is generated from pinned Optuna 4.9.0 in
`tests/matlab/fixtures/optuna49_public_api.json`; the MATLAB comparison lives in
`matlab/optuna49_api_coverage.json`. The radia-mcp contract checks the inventory
digest and reports missing, partial, and unmapped entries. It cannot report
complete compatibility until all three categories reach zero.

The native sampler set includes Random, TPE/MOTPE, CMA-ES, constrained joint
NSGA-II and NSGA-III, `GPSampler`, finite `GridSampler`, define-by-run
`BruteForceSampler`, `PartialFixedSampler`, and numeric Sobol/Halton
`QMCSampler`. Unscrambled QMC uses the same leading Sobol/Halton points as the
Optuna/SciPy definitions. `Scramble=true` delegates to the configured
`scipy.stats.qmc` Sobol/Halton implementation so its PCG64 scrambling and
fixed-seed sequence are the same as Optuna 4.9.0; this path fails loudly if
SciPy is unavailable. The MATLAB-native unscrambled Joe-Kuo Sobol direction
table is currently bounded to 32 dimensions. The pruner set includes Median, Percentile, Threshold, Patient,
Nop, asynchronous Successive Halving, Hyperband, and paired Wilcoxon. Rung,
QMC sequence, GP hyperparameter, and evolutionary sampler state is persisted.

`GPSampler` fits toolbox-free Matérn-5/2 ARD Gaussian processes over the
stable intersection search space. It supports numeric and categorical
variables, expected improvement, Monte-Carlo expected hypervolume
improvement, `c_i <= 0` feasibility probability, pending-trial repulsion,
deterministic/noisy objectives, Cholesky jitter recovery, and seeded resume.
Kernel hyperparameters use Optuna 4.9's default inverse-length-scale, kernel,
and noise priors. Acquisition search uses 2048 scrambled Sobol preliminary
points, ten bounded local starts, and 128 Monte-Carlo EHVI samples by default.
The local optimizer and random stream remain MATLAB-native counterparts of
Optuna's L-BFGS-B/Torch implementation, not numerical or random-stream bit
parity. Differential gates use identical Sobol startup trials for Branin,
Hartmann6, constrained search, multi-objective hypervolume, and pending batches;
they also compare deterministic kernel-fit and NSGA-III transform fixtures to
Optuna 4.9. `NSGAIIISampler` retains the constrained joint
NSGA-II generation/crossover contract and replaces cutoff-front crowding with
Optuna-style extreme-point normalization, reference-line association, and
niche preservation. `BruteForceSampler` exhausts finite conditional
define-by-run trees and prevents duplicate RUNNING branches by default.
The Simulink auto-policy uses GP for declared fixed-numeric budgets up to 250
trials, CMA-ES beyond that for unconstrained single-objective searches,
NSGA-II/III for larger multi-objective populations, and TPE/MOTPE for
categorical, conditional, or otherwise unknown spaces. The selected rule is
stored with the study as `auto_sampler_decision`.

The optimizer density/ranking kernels are native and live exclusively in the
independent `optuna_mex` gateway:
`optuna.pareto.rank_crowding`, `optuna.parzen.log_pdf_numerical`, and
`optuna.parzen.log_pdf_categorical`. The NumPy-compatible MT19937 stream is
also a checked native handle with create, uniform/normal/integer/permutation,
snapshot/restore, and destroy commands. They are discovered through
`radia.optuna.nativeStatus` and used automatically. Four fused TPE commands
(`best_numerical`, `best_joint`, and their observation-building variants) keep
candidate generation, Parzen-density evaluation, and acquisition selection on
the existing checked random-state handle. Grouped TPE additionally keeps an
incremental completed-trial history on that handle and generates every current
co-occurrence group through one `best_grouped_history` preparation call. The
fast path is deliberately limited to sequential in-memory studies with default
gamma/weights and no constraints, constant liar, or categorical-distance
callback; every other contract uses the high-level MATLAB orchestration path.
That selection occurs above the required `optuna_mex` boundary and is not a
missing-binary or legacy-gateway fallback. The high-level study, trial, CAE
artifact contract, and public MAT-table views remain in MATLAB. Internal study
state uses typed column arrays and materializes compatibility tables lazily.
This is deliberate kernel acceleration rather than an opaque second optimizer
implementation.

The 2026-08-25 LAB benchmark used MATLAB R2026a Update 3, Python 3.12.10,
Optuna 4.9.0, 100 sequential trials per workload, 11 repeats, and the median
of the final eight repeats after three warm-ups. Explicit seeds and final
checksums were identical on both sides. Table persistence and parallel
scheduling were disabled so this measures the shared sequential Optuna path.

| Workload | Python Optuna 4.9.0 | MATLAB Optuna | Python/MATLAB speed ratio |
|---|---:|---:|---:|
| Scalar TPE | 0.104030 s | 0.093358 s | 1.114x |
| Grouped conditional TPE | 0.327160 s | 0.198069 s | 1.652x |

MATLAB therefore exceeds upstream warmed throughput in both representative
workloads while preserving the upstream proposal checksum. The speedup comes
from removing `table` mutation and `datetime` construction from the hot path,
making public tables lazy views, and fusing repeated TPE estimator/candidate
work in MEX. Incremental grouped history removed the per-trial full-history
projection. Seven fresh MATLAB processes measured
the first `optuna_mex` call at a 0.007793-second median; individual measurements
ranged from 0.006980 to 0.037714 seconds as filesystem-cache and host activity
varied. The MEX is
151552 bytes and imports only MATLAB, VC runtime, Windows CRT, and Kernel32
libraries. MATLAB process launch remains outside both the warmed-throughput and
MEX-load measurements. The reproducible scripts and machine-readable result
live under `validation_test/optimization`.

Sampler constraint callbacks use Optuna's `c_i <= 0` feasibility convention.
The resulting vectors are persisted in `Study.ConstraintTable`. Multi-objective
TPE prioritizes feasible trials, ranks infeasible trials by summed positive
violation, and uses hypervolume contribution for boundary selection and good
density weights.

`NSGAIISampler` matches the Optuna 4.9 algorithm contract rather than the
former parameter-wise approximation. It uses joint parents, COMPLETE-only
generations, fixed parent caches, constrained dominance, cutoff-front crowding,
categorical uniform crossover, omission-based random mutation, and the six
built-in Optuna crossover families. `ConstraintCountTable` preserves the
semantic difference between an explicitly empty (feasible) constraint vector
and a missing vector. Sampler state schema v3 persists the crossover,
probabilities, strategy identities, generation assignments, parent caches, and
random state. The NumPy-compatible MT19937 stream, parent ordering, mutation
consumption order, and NSGA-III niche shuffles are locked to Optuna 4.9.0
fixed-seed proposal fixtures.

`radia.optuna.SheetMetalRunner` connects this lifecycle to the native
`optimizeHexSheetTopology`, `optimizeHCurlEddyBubbleHexSheet`, and
`optimizeHCurlEddyBubbleActivationHexSheet` drivers. The matching Simulink
block advances one complete trial per sample and sends incumbent-update and
Pareto telemetry to standard Scope/XY Graph blocks. Solver execution remains
MATLAB/MEX plus NGSolve/Cubit; there is no Python call in that optimization
loop.

The MATLAB TPE follows the Optuna 4.9 univariate TPE defaults and Parzen
construction, including startup trials, gamma/max-good selection,
truncated-normal kernels, history weights, an explicit prior, magic clipping,
and EI candidate sampling. `TPESampler` keeps `Multivariate=false` by default.
With `Multivariate=true`, it automatically infers the intersection search
space from persisted COMPLETE and PRUNED trial distributions and jointly
samples ordinary `suggestFloat`, `suggest_int`, and `suggest_categorical`
calls. Conditional parameters or changed distributions remain outside the
intersection and use their ordinary independent proposal. `Trial.suggestVector`
remains an explicit numeric convenience rather than a requirement for joint
sampling. The implementation is algorithmically aligned with the reference
Python optimizer. Its single-objective mixed-distribution and
multi-objective fixed-seed proposal sequences are checked directly against
Optuna 4.9.0. The MATLAB CMA-ES remains an independent MATLAB implementation,
but its RandomState draws, per-trial reseeding, pre-checkpoint seed
consumption, and standard full-covariance proposal sequence are likewise
oracle-locked. It follows the standard
full-covariance generation update: weighted recombination, cumulative
step-size and covariance evolution paths, rank-one/rank-mu covariance
adaptation, bounded resampling, and positive-definite covariance repair. It
automatically uses the numeric intersection search space and persists complete
and partial-generation state in `Study.SamplerStateTable`. Python-only CMA
variants such as separable CMA, margin handling, and restart strategies remain
separate future extensions rather than silent aliases. Use
`validation_test/optimization/validate_matlab_optuna_quality.m` to compare
fixed-seed, equal-budget regret, Pareto-front error/coverage, and ask/tell cost
before making sample-efficiency claims. For expensive CAE the field solve is
expected to dominate this bookkeeping, but that fact does not establish equal
search quality.
