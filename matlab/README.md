# Radia MATLAB and Simulink

Radia's final human-facing application interface is the single **Radia**
Simulink library. Its Electromagnet, PCB PEEC, Motor, Stream Function, and
Induction Heating blocks share one Library Browser entry. Electromagnet, PCB
PEEC, Motor, and Stream Function launch the validated Python/headless CLI once
on an explicit rising trigger. Induction Heating instead uses separate native
Eddy and Thermal MEX S-Functions and has no Python fallback in its simulation
loop.

The additional `Applications/Stream Function Optimization` subsystem keeps
that explicit batch boundary and adds MATLAB-native MMA/SQP over preassembled
stream-function response, regularization, and current-constraint operators.
It exposes the complete objective/constraint history as fixed-size Simulink
vectors for XY Graph display. Python is not launched from an optimization
sample step; the headless design/Pareto lane and analytic-adjoint lane have
separate explicit triggers.

The tracked standalone composition is `radia_streamfunction_optimization.slx`.
Regenerate it with `radia.simulink.buildStreamFunctionOptimizationModel`; the
`.m` generator, rather than manual SLX editing, owns its block layout and mask
contract.

When an application run computes a spatial field, the runner writes a checked
GMSH `.msh v4.1` post-processing artifact inside the run directory and lists it
in `result.json`. The GMSH path is runner-owned rather than a mask parameter.
Scalar/circuit-only modes report the artifact as not applicable instead of
inventing an empty field file.

Native MEX remains an optional promotion for the four batch application
blocks. It is mandatory for the IH runtime. Every promoted path must pass error
propagation, object lifecycle, repeated-run, and long-run stability tests.

## Optional MEX Interface

`radia_mex` exposes Radia C++ kernels directly to MATLAB without routing
solver calls through a Python process. NGSolve still owns TaskManager, mesh loading,
finite-element orientation, curved mappings, and HCurl/HDiv space updates.

On the current Windows distribution, the MEX gateway does not link
`python312.lib` and MATLAB does not launch `python.exe`, but the prebuilt
pip NGSolve `libngsolve.dll` transitively requires `python312.dll`. Full
Python-DLL independence requires rebuilding NGSolve/Netgen without Python
support; it is not claimed by this MEX target.

## Build the MEX Target

Use the checked MSVC environment wrapper for the focused MEX build:

```powershell
pwsh -ExecutionPolicy Bypass -File .\Build.ps1 -MatlabMexOnly -Verbose
```

The equivalent manual CMake route, from an initialized Visual Studio x64
developer environment, is:

```powershell
cmake -S . -B build-msvc `
  -DRADIA_BUILD_MATLAB_MEX=ON `
  -DMatlab_ROOT_DIR="C:\Program Files\MATLAB\R2026a"
cmake --build build-msvc --config Release --target radia_mex -j 8
```

The build copies `radia_mex.mexw64` into this directory. On Windows it also
places the required oneMKL dispatch DLLs beside the MEX; these local binary
artifacts are intentionally excluded from Git. MATLAB's existing OpenMP
runtime is reused, so `libiomp5md.dll` is not copied beside the MEX.

In MATLAB:

```matlab
addpath("<radia-repository>\matlab")
radia.setup()
info = radia.ngsolve.space_info("model.vol", 6);
[M, matrixInfo] = radia.ngsolve.matrix_dump( ...
    "model.vol", "hcurl", 6, "mass", no_grads=true);
```

`radia.ngsolve.matrix_dump` assembles a native NGSolve `T_BilinearForm<double>` and
returns a MATLAB sparse matrix with NGSolve's global DoF ordering unchanged.
Persistent spaces support `h1`, `vectorh1`, `hcurl`, and `hdiv`.
`vectorh1` is the native mesh-deformation space; matrix assembly supports
`h1`, `hcurl`, and `hdiv` with `mass` and `stiffness`, plus `curlcurl` and
`divdiv` for the corresponding vector spaces. The optional second output records the mesh dimension, order, DoF
count, structural nonzero count, numeric nonzero count, and TaskManager
thread count. The MEX gateway uses 1-based row/column triplets so no MATLAB
index conversion is needed by callers.
The previous `radia.spaceInfo` and `radia.ngsolveMatrix` spellings remain
available as compatibility aliases.

### Acoustic analytic kernels

The sphere-scattering references and convolution-quadrature primitives share
one C++ implementation between pybind11 and MEX. MATLAB wrappers live under
`+radia/+acoustic`:

```matlab
points = [0, 0, 1.4; 1.3, 0, -0.2];
scattering = radia.acoustic.softSphereScattering(2.3, 0.9, points, Terms=28);
grid = radia.acoustic.cqGrid(256, 1e-4, Method="BDF2");
```

The native commands cover soft, rigid, fluid, and elastic spheres, complex
wavenumbers, BDF delta, and CQ grid construction. Full CQ-BEM and acoustic FSI
still require NGSolve's Python API and therefore use the explicit in-process
fallback:

```matlab
result = radia.python.acoustic("cq", "cq_soft_sphere_scattering", ...
    {points}, Keywords=struct("radius", 0.9, "order", 2));
```

The fallback reports `backend="python-fallback"` and the loaded Python
executable and library. It requires Python 3.12 in `pyenv` InProcess mode and
is allowed for initialization, explicit update, artifact generation, and batch
solves only. It must not execute once per Simulink time step.

### Axisymmetric Henrotte element matrices

The native `axifem` slice exposes Q1 and Q2 axis-aligned magnetic element
matrices without reproducing NGSolve's finite-element object model:

```matlab
mu0 = 4*pi*1e-7;
element = radia.axifem.q1MagneticElementMatrices( ...
    1e-3, 2e-3, -0.5e-3, 0.5e-3, mu0, 5.8e7);
K = element.stiffness;
M = element.sigma_mass;

q2 = radia.axifem.q2MagneticElementMatrices( ...
    0, 1e-3, -0.5e-3, 0.5e-3, mu0, 5.8e7);
```

The Q1 node order is `(ra,za),(rb,za),(rb,zb),(ra,zb)`. Q2 appends the four
edge-midpoint nodes and the face-centre node, using the midpoint convention in
`s=r^2`. All DOFs are nodal `A_phi` values. The MEX command, pybind11 array
API, and NGSolve Q1/Q2 BFIs use one C++ source. Mesh, FESpace,
CoefficientFunction, BilinearForm, and global assembly remain NGSolve-owned.

For a compact transient teaching model, the Q2 matrices can be converted once
to an exact-ZOH state-space system and run by the native MEX S-Function:

```matlab
state = radia.simulink.makeAxiEddyElementModel( ...
    0, 1e-3, -0.5e-3, 0.5e-3, mu0, 5.8e7, SampleTime_s=1e-5);
radia.simulink.buildAxiEddyElementModel("axi_q2_eddy", state, Save=false);
```

Linear planar or axisymmetric coil matrices can use the same native state-space
MEX boundary. `makeCircuitFieldStateSpace` reduces the constrained field matrix
and signed branch source vectors to an RL model, while
`buildCircuitFieldStateSpaceModel` wires it to
`radia_state_space_mex_sfunction`. The field factorization is offline and no
Python call occurs per Simulink step.

Finite-element assembly and Python are not called per Simulink time step. A
global `.vol` assembly adapter is a separate promotion step.

### NGSolve CoefficientFunction and GridFunction handles

MATLAB can keep NGSolve field objects alive as checked native uint64 handles.
The finite-element space, element mappings, and expression tree remain native;
MATLAB transfers only the values explicitly requested:

` ` ` matlab
cf = radia.ngsolve.CoefficientFunction.constant(2.5);
cf2 = 4.0 * cf + cf;

gfu = radia.ngsolve.GridFunction.create( ...
    "model.vol", "h1", 2, NoGrads=true, Name="temperature");
gfu.setVector(zeros(gfu.DofCount, 1));
gfu.interpolate(cf);
coeff = gfu.asCoefficient();
values = gfu.vector();
` ` `

The current native slice supports real and complex constant coefficient
functions, NGSolve arithmetic expression trees (addition, subtraction,
multiplication, and scalar scaling), coefficient metadata, real/complex
GridFunction DoF vector access, coefficient interpolation, and
GridFunction-to-CoefficientFunction views. Coefficient functions can also be
evaluated at physical points through NGSolve's element search and mapped-point
evaluation. The delete method releases the native handle; a coefficient view
retains its GridFunction until the view is deleted. This is a native MEX path
and does not start a Python process. The current pip-provided NGSolve DLL
still requires the Python runtime DLL at load time.

### Persistent mesh, space, form, and matrix handles

For workflows that assemble once and iterate many times, MATLAB can retain the
NGSolve object graph itself:

```matlab
mesh = radia.ngsolve.Mesh.create("model.vol");
space = radia.ngsolve.FESpace.create(mesh, "hcurl", 6);
form = radia.ngsolve.BilinearForm.create(space, "mass");
matrix = form.matrix();

triplets = matrix.values();       % 1-based sparse triplets
A = matrix.sparse();              % optional MATLAB snapshot
x = matrix.vector();
x.setValues(ones(matrix.Cols, 1));
y = matrix.matvec(x);             % native BaseVector result
inverse = matrix.inverse();       % free-DoF inverse handle
```

`Mesh`, `FESpace`, `BilinearForm`, and `Matrix` use checked `uint64` handles;
dependent native objects retain the storage they need. `Matrix.matvec` and
`Matrix.inverse` do not rebuild the mesh or copy a MATLAB DoF array per native
operation. `values()` and `sparse()` are explicit observation boundaries for
inspection, export, and MATLAB-side diagnostics. The current persistent form
contract assembles one built-in volume form at a time (`mass`, `stiffness`,
`curlcurl`, or `divdiv`); scalar coefficient weighting and complex FESpaces
are supported. Tensor-valued/boundary forms and arbitrary Python callbacks
remain explicit unsupported boundaries.

An existing H1 space can also own a native GridFunction and a constant-source
RHS without reloading the mesh:

```matlab
temperature = radia.ngsolve.GridFunction.fromFESpace( ...
    space, Name="temperature");
rhsForm = radia.ngsolve.LinearForm.create( ...
    space, "constant", Value=2.0);
rhs = rhsForm.vector();
```

`LinearForm.createFromCoefficient` accepts a real or complex native
`CoefficientFunction` and assembles a volume right-hand side in real or complex
H1, HCurl, or HDiv spaces. NGSolve retains the test-function expression and
all element mappings. Boundary loads and general solver/preconditioner objects
remain explicit extension points.

For iterative solvers, keep the DoF vector native instead of fetching a
MATLAB array on every iteration:

```matlab
x = gfu.vectorHandle();
r = x.copy();
r.setZero();
r.axpy(2.0, x);
residual_norm = r.norm();
energy = r.dot(x, Conjugate=true);
```

`vectorHandle()` is a live view of one GridFunction component. `copy()` makes
an independent NGSolve `BaseVector`; `setZero`, `scale`, `axpy`, `dot`, and
`norm` operate in C++. Use `values()` or `setValues()` only at an explicit
MATLAB observation or control boundary. This keeps the iteration loop free of
per-step DoF-array copies while preserving a clear MATLAB escape hatch.

For the first native HCurl Eddy Bubble/CLN bridge, MATLAB can call the C++
NGSolve path directly. No Python process is involved:

```matlab
[~, info] = radia.ngsolve.matrix_dump("model.vol", "hcurl", 6, "mass");
ports = zeros(info.dof_count, 2);
ports(1, 1) = 1;
ports(2, 2) = 1;
basis = radia.ngsolve.hcurl_eddy_cln_native_basis( ...
    "model.vol", 6, ports, 4);
```

The native command loads the mesh, constructs the curved/high-order HCurl
space, assembles NGSolve's `mass + curlcurl` unit-shifted response operator,
factorizes it, and performs mass-metric orthogonalization in C++. MATLAB
receives only numeric arrays and metadata (`vectors`, `free_dofs`, `rank`, and
`orthogonality_error`). The result also carries `mass_gram`, `curlcurl_gram`,
and `port_rhs`, with the explicit projections `V'*M*V`, `V'*K*V`, and
`V'*ports`.

For a Python-free local HCurl diffusion CLN model, use the high-level builder:

```matlab
model = radia.ngsolve.hcurl_eddy_cln_model( ...
    "model.vol", 6, ports, 4, ...
    Conductivity=5.8e7, Reluctivity=1/(4*pi*1e-7));
```

This forms `R = Reluctivity*K_r` and `L = Conductivity*M_r`, then creates the
same passive discrete state-space contract used by Simulink. It is a local FE
diffusion projection. The VIM external Laplace/BEM inductance, topology-aware
surface/bridge eddy-bubble classification, and frequency-dependent SIBC DtN
rationalization remain separate explicit stages; they are not inferred from
the local mass and curl-curl matrices.

HCurl Eddy-Bubble topology optimization also has a direct MEX path. Construct
`radia.topopt.HCurlTopologyOperator` from the HACApK scalar Gram, reduced charge
maps, affine sub-tetrahedra, and their parent-cell map. Then use
`assembleHCurlResistanceShapeTangents` for NGSolve-owned Piola `R,dR`, or cache
`assembleHCurlCellCurlGrams` for activation optimization. The
`linearizeHCurlMultifrequencyJoule` and
`linearizeHCurlActivationMultifrequencyJoule` entry points solve all weighted
complex load cases and contract `dL` inside C++ without materializing a dense
directional inductance matrix. Their LP adapters feed the existing
`optimizeHexSheetTopology` GetTrafo/Cubit two-level loop; no Python call occurs
inside the MATLAB optimization iteration.

The regression suite is run with:

```powershell
matlab -batch "addpath('matlab'); results = runtests('tests/matlab'); assert(all([results.Passed]))"
```

`radia.setup` discovers the NGSolve runtime associated with the selected
Python installation and places its DLL directories before older Netgen DLLs.
This setup operation is cached; numerical calls enter C++ directly.

For continuous electromagnetic shape or material variables, use the checked
adjoint optimizer instead of treating every update as a black-box trial. The
evaluator returns one scalar objective, its design gradient, optional
inequalities in `c <= 0` form, and a design-by-constraint Jacobian:

```matlab
evaluate = @(rho) struct( ...
    "objective", jouleLoss(rho), ...
    "gradient", jouleAdjointGradient(rho), ...
    "constraints", volumes' * rho / volumeMax - 1, ...
    "constraint_jacobian", volumes / volumeMax);

qa = radia.topopt.checkAdjointGradient(rho0, evaluate);
result = radia.topopt.optimizeAdjoint(rho0, evaluate, ...
    Solver="mma", LowerBounds=zeros(size(rho0)), ...
    UpperBounds=ones(size(rho0)));
```

`checkAdjointGradient` is a QA-only directional finite-difference diagnostic;
the optimizer itself always consumes the supplied analytic derivatives. MMA
uses moving asymptotes and analytic convex-subproblem Jacobians, requires
finite bounds, and accepts inequalities. `Solver="sqp"` uses `fmincon` with
`SpecifyObjectiveGradient` and `SpecifyConstraintGradient` enabled and also
accepts equality constraints. Neither route silently estimates a missing
gradient.

`optimizeHCurlActivationAdjoint` is the first direct electromagnetic adapter.
It evaluates the native Eddy-Bubble multifrequency Joule MEX adjoint, adds a
dimensionless volume inequality, and sends the resulting material activation
problem to MMA or SQP without a Python call. A practical hybrid is TPE/CMA-ES
for discrete or global outer variables followed by MMA/SQP for continuous
field-sensitive refinement.

## Current production slice

- NGSolve TaskManager execution and HCurl/HDiv p-space construction
- native NGSolve H1/HCurl/HDiv mass, stiffness, curl-curl, and div-div matrix
  assembly with MATLAB sparse export
- native VectorH1 deformation handles with SetDeformation/UnsetDeformation and
  per-element GetTrafo Jacobian-ratio / spectral-condition sampling
- native NGSolve CoefficientFunction expression handles with constant,
  arithmetic, scaling, metadata, and GridFunction views
- native NGSolve GridFunction handles with real/complex DoF vectors,
  interpolation, metadata, and explicit lifetime management
- native NGSolve BaseVector handles with GridFunction views, independent
  work-vector copies, low-copy axpy/scale/dot/norm operations, and explicit
  value export/import
- native projected matrices, reduced block sums, diagonal preconditioners, and
  term-count diagnostics without a MATLAB sparse-array round trip
- `radia.RadiaField` B/H/A/M/Phi CoefficientFunctions with coordinate
  transforms, reusable point caches, and native trilinear voxel interpolation
- B-input `radia.EnergyStopMaterial` with checked native-handle lifetime
- complex dense solve and mixed-Galerkin Schur complement
- skin impedance and SIBC tail/termination kernels
- CLN Lanczos reduction, tridiagonal construction, impedance sweeps,
  loop-star coupling transforms, and ACA/SVD star compression
- PEEC filament geometry with a stateful HACApK inductance manager
- monopole HDiv charge-Gram H-matrix with regular, transpose, and symmetric matvecs
- persistent numeric HDiv/planar field evaluators for tet, curved-tet, cloud, and planar sources
- persistent ChargeGram charge-map/mass-map configuration, including direct
  `radia.ngsolve.Matrix` registration without MATLAB triplets, demag/mass
  matvecs, Riesz recovery, and diagnostics
- complete analytic TET/HEX/WEDGE ChargeGram directional derivatives, including
  self-volume, self-face, and TET charge-map row-rate terms
- persistent HACApK ChargeGram directional-derivative operators for topology
  optimization without materializing dense shape-derivative matrices
- EVRS/T-method algebra and tetrahedral HCurl reduced Gram integration
- complex filament/surface Biot-Savart kernels and P1/P2 SLDL Galerkin assembly
- stateful HACApK scalar BEM construction, build, matvec, and statistics
- shared planar charge-cloud field, vector-potential, Maxwell torque, and force kernels
- closed-form cuboid average-B and average-demagnetizing-tensor kernels
- static and harmonic equivalence-surface reconstruction with separate real/imaginary outputs
- HACApK H-LU truncation/parallel/accumulator controls, timing/materialization
  statistics, mixed-operation breakdown, cluster strategy, and self-tests
- legacy/pybind object creation: tetrahedron, hexahedron, wedge, pyramid,
  extruded polygon, cylinder magnet, rectangular/arc/racetrack/filament current,
  arc polygon magnet, background field, and containers
- legacy/pybind field utilities: field lists and integrals, center-field,
  Maxwell-tensor force, force surfaces, precision criteria, and length tolerance
- legacy/pybind nonlinear material and solver entry points: formula/table
  saturation, anisotropy, lamination, `MatMvsH`, `SolveNonl`, and solve statistics
- B-input hysteresis materials: Play/Energy construction from MATLAB cell tables,
  state save/restore/commit, batched forward/commit stepping, reversible
  reluctivity, and irreversible field output
- solver configuration parity through scalar MATLAB structs, including relaxation,
  Newton damping, B-input Newton/Hantila, and magnetization retention
- HACApK cluster strategy inspection through `radia.getClusterStrategy` and
  `radia.hluClusterStrategy`
- legacy/pybind overlap: linear/saturating material, solve, field, transforms,
  object state, and utility calls

`MatSatAniso` follows the native Radia coefficient-vector contract. Valid
parallel/perpendicular coefficient lengths are selected by the C++ material
implementation; an invalid combination is reported consistently by the MEX
boundary. `SolveNonl` is exposed for parity, while mesh-less soft-iron objects
remain rejected by the current solver and must use the NGSolve/Netgen FEEC path.

`ObjMltExtPgn`, `ObjMltExtRtg`, and `ObjMltExtTri` are not public MATLAB or
Python APIs. The first two use process-unsafe legacy geometry paths, and the
third depended on the removed Triangle bundle. Construct curved or extruded
geometry through Netgen or Cubit and pass the resulting mesh to the
NGSolve-native path.

`radia.BuildMatrix`, `radia.GetInteractMatrix`, and `radia.GetFaceGeom` expose
the legacy interaction-matrix inspection contract. Matrix construction may be
rejected for mesh-less objects; the NGSolve/Netgen finite-element matrix path is
the production route for HCurl/HDiv and curved high-order elements.

The MEX boundary uses MATLAB column-major arrays externally and converts them
to Radia's row-major `[target][source]` convention internally.

## Simulink and induction heating

Open the native IH model with:

```matlab
addpath("<radia-release>\matlab")
install_radia_ih()
```

`radia_ih.slx` exposes current, workpiece angle, and ambient temperature as
ordinary Simulink sources. Eddy produces a spatial heat-density vector;
Thermal owns the discrete temperature state and feeds it directly back to Eddy.
The update order is fixed as Eddy at the accepted temperature, conservative
transport from the previous angle to the current angle, then one implicit
thermal step. A linear magnetic configuration reuses its Eddy solution under a
current-only change and rescales heat by current squared.

The preview accepts linear or temperature-linearized operator data only.
Nonlinear BH iteration is not implemented and `bh_mode="nonlinear"` fails
loudly rather than substituting the dense linear solve.

The first distributable native IH package is a preview runtime. Its model ships
with a one-DOF diagnostic configuration and accepts checked MAT/JSON files that
contain preassembled native operators. `radia.simulink.makeIHNativeConfig`
validates row-major matrices, strict `.vol` reports, rotation mode, and thermal
state dimensions. Cubit `.vol` to PEEC/BEM-A/BIM/FEM operator assembly is the
remaining production gate and is not claimed by this preview.

LUT and lumped-state-space IH helpers are not reachable from `openIH`, are not
included in the IH release package, and are not supported production
interfaces.

### Reduced state-space MEX blocks

Separate reduced-order blocks use `simulink.state_space.create`, `info`,
`step`, `reset`, and `destroy`. They are intended for fixed Motor, MagLev, and
HCurl Eddy Bubble/CLN models, not for the native IH field runtime above. Moving
height-family interpolation remains a MATLAB S-function because its parameter
field changes at runtime.

The validated TEAM 28 levitation curve is also available as a fixed-50 Hz
control model. `force_N` follows the benchmark sign convention (negative is
upward), while `upward_lift_N` is positive upward. The current scaling is
quadratic around the 20 A reference current:

```matlab
team28 = radia.simulink.makeTeam28CLNLUT();
[force_N, lift_N] = radia.simulink.evaluateTeam28CLNForce( ...
    team28, 0.004, 20.0);
radia.simulink.buildTeam28CLNModel("radia_team28_cln", team28, ...
    SampleTime_s=1e-3, StopTime_s=1.0);
```

This is the fast LUT/LTI control path for Simulink and RL experiments. At
50 Hz the 3 mm aluminium disk is volumetric rather than SIBC-applicable, so
the high-fidelity path uses the reduced HCurl-VIM contract below. The Team28
specific mesh/matrix export and height-dependent force operator are still
required before replacing the force LUT in a full mechanical plant.

### Reduced HCurl Eddy Bubble/CLN bridge

NGSolve remains responsible for the mesh, curved high-order HCurl basis,
orientation, and assembly. MATLAB receives the reduced numeric matrices
`R`, `L`, and `P`, so the same p=6 HCurl Eddy Bubble/CLN model can be used by
MATLAB and Simulink without duplicating finite-element plumbing:

```matlab
reduced = radia.simulink.makeHCurlEddyCLNModel(R, L, P, ...
    SampleTime_s=1e-5);
coefficients = radia.simulink.solveHCurlEddyCLNHarmonic( ...
    reduced, 50, 20);
radia.simulink.buildHCurlEddyCLNModel("radia_hcurl_eddy_cln", reduced, ...
    StopTime_s=0.02);
```

The reduced equation is `(R+sL)c=-sP*i`. Its Simulink input is
`u=-di/dt` and its output is the reduced port response `P'*c`. The harmonic
solve calls the existing `hybrid_vim.solve` MEX kernel, while the Simulink
block uses an exact zero-order-hold discretization of the same passive model.
SIBC termination is not silently approximated: it must first be rationalized
before it can be included in this state-space contract.

The NGSolve/Python side can write the same contract without SciPy:
`radia.vim.ExportHCurlEddyCLNJSON(model, "hcurl_exchange.json", ... )`.
Load that file with `radia.simulink.loadHCurlEddyCLNModel`. When the exporter
also provides the reduced force operator, evaluate the time-average force by
calling `radia.simulink.evaluateHCurlEddyCLNForce(model, coefficients, i)`.
The Team28 validation driver accepts `--export-model` to create this exchange
file from its first p=6 HCurl-VIM case. The resulting force helper is a
height-fixed reduced operator; the family path below supplies the validated
height-indexed extension for a moving mechanical plant.

The moving family is now supported directly:

```matlab
family = radia.simulink.loadHCurlEddyCLNFamily( ...
    "team28_hcurl_eddy_cln_family.json");
snapshot = radia.simulink.interpolateHCurlEddyCLNFamily( ...
    family, height_offset_m, BuildStateSpace=true);
coefficients = radia.simulink.solveHCurlEddyCLNHarmonic( ...
    snapshot, 50, coil_current_A);
force_N = radia.simulink.evaluateHCurlEddyCLNForce( ...
    snapshot, coefficients, coil_current_A);
radia.simulink.buildHCurlEddyCLNFamilyModel( ...
    "radia_team28_hcurl_family", family, StopTime_s=1.0);
```

The family block input layout is `[minus_dI_dt; height_offset_m;
coil_current_A]`; its output is `[port_response; force_x; force_y; force_z]`.
The Python/NGSolve driver is needed only to assemble and export the family;
MATLAB/Simulink execution thereafter uses no Python process. The exporter
requires a common reduced state basis, which is why the Team28 sweep keeps
the p=6 disk basis fixed and varies the incident coil fields and force
operator with height.

The validation driver creates the family on the NGSolve/Python side with:

```text
python validation_test/maglev/team28_hcurl_vim_force.py \
    --export-family C:\temp\team28_hcurl_eddy_cln_family.json
```

The full 25-height sweep is validation-class work and should be run on the
designated compute host; a short `--height-offsets-mm` list is suitable for a
developer smoke.

## Table-backed optimization

The MATLAB side also provides an Optuna-style, table-backed study. The main
trial table is accompanied by normalized parameter, intermediate-value, and
user-attribute tables, and can be persisted to a MAT file:

```matlab
study = radia.optuna.create_study( ...
    study_name="ih-control", direction="minimize", ...
    sampler=radia.optuna.TPESampler(Seed=42), ...
    storage="C:\temp\ih-control.mat");
results = study.optimize(@(trial) ihObjective(trial), 30);
best = study.best_trial();
```

The best completed single-objective trial is retained in the MAT-backed tables
and can be reused after reopening the study:

```matlab
bestValue = study.bestValue();
bestParams = study.bestParams();
bestSnapshot = study.bestSolution();
```

`bestSolution` contains `trial_number`, `value`, `params`, `state`, and
`end_time`. With no completed trial, `bestValue` is `NaN` and
`bestSolution.available` is false. Multi-objective studies use
`paretoFront()` instead. The snapshot is derived from the persisted tables, so
it cannot become stale when a study is reloaded.

`Trial.suggest_float`, `suggest_int`, and `suggest_categorical` are compatible
spellings for the MATLAB camelCase methods. `Trial.report` and
`Trial.should_prune` support intermediate-value pruning through
`radia.optuna.MedianPruner`. `RandomSampler`, `TPESampler`, and
`CmaEsSampler` share the same ask/tell lifecycle. TPE uses the persisted
MATLAB tables for good/bad density proposals. CMA-ES infers the numeric
intersection search space, samples a full population jointly, and applies
full-covariance rank-one/rank-mu adaptation with cumulative evolution paths;
integer variables are quantized after sampling. Its mean, covariance, step
size, evolution paths, partial population, and random state are retained in
`Study.SamplerStateTable` and survive MAT-file reloads. The workflow and tables are
compatible with Optuna, while sampler random streams and optimizer internals
are not promised to be bit-for-bit identical to Python Optuna.

Set `Multivariate=true` to let TPE infer the intersection search space from
the persisted COMPLETE and PRUNED trials. Ordinary `suggestFloat`,
`suggest_int`, and `suggest_categorical` calls then participate in one mixed
joint proposal without requiring `suggestVector`:

```matlab
sampler = radia.optuna.TPESampler(Seed=42, Multivariate=true);
study = radia.optuna.create_study( ...
    direction="minimize", sampler=sampler, ...
    storage="C:\temp\multivariate-study.mat");
study.optimize(@(trial) objectiveWithScalarSuggestions(trial), 40);
```

Parameters that are absent from any eligible trial, or whose persisted
distribution changed, stay outside the intersection and use their ordinary
independent proposal. This preserves define-by-run conditional spaces across
study reloads. `suggestVector` remains available as an explicit numeric
convenience, but it is not required for multivariate TPE.

For electromagnetic CAE, multivariate TPE and CMA-ES are the primary
production samplers. TPE owns mixed, conditional, and discrete search spaces;
CMA-ES owns continuous correlated geometry and control variables. Random
sampling remains the quality baseline, while NSGA-II and the separate MOTPE
entry point remain alternative multi-objective validation routes rather than
additional production defaults.

Simulink is evaluated through `radia.optuna.SimulinkRunner`. Its
`ConfigureFcn` receives a `Simulink.SimulationInput` and a Trial, so the same
table-backed study can tune Radia/VIM parameters, LUT/LTI controller values,
and model variables without placing the optimizer inside a time-step block:

```matlab
runner = radia.optuna.SimulinkRunner("radia_ih_model", ...
    ConfigureFcn=@(in,trial) in.setVariable("frequency_Hz", ...
        trial.Params.frequency_Hz), ...
    ScoreFcn=@(simOut,trial) ihLoss(simOut));
study.optimizeSimulink(runner, 30);
```

This keeps optimization history queryable as MATLAB tables while Simulink
remains the plant and dynamic-system evaluator. For multiple objectives, pass
`directions=["minimize","maximize"]`, inspect `study.paretoFront()`, and attach
`radia.optuna.LiveMonitor` through `ProgressFcn`. Use
`radia.optuna.MOTPESampler` for Pareto-ranked multi-objective TPE or
`radia.optuna.NSGAIISampler` for non-dominated sorting, crowding, crossover,
and mutation. Use
`SimulinkRunner.optimizeParallel` (`parsim`) or
`LTspiceRunner.optimizeParallel` (`parfeval`) for parallel engineering trials.

LTspice coupling does not require Simscape. `radia.ltspice.run` accepts both
netlists and `.asc` schematics; `SchematicEditor` edits component values, and
the conversion path stages recursive local include/model/hierarchical files.
`radia.ltspice.netlistToSchematic("circuit.cir")` is the thin MATLAB wrapper
over the canonical `radia-spice-lab` Python converter. By default it converts
back with hidden LTspice and rejects the result unless the node-rename-invariant
Python topology check proves connectivity equivalence.
Use `runNoise`, `analyzeFFT`, and `runIntervals` for noise spectra, transient
FFT, and explicit capacitor-voltage/inductor-current handoff between long
windows. Run `radia.simulink.buildLibrary` once after installation and
`sl_refresh_customizations` to expose the single **Radia** entry in Simulink
Library Browser. Its `Applications` subsystem contains Electromagnet, PCB PEEC,
Motor, Stream Function, and Induction Heating blocks. The first four
batch-analysis blocks delegate to the same `DesignSpec` and validated headless
Python CLI used by MCP; a boolean rising trigger runs once and leaves `run.log`
and `result.json`, with solver artifacts when execution reaches that stage.
Create the versioned settings file with
`radia.simulink.writeApplicationConfig`. The Python backend is intentional
for those four blocks while their MEX routes mature. IH is the native-only
exception described above. Notebook workbenches are retired for all
applications, including IH.

The remaining library blocks delegate to the same tested MATLAB APIs.
The `Material Dictionary` block is the normal material entry point for a
`.vol`-based model. Define `radia_materials` as a MATLAB `dictionary` whose
keys are reusable material names and whose scalar-cell values come from
`radia.simulink.makeMaterialSpec`. Define `radia_region_materials` as a second
dictionary from exact `.vol` region names to material names. Model
initialization runs `check-vol`, verifies complete region coverage and the
optional mesh SHA-256, then compiles the variable-size setup data to a
fixed-width numeric `RadiaMaterialBus`. No string lookup, dictionary lookup,
Python call, or finite-element factorization occurs per Simulink step.

```matlab
air = radia.simulink.makeMaterialSpec();
steel = radia.simulink.makeMaterialSpec( ...
    BH_B_T=[0; 1.0; 1.5], BH_H_A_per_m=[0; 250; 1200]);
copper = radia.simulink.makeMaterialSpec( ...
    Conductivity_S_per_m=5.8e7, Density_kg_per_m3=8960, ...
    SpecificHeat_J_per_kgK=385, ThermalConductivity_W_per_mK=401);
radia_materials = dictionary(["air","steel","copper"], ...
    {air,steel,copper});
radia_region_materials = dictionary( ...
    ["air_domain","stator","phase_a"], ...
    ["air","steel","copper"]);
```

This keeps the human workflow at the Simulink level: geometry and labels come
from `.vol`; materials, windings, circuits, excitation, motion, control, and
requested outputs are composed in the model. Radia-NGSolve-specific Python or
finite-element-space construction is not required for normal use.

Define windings with `radia.simulink.makeWindingSpec` and a second MATLAB
dictionary. `Regions` selects exact `.vol` coil-side regions and
`RegionPolarity` records the signed current direction of every side. Turns,
parallel paths, resistance, and positive/negative circuit terminals compile to
`RadiaWindingBus`. The same terminal identity is consumed by the native
field/circuit MEX route and the LTspice interval route.

`RadiaMachineCommandBus` carries voltage/current excitation, rotor angle and
speed, translational position/velocity, load torque, and load force.
`RadiaMachineResponseBus` returns current, flux linkage, back EMF,
electromagnetic torque/force, and losses. Consequently Simulink or Simscape
Multibody owns the mechanical equations and controller; the electromagnetic
solver does not hide motion in a private application loop.

The `LTspice Circuit` block accepts vector `InputNames` and `OutputTraces`.
Each sample advances one reset-time LTspice interval and hands saved node
voltages and inductor currents to the next interval. Block state is isolated
per instance and cleared by the Simulink `Terminate` callback, including model
restart and parallel-process use.
For magnetic power-electronic plants, use
`radia.simulink.buildHystereticLTspiceBlock`. It performs waveform relaxation
inside every sample interval: LTspice supplies coil current, the Radia
Play/Energy material advances flux, and differentiated flux linkage returns as
a PWL back-EMF source. Circuit and hysteresis states are committed only after
the coupled residual converges.
The native reduced field/circuit MEX path and LTspice are alternate circuit
backends behind the same winding terminal identity. This allows a controller,
power converter, rotating mechanical plant, and magnetic model to be composed
as ordinary Simulink subsystems while LTspice remains available for detailed
switching-circuit cross-checks.

## Field Study: completing the FEMM physics surface

Magnetics alone is not a FEMM replacement. The default Simulink interface also
exposes electrostatics, DC or harmonic current flow, steady heat conduction,
and time-harmonic magnetic eddy current through `Applications/Field Study`.
Geometry and region/boundary names come from one checked Netgen `.vol` file;
the existing material and winding dictionaries remain the only property and
coil identities.

```matlab
study = radia.simulink.makeFieldStudySpec( ...
    Physics="steady_heat",Formulation="axisymmetric",ElementFamily="P2", ...
    DirichletValues=struct("cooled",293.15), ...
    RobinBoundaries=struct("outer",struct( ...
        "transfer_w_per_m2_k",20,"ambient_k",293.15)), ...
    VolumetricSources=struct("workpiece",2.5e6));
contract = radia.simulink.compileFieldStudy(study,radia_material_contract);
requestFile = radia.simulink.writeFieldStudyRequest( ...
    contract,"C:\temp\radia_field_study\heat.json");
radia.simulink.writeApplicationConfig("field", ...
    struct("request_file",requestFile), ...
    "C:\temp\radia_field_study\run.json");
```

`harmonic_eddy` uses `(K+j*omega*M_sigma)a=S*i` and accepts complex RMS
winding-current phasors. It reports flux linkage, induced branch voltage,
magnetic energy, Joule loss, residual, and the independent real-power closure.
The current implementation deliberately rejects nonlinear or hysteretic
harmonic permeability instead of silently using an initial slope. Every field
study writes versioned JSON, four-part timing, Gmsh `.msh v4.1`, and `.geo/.opt`
launch companions. These lanes establish the default interface; final FEMM
retirement still requires frozen live artifacts for each physics/formulation
combination.
On Windows, the launcher uses `pwsh Start-Process -WindowStyle Hidden` rather
than the `CREATE_NO_WINDOW` process flag: LTspice 26 can fail with
`0xC0000409` under an SSH session when that flag is used. This hidden-window
route is verified on LAB and 100号機 and remains compatible with batch `-Run
-b` execution.
The executable details and limitations are reported by the
`matlab_optuna_simulink_contract` MCP tool.

CAD material topology optimization uses Cubit element IDs and volumes,
`radia.topopt.linearizeVIM` for analytic VIM response linearization, and
`radia.topopt.optimizeVIMLP` for volume/move-limit constrained LP updates.
`radia.topopt.writeCubitJournal` returns the final solid/void assignment to
Cubit. Generate the executable workflow with the `matlab_cad_topology_build`
MCP tool.

For sheet topology with geometry motion, use the two-level production driver:

```matlab
result = radia.topopt.optimizeHexSheetTopology( ...
    initialState, @linearizeStep, @makeVectorH1Deformation, ...
    @evaluateObjective, @rebuildHMatrix, cubitBackend, elementSizes, ...
    InnerIterations=10, MaxOuterIterations=10, ...
    ActivationRemoveThreshold=0.35, ActivationRestoreThreshold=0.65, ...
    CubitBatchInterval=5, CubitBatchFraction=0.05);
```

The inner batch is restricted to 5--20 iterations. Activation remains
continuous throughout that batch while native `Mesh.setDeformation` and
`Mesh.trafoQuality` sample NGSolve `GetTrafo`. The 0.35/0.65 activation
hysteresis prevents remove/restore chatter. Pending topology changes remain
continuous until either their fraction or age reaches the Cubit batch limit.
Cubit is called only when the batch must commit those changes, an explicit rebuild is
requested, a converged deformation needs a CAD commit, or the Jacobian gate
rejects in-place deformation. `rebuildHMatrix` is then called exactly once
after Cubit and never inside the inner loop.

## Reinforcement learning

The same MEX boundary can be used as a fast reinforcement-learning
environment. `radia.rl.Environment` is independent of Reinforcement Learning
Toolbox, while `radia.rl.toFunctionEnv` adapts it to `rlFunctionEnv` when that
toolbox is installed:

```matlab
environment = radia.rl.Environment( ...
    ResetFcn=@resetNativeState, ...
    StepFcn=@stepNativeSolver, ...
    MaxSteps=2000);
[observation, info] = environment.reset();
[nextObservation, reward, isDone, info] = environment.step(action);
```

The callbacks own an explicit native solver state and advance it exactly once
per environment step. IH production models use the distributed Eddy and
Thermal MEX S-Functions; the retired lumped thermal plant and IH LUT helpers
are not RL backends. `toFunctionEnv` connects the same checked callback
contract to MATLAB Reinforcement Learning Toolbox when it is installed.

## Native NGSolve field and matrix handles

Radia source fields can be passed directly into NGSolve assembly without
sampling through Python or MATLAB on every evaluation:

```matlab
B = radia.RadiaField(objectId, "b", Units="m");
B.prepareCache(observationPoints);
Bvoxel = B.asVoxelCoefficient(mesh, 61);
form = radia.ngsolve.BilinearForm.createFromCoefficient( ...
    space, "mass", Bvoxel);
```

The same persistent matrix handle supports low-copy model reduction and native
preconditioning:

```matlab
Mr = radia.ngsolve.Matrix.projected(form.matrix(), projection);
P = Mr.diagonalPreconditioner(1e-14);
x = Mr.vector();
x.setValues(reducedValues);
y = Mr.matvec(x);
```

For geometry optimization, a configured `radia.HACApKChargeGram` evaluates the
complete directional derivative in C++:

```matlab
dGtet = gram.tetChargeGramDirectionalDerivative(tetVelocity, faceVelocity);
dGhex = gram.hexChargeGramDirectionalDerivative(hexVelocity, quadVelocity);
dGwedge = gram.wedgeChargeGramDirectionalDerivative(wedgeVelocity, faceVelocity);
dGop = gram.directionalDerivativeOperator( ...
    "hex", hexVelocity, quadVelocity, AcaEps=1e-8);
directionalProduct = dGop.matvecSym(chargeVector);
products = gram.directionalDerivativeContractions( ...
    "hex", batchedHexVelocity, batchedQuadVelocity, left, right);
```

The batched contraction returns one `left.'*dG(k)*right` value per deformation
mode. It traverses the already-built parent H-matrix leaves and contracts ACA
factors immediately, so neither dense derivatives nor derivative H-matrices are
materialized in the optimization loop.

## Binding policy

The executable parity audit compares three pybind11 surfaces with the
`radia_mex` command table: 94 public top-level names, 27 underscore-prefixed
numerical kernels, and 111 stateful class members. All 232 entries are covered
by the current 311-command gateway. Three internal mesh/test helpers are
classified explicitly rather than silently omitted. The remaining `radentry`
C ABI is not a backward-compatibility contract: dead or unsafe entries are
deleted rather than retained. These families are represented as follows:

| Family | MATLAB representation |
|---|---|
| Legacy object, transform, material, solve, and field calls | Thin numeric wrappers over `radentry` |
| Pure array kernels such as SIBC, Schur, CLN, EVRS/T, Biot-Savart, and Galerkin assembly | Direct MEX commands |
| Stateful EnergyStop, HACApK, and field evaluators | Typed `uint64` handles with `mexLock` lifetime |
| NGSolve mesh and finite-element work | Persistent Mesh/FESpace/BilinearForm/Matrix handles plus numeric snapshots |
| NGSolve `CoefficientFunction` / `GridFunction` | Native MEX handles with explicit vector and metadata contracts |
| Python callbacks and interactive adapters | Keep in Python; expose numeric or handle-based equivalents |
| Underscore numerical kernels required by production methods | Direct MEX commands, covered by the parity audit |
| Mesh-plumbing and Python-only regression helpers | Explicit audited exclusions; use native NGSolve handles or tests |

Module-level parity is checked separately by
`tools/check_matlab_python_parity.py` against
`matlab/python_api_parity_manifest.json`. Every production `.py` module is
classified as `native-mex`, `matlab-native`, `python-fallback`, or
`private/not-applicable`; a complete low-level MEX inventory does not imply
that every high-level Python composition API has already become native MATLAB.
The current gaps and measured backend comparison are recorded in the
[MATLAB/Python API parity audit](../docs/api/MATLAB_PYTHON_API_PARITY.md).

This is deliberately not a mechanical one-MEX-file-per-`.def` translation.
One gateway avoids duplicating the Radia/NGSolve runtime and gives every
stateful class the same checked lifetime and exception boundary.

The parity audit also rejects any reappearance of the retired
`ObjMltExtPgn`, `ObjMltExtRtg`, or `ObjMltExtTri` names in pybind11, MEX, or the
C ABI declarations, definitions, and export table. This contract covers
Radia's pybind numerical surface and the NGSolve operations required by Radia;
it does not duplicate the complete upstream NGSolve Python API.
