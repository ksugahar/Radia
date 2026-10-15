# IH resident CUDA development validation

## Research release acceptance (2026-10-13)

The user explicitly permits an experimental temperature-dependent PEEC-BEM
(ESIM) release without independent full-physics certification. The historical
requirements below do not make that certification a research-release blocker.
Implementation consistency, checked convergence/true residuals, per-element
temperature and heat state, complete reaction coupling, Simulink repeated-run
and cleanup checks, exact-artifact distribution and documented limitations
remain required. This permission does not certify the physical model or waive
unsupported-input errors. Already-running NGSolve A100 evaluation is authorized;
this task never starts or stops that host.

`ih_temperature_material.cuh` stages an explicit piecewise-linear material
table in kelvin at each element's evaluated temperature. Temperature FE
coefficients pass through an explicit evaluation map; coefficients are not
assumed to be nodal temperatures. Resistivity, complex initial permeability,
saturation and common-knot cubic B(H) coefficients have declared interpolation.
Out-of-table temperatures are rejected; there is no implicit extrapolation.
`finiteElementESIMCells` reuses the existing finite-cell arithmetic with
element-specific parameters and checked workspace stride.

`run_temperature_material.py` checks three independently authored generic
material cells, distinct mapped temperatures, repeated calls and rejected
out-of-table input. It is a component check, not complete PEEC-BEM, nonlinear
temperature-dependent B(H), Simulink or release acceptance. The material
primitives are not yet connected to the persistent public MEX context.

`run_resident_smoke.py` builds a validation-only mexcuda gateway and checks the
resident complex CSR action, scalar-drive RHS and coherent quadratic loss
against independent NumPy complex algebra. Ten iterations retain operators and
state on device, with no intermediate vector transfers or explicit stage
synchronization. The loss includes basis cross terms. This is not a physical
PEEC-BEM(ESIM) solve, pybind11 physics parity, thermal time integration or
Simulink acceptance. Those remain mandatory before the radia-gpu release.

Small builds and correctness checks may run on 100 under low load with a
below-normal-priority owned MATLAB session and one compiler invocation, without
replacing installed runtimes. Do not compete with active heavy CI computation.
Formal transfer-inclusive time and memory comparisons require an otherwise idle
100 and are diagnostic provenance, not formal idle-historical_cpu_worker performance acceptance.
Results are always written as JSON. A100 is never started or connected to.

Next connect checked persistent context ownership, complete complex PEEC-BEM
solve/residual checks, thermal state updates and shared Simulink block state.
Keep electromagnetic/thermal vectors resident between updates; transfer only
required controls and observed outputs. Preserve block output/update ordering,
SI units, ESIM physical contracts and independently validated CPU results.

Device-scalar BiCGStab stages now compile and pass the identity-system first-
correction convergence check, including bypass of the zero omega denominator.
Full cuBLAS device reductions, the persistent solve loop, true-residual gate,
preconditioning, thermal evolution and Simulink context are still pending.
The stage check is not iterative-solver or physical IH acceptance.

The next compiled validation connects complete complex BiCGStab (cuBLAS device
reductions, eight-iteration CUDA graph batches, final true residual) to coherent
loss, device heat-load projection and a backward-Euler temperature increment
solve. The 3-unknown synthetic system agrees with NumPy direct algebra and
passes true-residual gates for zero, positive and negative excitation. This
checks numerical plumbing, not a PEEC-BEM physical model. The context is currently
retained only during the validation call; checked persistent MATLAB ownership,
Simulink output/update lifecycle and actual CPU/pybind11 physics comparison are
still required. No timing/performance claim follows from this result.

`persistent-state-20261012.json` checks ten coupled updates with one uploaded
operator set, device-only electromagnetic/loss/thermal stage exchange, accepted
state preservation on invalid controls, reset and rejection of a missing dynamic
material implementation. The device state is currently scoped to a single
validation MEX invocation. Checked persistent MATLAB handles, actual ESIM cell
and reaction-field updates, independent physical CPU parity and Simulink
acceptance remain pending. No full IH or release acceptance is claimed.

## Native nonlinear ESIM cells

`run_esim_cell_smoke.py` differentially checks the unmodified Python
`ESIMFiniteSlabSolver` against device finite differences for slab and cylinder
cells, constant real/complex permeability and monotone B(H) piecewise cubics.
Each CUDA thread solves one cell with the original boundary stencils, Picard
relaxation, convergence criterion and midpoint active/reactive power integration.
This is finite-difference cell arithmetic, not a replacement for NGSolve FE
semantics. Nonconvergence and linear true residual above 1e-6 fail loudly.

The validation supplies interpolation coefficients from the independent CPU
oracle. A MATLAB/native B(H) coefficient constructor, full material-state and
reaction-field coupling, checked persistent MATLAB handles, pybind11 model
parity and Simulink acceptance are still required. This validation is not a
Python prerequisite for the eventual GPU product, and its gateway is not shipped.
The official MATLAB MCP currently reports `failed to attach to MATLAB session`;
no production SLX acceptance has been performed. Engine-owned sessions were used
only for these standalone CUDA development checks.

## Current checkpoint: FP32 preconditioning and persistent Simulink state

This checkpoint supersedes the earlier pending implementation notes above.
The checked persistent mexcuda context now owns electromagnetic, loss and
thermal device state across MATLAB calls. Fixed right preconditioning can store
and apply its inverse-action matrix in FP32 while the original operator, outer
BiCGStab vectors/reductions and final true residual remain FP64. The supplied
inverse-action matrix is not yet a native IC/AMG builder.

`results/fp32-preconditioner-20261012.json` records nine algebra cases (no,
FP64 and FP32 preconditioning, three drives, ten thermal updates), independent
NumPy comparisons, true residuals, iteration counts and owned buffer sizes.
It does not contain timing acceptance or sampled peak memory.

Official MATLAB MCP attachment was restored in an operation-owned session.
The official Simulink Toolkit checked the disposable resident-state model;
save/close/reopen and repeated simulation passed with zero live contexts after
simulation. `results/fp32-simulink-algebra-20261012.json` records the FP32
preconditioned temperature and heat trajectories against independent NumPy
backward Euler algebra. This is synthetic orchestration validation, not physical
PEEC-BEM/ESIM parity. Actual native geometry assembly, nonlinear ESIM/reaction
coupling and independent pybind11 physical acceptance still block release.
A100 was not used. No speedup is claimed from these checks.

## Constitutive iteration checkpoint (2026-10-13)

`run_esim_coupling.py` exercises temperature material updates, ESIM cells,
complex solves and thermal updates through `CoupledState` owner callbacks.
Five updates, rejected-candidate state retention and reset pass on 100 and
the same binary on the already-running NGSolve A100. The generic diagonal
block does not include complete coil/body reaction assembly or nonlinear B-H.
The production persistent MEX gateway and Simulink connection remain pending.
See `results/esim-coupling-20261013.json` for residuals and memory limits.

## Persistent element-temperature context (2026-10-13)

`checkPersistentElementESIM.m` verifies the `element-esim-operators` frontend
in the production checked-handle gateway. Its explicit material tables use
kelvin and one row per element; impedance-map rows follow the eddy matrix CSR
storage order. `esim_cell_parameters` columns are node count, cell iteration
limit, cylinder flag, half thickness (m), angular frequency (rad/s), cell
tolerance and relaxation. Nonlinear temperature-dependent B-H inputs fail
explicitly until their gateway is accepted. Loss weights contain the existing
quadrature/one-half factors; the runtime multiplies them by real impedance.

The generic diagonal block passes five updates, candidate rejection with all
accepted states retained, reset, Simulink repeat, save/close/reopen and zero
live handles on 100 and the same native artifact on NGSolve A100. Official
MATLAB MCP execution and the official Toolkit structural check on 100 passed.
This does not certify complete reaction geometry assembly, physical validity,
or research release completion. That assembly must preserve the existing
SL/DL singular quadrature, mass inverse, surface-Poisson transpose work, coil
saddle constraint and total-field heat; incident-only coupling is insufficient.

## Temperature-dependent amplitude B-H research model

The dedicated coupled candidate accepts `bh_knots` (common H axis in A/m)
and `bh_values_T` (element x temperature x H knot, tesla). It builds checked
piecewise-linear B(H) coefficients in the same mexcuda module. Temperature
interpolation updates each element's resistivity, initial permeability and
amplitude law before its ESIM solve. No Python solver is invoked by MATLAB or
Simulink. Unsupported hysteresis history is not inferred from these tables.

Build an isolated candidate on the build host using
`Build.ps1 -IhCudaOnly -RequireNativeProvenance`. In MATLAB, add the resulting
artifact directory, `packages/radia-gpu/matlab` and this validation directory
to the path. Run `checkTemperatureBH(artifactDirectory, outputDirectory)`.
It creates `ih_temperature_bh_public.slx` with current/ambient input and
temperature, heat and separate real/imaginary element-Zs observation outputs. The model
workspace owns `ihConfig`; replace the synthetic material tables only with
measured, unit-checked data covering the operating temperature range.

Run the independent CPU equation check with
`python validation_test/ih_cuda/check_temperature_bh_cpu.py outputDirectory`.
Its original CPU finite-difference cell solver and dense block solve are
independent of the CUDA solution; geometry/BEM assembly is shared. This is
not independent full-physics certification or experimental validation.
The original public model uses a thin surface thermal layer, not resolved
bulk thermal FE. The released 0.1.3/0.1.4 wheels remain unchanged.

`TEMPERATURE_BH_EXPERIMENT.md` defines measurement items, staged validation,
uncertainty, held-out conditions and the distinction between amplitude B-H
and hysteresis loss. Raw experimental data and real-device cases remain
private. See the candidate result JSON for actual numerical acceptance;
documentation or a successful CUDA build alone does not establish it.

Set `observe_element_impedance=true` for the two optional Zs observation ports.
The default two-port temperature/heat interface remains available. Zs is the
last accepted electromagnetic step's impedance, evaluated at that step's
starting temperature; it is zero before the first accepted update and after
reset. Observation downloads occur only when requested. Intermediate heat and
temperature exchange between the electromagnetic and thermal stages remains
on the device. `element_impedance(handle)` likewise returns the accepted Zs,
including after a rejected candidate. Supply `solver_relative_tolerance` when
a stricter inner solve is needed for an ill-conditioned coupled block.

## Integration review checkpoint (2026-10-15)

The 2026-10-13 temperature-BH JSON is historical acceptance of the recorded
native artifact and MATLAB harness, not acceptance of a new wheel or current
main. The new package staging is opt-in via RADIA_IH_BUILD and requires a clean
exact-source receipt with the declared research features. Old fixed-operator
IH receipts cannot be used to populate the release payload.

A new wheel still needs the CI IH build/test lane, all-module installed
qualification including IH CPU equations and residuals, official Toolkit
save/reopen/repeat/failure checks on its installed sample, independent Claude
review and the owning publication gates. Independent full FEM and experimental
accuracy remain unverified under the research release scope; they are not
replaced by algebraic CPU parity.
