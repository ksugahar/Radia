# Mixed scalar potential backend probes

Overall formulation parity remains **OPEN**. Per-run pass flags apply only to
the declared stencil and gate; they do not establish whole-interface equivalence
or general solver qualification.

These executable investigations preserve the matching-trace Python formulation
while exercising the optional direct pybind11 CSR adapter in an isolated
NGSolve 6.2.2608 environment. They do not qualify a complete native package or
claim performance of a different formulation or preconditioner.

`mixed_omega_backend_probe.py --out result.json --order 2 --maxh .14` assembles
the production mixed total/reduced operator with nonzero physical boundary
data and relative permeabilities 100/1. The production CPU CG/local solution is
the field reference. The caller eliminates constrained coefficients and passes
the identical free-row SPD system to explicit CPU sparsecholesky and GPU
Jacobi-PCG selections. GPU timings include CSR preparation, transfer, solve,
download and the adapter's original-system residual check. The first call and
subsequent calls are separate. CPU CG phase timings and CPU direct timings are
reported separately; direct versus iterative timings are not interchangeable.
The optional NPZ is solely independent cross-version validation evidence.

`mixed_omega_cpu_regression.py --input result.npz --out old-runtime.json`
checks that saved system in a separately selected CPU runtime. It does not
test cross-version FE assembly or regenerate the mesh.

`mixed_omega_kelvin_sphere_probe.py --out sphere.json --maxh .16 --order 3
--geometry-order 3` studies a radius 0.3 m permeable sphere in a unit applied
field, with a 1 m Kelvin radius and a translated image ball. The applied
uniform field is analytic; only the decaying reaction is solved in the Kelvin
ball. Removing the applied source's outer normal flux from the physical
volume load is essential for that split. The declared five-point field stencil
is compared with the analytic permeable-sphere solution with a 1% relative
RMS gate. This does not certify all interface traces or arbitrary source fields.
Mesh, polynomial order and geometry order are independent inputs. Earlier
coarse failures are retained; an initial exploratory 5% gate is superseded by
the current 1% gate. In particular a 2.14% result is not accepted under the
current gate.

`mixed_omega_mmm_sphere_probe.py --runtime-source SOURCE/src --kelvin-result
sphere.json --out mmm.json` loads an explicitly supplied matching native/source
runtime and independently solves the same sphere using HDiv-MMM. It compares
average magnetization against the analytic solution and the Kelvin stencil's
interior magnetization, with 1% gates and an original-operator residual check.
It supplies the applied field to both routes: feeding the iron's already
demagnetized total field to MMM would count self-interaction twice. This is
an independent shared-physics cross-check, not an implemented coupled
optimization loop. MMM native provenance remains separate from CSR provenance.

All probes save machine-readable JSON. Keep operational paths, native runtime
provenance, raw failures and diagnostic timings outside published evidence.
Shared-host timings are diagnostic, not a controlled performance certificate.

## Two-way separated-body investigation

`mixed_omega_kelvin_coupled_probe.py` implements a bounded, linear two-body
investigation. A radius 0.2 m FEM sphere is inside a radius 0.5 m physical
Kelvin ball. An independently meshed MMM sphere is centered at (0,0,0.85) m;
its complete volume must remain outside that physical ball. Both relative
permeabilities are 100. The global unit applied field is supplied once.

The MMM field is harmonic throughout the physical FEM ball. Consequently its
volume load minus its outer flux is exactly its inner interface flux. The
Kelvin image represents only the FEM reaction, not the remote MMM source.
The periodic scalar H1 unknown is the reduced potential in air/image; adding
the projected source lift in iron reconstructs the total potential. Its SPD
operator can be compared with the existing multiplier formulation in the
separate parity mode; coupled and repeated modes do not perform this comparison.
NGSolve owns the periodic space, geometry, trace bases, assembly and source
projection. Quadrature sensitivity is demonstrated and remains under
investigation; the declared quadrature is not a proof that the near-image
metric is resolved.

Gauss-Seidel iteration sends the MMM source field to FEM and only the mapped
FEM reaction back to MMM. A Kelvin inverse/Householder transform evaluates the
reaction outside the physical ball. A declared voxel transfer is independently
checked at 26 points. Self-interaction is never sent back as applied field.
An unchanged small update is insufficient: the converged MMM coefficients are
checked against a fresh applied field and the original mass-plus-demagnetizing
operator for the voxel-approximated external field. This residual does not
certify the exact transferred coupled field; the voxel holdout error gate is
separate. The optional `--crosscheck-mmm` uses a genuinely separate two-body
all-MMM solution and checks both bodies and original residuals.

Source choices are explicit. `--harmonic-source-degree 8` fits a local harmonic
potential using native source samples and checks separate held-out samples at
an absolute RMS gate of 1e-5 A/m. It is an additional controlled approximation,
not exact source evaluation or a whole-volume error bound. Lower-order
failures are retained. `--interpolate-source` instead uses the existing P2
NGSolve source interpolant with an explicit native batch evaluator. Neither
choice silently replaces a failed exact request. Their gates are separate from
linear residuals and the 1% physical comparison.

`--reuse` retains the CPU factor or prepared host CSR arrays. GPU transfers
still occur on every call. `--reuse-surface-load` reuses constant/mass data for
this fixed, single connected interface, with the same flux compatibility gate;
it is not a general multi-component cache API. The MMM geometry operator is
retained through public `HDivSolver` calls. Explicit CPU CG choices include
local, H1-AMG and BDDC/H1-AMG, with true original-system residual checks.

The GPU selection applies to the scalar FEM solve. MMM and source preparation
remain explicitly CPU work. Reports distinguish scalar transfer-inclusive
solve time, initial preparation, complete coupling time and independent
reference-validation time. This is solver optimization and coupling validation,
not design/topology optimization or a fully device-resident coupled solver.


CPU Krylov preparation now keeps the original NGSolve operator instead of
exporting and rebuilding an unused compact matrix. Every backend checks the
free-equation residual against that original operator. The compact export is
still used for explicit direct/GPU routes.

The coarse P2 parity PASS (1353 free DOFs, q20, about 6.95e-9 A/m)
is not a baseline for the P3 FAIL (4425 free DOFs, about 1.58e-7 A/m).
With the historical independently integrated lift load, P3 had not passed the
1e-8 A/m gate; no native-copy regression was established.
P2 q10 already failed at about 1.1e-7 A/m. Quadrature sensitivity is therefore
open, and the absolute gate was retained from the P2 q20 experiment rather
than derived from an error bound. Analytical accuracy, unchanged-operator
agreement and formulation parity are separate acceptance signals.

GPU dependencies are imported only with `--backend gpu`. Supply any required
GPU DLL search directories explicitly with repeated `--gpu-dll-directory`;
CPU native runtime dependencies can be supplied with `--native-dll-directory`.
No installation path is assumed. The coupling probe's `--bonus-intorder`
controls both compared assemblies; each sweep run records its order and q.

The direct-solver quadrature sweep at q=10/20/30/40 reproduces the P2 trend,
while P3 differences are nonmonotonic (about 1.41e-8, 1.58e-7, 4.11e-8,
2.03e-7 A/m). At q40 the copied interface trace has zero L2 difference;
removing mean-flux correction changes the maximum difference by less than
2e-16 A/m. Neither trace copying nor that compatibility correction explains
this discrepancy in the tested independently integrated-load case. Those runs
remain failed; increasing quadrature alone did not establish equivalence. Every run retains its original
pass/fail gate and separate JSON; no threshold was relaxed.

For linear-algebra diagnostics, `--refinement-steps` retains the full assembled
reference residual, including multiplier rows, and applies corrections through
the existing matching-trace inverse. This inverse eliminates constraints and
factors only trace-mass and constrained primal SPD matrices; it does not apply
Cholesky to the indefinite block matrix. `--condition-estimate` uses one-norm
estimates and the implemented inverse; these are not certified condition bounds.
On the coarse P3 q20 case, refinement reduces the full free residual from about
1.46e-12 to 5.53e-14 without reducing the field discrepancy. Estimated unscaled
conditioning alone therefore does not establish the cause. Both fields use the
same mesh point evaluator and physical reconstruction convention.


## P3 parity resolved for the declared stencil

Mechanism: load/operator quadrature consistency.

The remaining coarse P3 discrepancy is localized to the independently
integrated iron lift load. On the same mesh the condensed reference matrix
matches the common-H1 matrix to about 1e-16, and the copied lift satisfies the
multiplier constraints to roundoff. The right-hand sides differ by about
4.24e-10 in Euclidean norm. Building the lift contribution as
`-A_iron * lift`, with the same assembled iron operator used by the solve,
reduces the reference right-hand-side difference to about 1e-21 and the
four-point maximum field difference from 1.58e-7 to 3.76e-14 A/m. The original
1e-8 A/m gate is unchanged. Both comparisons use the same in-memory mesh,
source, projection and point-evaluation operator.

An additional check gives both the linear and bilinear lift expressions the
same explicit NGSolve integration rule. Their relative differences then fall
to approximately 1e-15. Equal `bonus_intorder` alone had not guaranteed an
identical discrete lift integral in the two forms. This is a load/operator
consistency issue in this probe, not evidence of defective interface matching
or an indefinite Cholesky factorization.

The probe now defaults to `--lift-load-route assembled`; the independently
integrated route remains an explicit attribution experiment. `--operator-audit`
records matrix, load and constraint defects and a paired solve on one mesh.
`--quadrature-audit` additionally compares shared explicit integration rules.
The declared P3 stencil passes for direct, BDDC/H1-AMG and explicit GPU solves;
this does not establish a whole-interface error theorem, analytical accuracy
on all meshes, or a new timing certificate. Earlier failures remain preserved.


The production matching-trace driver's reduced lift and the general driver's
Kelvin lift also use original-operator actions. Curved P2/P3 regression tests
pin their assembled right-hand sides against the corresponding stiffness
contributions; both old implementations fail these tests. Physical source
fields, Hodge-projection data, prescribed fluxes and weak jump data retain
their own integrals and accuracy gates. They are not silently interpolated or
replaced by a discrete lift. Reassembling the regional lift operator is counted
in RHS preparation, so previous preparation timings do not qualify this update.


The separated-body adapter accepts `--threads` from 1 through 16 and records
the requested count in each JSON result. The default remains one thread.
Thread affinity and topology are supplied and recorded by the external run
supervisor; selecting a thread count alone does not establish affinity.



## Experimental resident algebraic AMG through Python

`Build.ps1 -CudaMexOnly -ResidentAmgPython` stages an optional pybind11
module using the admitted MATLAB-bundled CUDA compiler and supported MSVC.
It does not install the module or change a running backend. The module has
no MATLAB Engine or NGSolve runtime calls; the caller owns FE assembly.

The standalone validation adapter is `mixed_omega_amg_resident.ResidentAMG`.
Pass `directory=DIR` explicitly; it never changes the production backend.
PyAMG constructs a smoothed-aggregation Galerkin hierarchy on the CPU from
the assembled free system and an NGSolve-interpolated constant candidate.
A bounded coarse inverse action is prepared using NGSolve sparsecholesky.
The validation adapter owns a Python TaskManager region only for this coarse
CPU preparation; FE assembly remains in its caller-owned region. Device
construction and execution occur outside this adapter-owned CPU region.
Matrices, transfers, coarse action, CG vectors and level workspaces then
remain on the GPU. Symmetric pre/post damped-Jacobi smoothing and the coarse
action form a fixed V-cycle. CG batches eight iterations in a CUDA Graph;
the original-system residual is checked independently on the host as well.
This route is smoothed aggregation, not a GPU port of h1amg or AMS.

The `strength` and `damping` constructor keyword arguments specify separate experimental
configurations. The damping is scaled by a certified normalized row-sum
spectral bound and must lie strictly between zero and two. Unsupported or
singular coarse hierarchies fail; there is no automatic configuration or
backend fallback. CPU references remain independent of these optional imports.

Preparation time includes hierarchy construction, coarse preparation and
device upload. Solve time includes right-hand-side upload and solution
download. The process-level bounded-job receipt measures startup-inclusive
time. Setup amortization and full workflow time must accompany a warm solve
comparison. Run `resident_amg_contract_probe.py` to check repeated changing
right-hand sides, zero input, nonfinite rejection and released contexts.



HDF5 transfer reuse is explicit: `cache_path=PATH` with `cache_mode="build"`,
`"reuse"` (fresh Galerkin operators) or `"reuse-hierarchy"` (unchanged fine
operator within the recorded roundoff guard; retained coarse blocks serve
only as a preconditioner). Candidate coefficients, settings, sparse structure
and all payload hashes are checked. Existing files are never overwritten.
The arrays are algebraic coefficients in the caller's consistent SI system;
no coordinate conversion occurs. `resident_amg_cache_probe.py --out result.json`
checks this contract without CUDA or NGSolve imports.

For changing SPD operators, `nonlinear_cache_path=PATH` stores interpolation
only. `nonlinear_cache_mode="reuse"` checks the sparse layout and candidate;
all fine and Galerkin operator values are current. This does not certify
nonlinear FEM or material coupling. `nonlinear_amg_cache_probe.py` checks the
updated-operator solves using explicit native and CUDA DLL directories.

The native PCG always checks the original resident system. If recursive
convergence leaves an excessive original residual, at most three resident
corrections use the same requested tolerance and the same total iteration
budget. The requested PCG tolerance is distinct from the native original-system
certificate: `min(max(10*tolerance, 1e-12), 1e-6)`, returned as `certificate_gate`.
The validation adapter separately evaluates `||Ax-b||/||b||` with SciPy on
its retained original fine matrix, raises above `1e-6` (or on nonfinite
results), and records `host_true_relative_residual`, `host_certificate_gate`
and verification time. For a zero RHS the absolute residual is checked.
There is no CPU solve fallback. The
contract probe includes independently solved SPD systems, scaled right-hand
sides, repeated solves and lifecycle failures. It records any corrections
in its JSON result. No device build artifact is distributed by these probes.

The standalone native module requires compute capability 7.5 or newer and
rejects older devices before allocating the resident context.
