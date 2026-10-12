# RadFld CUDA validation and next acceleration boundaries

The first experimental CUDA backend implements **the existing retained point
cloud representation** of the HDiv RadFld evaluator. It does not replace
analytic polynomial or curved-element sources with sampled points, implement
generic `rad.Fld` object dispatch, or claim that all RadFld source families are
accelerated. Analytic/curved evaluators reject cloud export explicitly.

Both CPU and CUDA share the point-source, compensated-sum and quadrupole
kernels. The existing CPU tree constructor supplies exactly its source order,
bounding boxes and moments. CUDA follows the same opening criterion, leaf
order, and mirror/rotation images through a threaded traversal without a fixed
device stack. Both `direct` and `tree` are explicit choices; there is no hidden
CPU evaluation fallback or automatic change of approximation.

`radia.gpu.FieldCloud` takes numeric source positions and strengths. Its result
uses the existing cloud convention `sum q*(r-x)/|r-x|^3`, without `1/(4*pi)`;
the caller retains the original physical conversion. The source geometry and
tree stay on the device between calls. No Python adapter, Python-generated
solver snapshot, NGSolve DLL or TaskManager region is required at runtime.
This is a staged experiment, excluded from wheel staging until its own
artifact/provenance and release gates are added. Published runtimes are unchanged.

Build on 100 in an isolated clean checkout with `Build.ps1 -FieldCudaOnly
-RequireNativeProvenance`. MATLAB-supported MSVC/CUDA configurations compile
sm_75, sm_80 and compute_80 PTX. CUDA multiply-add contraction is disabled;
there is no fast-math override. Exact device buffer counts are admitted only
when max(2 GiB, 20% total GPU memory) remains free. A100 stays off and is never
connected to by this workflow.

## Independent CPU reference and measurements

`build_reference.py --out <job-reference-dir>` compiles a validation-only
pybind11 adapter against **unmodified pre-CUDA CPU sources** recovered from the
explicit Git revision in its receipt. This is an independent CPU route using
NGSolve scheduling, not an input prerequisite for MATLAB users. No reference
binary is installed into another environment.

`validate.py --reference-dir <job-reference-dir> --out <job-results-dir>` uses
HDF5 inputs at every size (including inputs beyond 2 GB). It checks
CI/preflight queues, other numerical/build processes and GPU load before
starting MATLAB and between every sample. A conflict stops the run, closes
only the task-owned Engine session and saves a failed/pending JSON; it never
kills another job. Default sizes are 4096, 16384 and 65536 sources with 1024
targets. Near, far, coincident-cloud, signed cancellation and image rotation
cases are included. Increase sizes only within device and host memory budgets.

Warm-up precedes two ABBA rounds for each algorithm, comparing original CPU
construction/evaluation/result recovery against MATLAB source/tree setup,
upload, CUDA evaluation, download, destruction and synchronization. Input
reading and numerical checks are outside the timing. The public API completes
within MATLAB. Startup time is recorded separately. CPU threads and all
source, adapter, binary and input hashes are retained.

GPU versus CPU **for the same algorithm** must meet relative vector and
maximum pointwise error <= 1e-12. The tree also must meet <= 1e-5 versus CPU
direct, with an explicitly recorded near-zero floor for pointwise ratios.
Field comparisons use vector differences, never differences of magnitudes.
Failure stops acceptance and retains the JSON. No linear solve occurs here,
so a solver true residual is not applicable; HDiv-MMM solver validation still
requires its existing true-residual and physical gates.

JSON records exact simultaneously owned CUDA bytes, sampled CUDA free memory,
MATLAB available/host memory, and Python process working-set/private bytes.
These are not a whole-process native allocation peak; Python's reported peak
is process-lifetime, not per-case. 100/RTX 5000 timing is diagnostic provenance,
never formal idle-mdx performance acceptance. Raw MAT files, inputs, logs,
compiler receipts and operational identities remain in private durable storage;
only reviewed public JSON belongs in Git.

## Following work authorized with this acceleration program

1. Complete the staged cloud numerical/memory measurements after 100 is idle.
   Connect supported RadFld source packing explicitly; do not present this
   numeric backend as complete public RadFld parity. Port analytic/high-order
   and curved near-field kernels without changing their physical contracts.
2. Matrix Build: profile the existing `BuildHMatrix` and charge-basis timings.
   Batch well-separated entry/block evaluation for CUDA while preserving
   HACApK compression, symmetric host blocks, images, singular/Duffy near
   integrals, and NGSolve-owned geometry/orientation/quadrature. Accelerating
   a sampled-Laplace test strategy alone does not certify HDiv-MMM.
3. HDiv-MMM: profile repeated H-matrix matvec, constitutive updates and solve
   orchestration. Keep reusable data resident and measure transfer-inclusive
   application time. Preserve NGSolve `sparsecholesky` for direct FE subsolves,
   nonlinear convergence and true residual <= 1e-6. TaskManager is not started
   in mexcuda. Record numerical, physical and memory evidence separately.
4. A RadFld-specific FMM remains a candidate after direct/tree measurements.
   Core compact interactions retain HACApK; do not restore core ExaFMM.
   NGSolve native multipoles/FMM may serve the policy-approved smooth
   free-space BEM/Biot-Savart scope. No GPU FMM backend is claimed here.

## Analytic and curved CUDA migration in progress

The production scalar formulas now live in `rad_hdiv_field_math.inc`, instantiated
with the existing degree-18 CPU workspace and a degree-6 CUDA workspace.
`rad_field_source_kernel.inc` shares exact source dispatch, origin translation,
and the existing signed closest-reference Duffy integration. Far-field quadrature
tables are generated by the original CPU routine and uploaded unchanged.
`ExportSnapshot` retains analytic, high-order and curved source records; the
cloud-only export still rejects element sources. CPU pybind11 scheduling remains
caller-owned TaskManager; mexcuda uses CUDA or OpenMP, without starting TaskManager.

This migration is **not yet built or numerically accepted**. The native `create_elements` command and MATLAB `FieldElements` wrapper now
accept a version-1 numeric struct with `representation`, `volume`, and `surface`.
Representations are `tet` (16/22 columns), `polynomial_tet` (32/22), `cubic_tet`
(32/29), `sextic_quartic_tet` (96/44), and `curved_tet` (34/24). Affine origins
are optional Nx3 arrays. Curved sources require Gauss points/weights as column
vectors. Records follow the production CPU header exactly; this numeric API does
not assemble NGSolve finite elements. Both experimental wrappers remain excluded
from wheel staging pending release acceptance.

Run `validate.py --representation <kind> --sizes 4 16 --targets 16 ...` separately
for every element family before scaling. The independent pre-CUDA oracle adapter
must first be rebuilt into a new job-owned directory. Original synthetic cases
exercise materially scaled high-order terms and actual non-affine P2 geometry,
near/interior/coplanar/far targets, and image rotations. The same strict parity
limits apply initially; no element family has yet passed GPU validation. Distribution
contract tests do not certify these mathematical changes. Rebuild and compare
against the pre-migration CPU reference before any release integration. The
session deliverable is a radia-gpu release after the requested RadFld, Matrix
Build and HDiv-MMM work and their applicable release gates have been completed.

## Required IH Simulink release comparison

The release scope also includes the radia-ih GPU implementation running in
Simulink, compared against an independent Python/pybind11 CPU reference for the
same publishable physical IH case. A one-DOF smoke configuration, a frozen unit-
current loss field, or GPU PCG alone does not certify the coupled IH workflow.
The existing IH_RELEASE.md preview limitations must be checked against the final
implementation rather than assumed resolved.

Record the supported electromagnetic/thermal physics, mesh/label and material
contracts, excitation and time-step history, numerical limits, true residuals,
loss/power and temperature comparisons, cold start and warmed transfer-inclusive
full execution time, memory indicators and limitations, and exact source/input/
artifact identities in JSON. Use identical supported physical assumptions on
both routes; do not compare nonlinear CPU physics with a frozen GPU model.
Python is an independent accuracy oracle, not a GPU-user setup dependency.

Apply simulink-app-health using the official MATLAB MCP and Simulink toolkit.
Check the actual model, save it, close and reopen the exact saved path, repeat
execution, and verify failure cleanup. Do not edit SLX ZIP/XML. Integrate the
latest accepted IH work without overwriting another session's changes. CI and
same-artifact distribution acceptance remain release gates. A100 remains off.

The IH implementation now starts with resident complex CSR matvec, scalar-drive
RHS generation, and a coherent quadratic loss projection in
`packages/radia-gpu/native/ih_resident_kernels.cuh`. These are unbuilt primitives,
not a working coupled IH backend. The loss projection sums complex basis fields
before squaring, retaining cross terms; volume and surface measures must remain
explicit and distinct. Next connect checked operator/state ownership, the
complex solve and thermal update, residual gates and a shared Simulink context.
Upload static operators once, exchange only required controls/observations at
step boundaries, and avoid downloading heat/temperature solely to re-upload it
between the coupled device computations. Preserve Simulink output/update order.

## Native FEM projection bridge checkpoint

`Build.ps1 -ProjectionCudaOnly -RequireNativeProvenance` now builds an
experimental mexcuda gateway that reuses NGSolve SetValues. Host NGSolve code
is compiled by supported MSVC (C++20, matching AVX2 ABI); mexcuda links the
entry point. This compilation alone is not GPU acceleration.

The context captures the actual mapped evaluation targets through an initial
SetValues pass, then accepts bulk field values as a native vector CF and uses
SetValues again. `projectBatchedField.m` connects MATLAB/CUDA field callbacks
and exact-target/source-revision reuse. No nodal-sampling substitute is used.
The initial capture pass and source construction must be included in cold
whole-workflow timing; the helper's timer alone excludes those construction costs.
NGSolve owns all FE interpolation, orientation and element transformations.
No TaskManager is started. Only real 3D HDiv/HCurl/VectorH1 are supported now.

`check_projection_bridge.py` independently constructs Python/NGSolve analytic
Set references and optionally compares retained-cloud CPU pybind direct fields
with the identified earlier CUDA cloud artifact through the same native
projection. Python is validation-only; MATLAB performs the batched GPU route.
`results/native-projection-cloud-20261012.json` records the clean projection
build and the original cloud artifact/reference hashes. Three order-1 spaces
passed, with analytic coefficient-vector error below 9e-16 and zero observed
CPU/CUDA cloud projection difference for this small case. This is not physical
FE L2 acceptance or general Radia-object/analytic/high-order/curved acceptance.
Full legacy command migration, profiling counters, lifecycle across clear/reload,
source-strength update tests, physical field gates and idle timing/memory are
still pending. No speedup or complete RadFld migration is claimed.
