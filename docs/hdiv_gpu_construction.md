# HDiv resident GPU construction

Status: implementation design. The current public workflow still constructs
its charge Gram on the CPU and exports leaves. The private affine-TET input
preparation is a first seam, not device assembly or performance qualification.

## Ownership and supported first route

NGSolve owns HDiv spaces, orientation, Piola maps, geometry, FE assembly and
quadrature semantics. Radia owns the charge Green operator, compression and
resident material solve. The first device-construction route targets affine
pure tetrahedra, HDiv orders 1 and 2, without live mesh deformation or images.
Unsupported geometry must raise before device allocation; no CPU construction
fallback is permitted for a selected device-build route. Keep the current
DOF, charge, element, workspace and reserve guards.

The receipt-checked mesh snapshot remains the sole input to workflow-owned
FE staging. Do not reopen mutable caller files after checking their bytes.

## Preserve the actual CPU entry contract

`_flat_tet_charge_inputs` preserves charge numbering, the NG charge map,
geometry, and separate near, inner and low-order far rules from the existing
CPU builder. The private GPU preparation uses the current workflow defaults;
it must not silently substitute the nonlinear builder's different defaults.

For supported affine polynomial charges, CPU `PhiInner` selects analytic
moments for volume degree <= 2 and surface degree <= 3. The near entry integrates
that potential over the target rule. The far entry uses the existing low-rule
double quadrature and geometric separation criterion. Preserve ordered symmetric
averaging and diagonal normalization. The existing curved CUDA entry kernel
is therefore not an unconditional replacement: it has no equivalent flat
near/far dispatch. Merely carrying inner-rule arrays does not mean this CPU
route executes a Duffy/subtraction rule for these supported polynomial degrees.

Reuse the shared affine analytic moment kernels after checking vertex/reference
conventions against the CPU route. Audit mapping, target weights, far product,
separation and symmetry independently. Fixed constant-charge assembly is not
a substitute for the order-1/2 polynomial charge basis.

## Device data path

1. Stage NG-owned charge mapping and geometry once, retaining original numbering.
2. Generate the HACApK block frame and permutation without evaluating coefficients.
   Existing pre-fill admission borrows unfilled leaves; a durable frame export
   must also bind the control/permutation lifetime and copy metadata safely.
3. Upload geometry, rules, sparse maps and frame metadata once. Compute diagonal
   normalization and matrix coefficients on device.
4. Construct compressed factors on device. Keep pivot rows, columns and factors
   resident. CPU-side ACA calling a GPU entry oracle and receiving coefficient
   arrays is an intermediate route, not resident GPU compression.
5. Retain coefficients, factors, reference mass, material maps, preconditioner and
   solver workspace across nonlinear iterations. Transfer convergence scalars
   and explicitly requested final outputs rather than repeated full states.

A new explicit device constructor must consume a validated owned payload,
fail before publishing partially initialized owners, and support deterministic
close and exception cleanup. It must not reuse the leaf-upload constructor by
first creating CPU coefficient arrays.

## Shared compilation boundary

Compile Radia-owned CUDA with nvcc and package the same compiled objects in a
CMake STATIC numerical archive. Python/pybind11 is the first development and
qualification gateway. A thin plain MEX gateway can link the same archive;
mexcuda is not a fallback. NG staging remains host C++, outside the CUDA unit.

Keep creation/destruction and diagnostics inside each linked gateway module;
do not exchange native handles between gateways. Static linking shares the
compiled numerical implementation, not process state. A runtime DLL and generic
opaque ABI are not required by this design. The current child adds the static
archive and provenance seam only; a MEX consumer and device constructor are
still separate work, not certified here.

Bind exact source, compiler/CRT, architectures and object/archive hashes. The
archive is a build artifact, not a wheel payload or a new runtime dependency.
Missing/invalid artifacts fail explicitly. Preserve the existing runtime,
packaging, fault/lifecycle and numerical checks.

## Bundled Windows CUDA SDK layout

The selected MATLAB-bundled nvcc uses its own CICC/profile configuration, while
libdevice is located in the MATLAB binary directory. Pass that directory
explicitly using nvcc's documented `--libdevice-directory` option. Require the
bundled compiler, profile, CICC, libdevice and matching-version CUDA header;
record their paths/hashes privately. Do not select another toolkit via PATH or
environment fallback. The profile is retained; a broader profile bypass requires
separate compile evidence. Source fixtures do not certify this invocation.

NVIDIA documents the option in its [CUDA 12.8 compiler driver manual](https://docs.nvidia.com/cuda/archive/12.8.0/cuda-compiler-driver-nvcc/index.html).

## Memory accounting and measurements

Admit the maximum simultaneously live bytes, including geometry/rules,
coefficient/factor capacity, sparse maps and transposes, compression scratch,
preconditioner, solver vectors, and library/graph resources. Account separately
for construction and solve peaks and preserve the existing device reserve.
A matrix fitting in memory does not imply the complete workflow fits or is fast.

Record host FE preparation, H2D input, GPU construction, compression,
preconditioner, solve, D2H result, and end-to-end times. Use stream events for
GPU stages and wall time for complete workflows; distinguish cold/warm runs.
Record transfer bytes/counts, live/peak/closed owned device bytes and device-wide
free-memory observations separately. Do not sum overlapping intervals or
attribute another process's allocations to this owner.

## Qualification

Compare against the independent CPU entry/solve route using identical inputs
and existing adopted numerical criteria. Include near/self, separated pairs,
polynomial order, permutations and geometry rejection controls. Inject failures
at construction/allocation/publication and verify cleanup and no CPU fallback.
Keep Python/MEX same-kernel bitwise parity distinct from CPU/GPU tolerance checks.
Source tests or import/compile success do not certify numerical GPU acceptance.
Scaling and speed require measured end-to-end evidence; timings are diagnostic.
