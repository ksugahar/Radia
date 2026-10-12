# Resident Matrix Build and HDiv-MMM CUDA subgates

The experimental mexcuda runtime retains symmetric HACApK leaf coefficients,
permutation, physical charge normalization and vector state on the GPU. It
applies a Gram operator or the charge-map composition `B^T G B` without vector
transfer or a host reduction between applications. Dense and low-rank upper
leaves mirror their transposes; diagonal leaves apply once. Deterministic row
reduction avoids cross-block atomics. General Gram application has four CUDA
launches, and the charge-map composition has six, independent of leaf count.

`fromHacapk` consumes native HACApK storage directly inside C++. MATLAB's numeric
leaf schema exercises the same packing adapter; native pointers never cross
that interface. MATLAB uses `radia.gpu.HMatrixOperator` for checked ownership.
The native module has no Python runtime or NGSolve TaskManager requirement.

## Bounded native Matrix Build

`assembly='analytic_constant_block'` constructs a block from flat row-major
coordinate column vectors `cells` (12 coordinates per linear tetrahedron) and
`faces` (9 per linear triangle). The prototype accepts at most 128 charges and
uses the existing all-analytic constant-density route: shared exact `PhiTet` /
`TriPotential`, the unchanged 64-point Gauss-Duffy volume rule and 7-point
Dunavant surface rule, symmetric off-diagonal evaluation, and charge
normalization. Coefficients remain on the GPU after generation. No CPU-generated
matrix or Python snapshot is needed for this entry point. Images, high-order,
curved geometry, alternative near/far choices and larger blocks are unsupported
and must not be substituted by centroid or point-cloud approximations.

This is a bounded block-generation frontend, not a replacement for HACApK's
cluster construction and ACA+ compression. Connecting these generated blocks
to full adaptive assembly, real NGSolve charge maps, constitutive terms,
preconditioners and the nonlinear HDiv-MMM solve remains required. FE direct
subsolves retain `sparsecholesky`; linear true-residual and physics gates remain
unchanged. Repeated operator application is not a linear solver and therefore
has no solver residual to certify.

## Reproduce

Build only on 100 from clean reviewed source:

```powershell
pwsh -NoProfile -File ./Build.ps1 -HMatrixCudaOnly -RequireNativeProvenance
python validation_test/hmatrix_cuda/validate_resident.py --artifact <artifact-dir> --out <new-job-dir>
python validation_test/hmatrix_cuda/measure_resident.py --artifact <artifact-dir> --inputs <numeric-job-dir>/inputs.mat --out <new-idle-job-dir>
```

Numerical checks run one below-normal MATLAB session and a one-thread,
caller-owned NGSolve TaskManager CPU reference. Synthetic dense/low-rank leaves,
uniform and unequal block sizes, permutations, normalization, dense/sparse
charge maps, real public pybind integrated-tetrahedron Gram evaluation and
native CUDA tetrahedron/triangle assembly are distinguished. Relative vector
differences must be <=1e-12. Invalid input, wrong/stale handles, reload, repeated
execution, wrapper cleanup and final live count are checked. Device reset
acceptance remains pending; do not reset a shared GPU to test it.

Timing checks queued/running CI, other computation, GPU and host budgets before
and after measurements. After warmup, two ABBA rounds compare rebuilding and
transferring on every application against one resident context for three
applications. Each continuous interval includes numeric leaf packing, context
construction, transfers, apply, recovery and release. It excludes geometry/FE
assembly and is an operator-residency diagnostic only. The first-payload labels
are after MATLAB/MEX initialization, not cold process/module startup.

100 timing is diagnostic, not formal idle-mdx performance acceptance. Exact
owned CUDA buffers, setup scratch and sampled free memory are distinguished;
no host/process/global allocation peak is claimed. Allocation checks leave
max(2 GiB, 20% device memory) free. A100 is never started or connected. The gateway rejects direct payloads over 2 GiB. The required HDF5 frontend is
not connected yet; these bounded fixtures are small.
The experimental class/module remain excluded from published wheel payloads
until their physical and release acceptance is complete.

## Recorded checkpoint (2026-10-12)

`results/resident-numerical-20261012.json` accepts ten small cases against
independent NumPy algebra and the installed public Python/pybind CPU route.
The latest native artifact is built from `0863b258a`; maximum relative vector
difference is 7.223e-16. Native geometry-to-CUDA constant-charge assembly passes
for tetrahedra and tetrahedron/triangle combinations. This is not full adaptive
Matrix Build or complete HDiv-MMM physical acceptance.

`results/residency-diagnostic-20261012.json` retains an earlier `b7be1da75`
artifact's guarded idle-100 diagnostic: a 275-dimensional synthetic mixed-leaf
operator, three applications, continuous packing/setup/transfers/recovery/release
median 7.718 ms rebuilding versus 3.696 ms resident (2.09x). Exact owned device
buffers are 573312 bytes in both modes; no buffer-memory reduction is claimed.
These times do not certify the latest artifact, a large physical problem, cold
process startup or mdx performance. A later timing attempt was rejected by the
other-computation guard, so latest-artifact timing remains pending.

`results/quadrature-sharing-20261012.json` verifies that the extracted shared
quadrature helper is textually identical to the previous CPU implementation.
The full CPU solver build remains a separate CI gate.

## Indexed assembly checkpoint

`hmatrix_entry_oracle.cuh` retains linear constant-charge geometry and its exact
quadrature on the GPU. It evaluates bounded indexed row/column requests without
allocating the global dense matrix. The native caller supplies a stream and
reuses one workspace; shared host callers serialize that workspace. Values are
raw physical entries, so the existing HACApK normalization remains the caller's
responsibility. The diagnostic `analytic_entry_batch` command uploads geometry
once per command and exercises three 256-entry-capacity batches; it is not a
persistent MATLAB assembly handle or an adaptive compression implementation.

`results/indexed-entry-20261012.json` records the clean `712fa3960` build,
600 requested entries with CPU relative difference 5.924e-16, and the ten
existing numerical/lifecycle cases. Full adaptive ACA+ callback wiring,
OpenMP leaf scheduling, real NGSolve charge-map assembly, material mass,
constraints, iterative solver and true residual remain required. Large-scale
assembly, whole-solve timing and peak-memory acceptance are pending. CI was
running and queued, so this checkpoint contains no new performance comparison.
