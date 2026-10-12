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

The shared ACA active-vector extraction now has an opt-in compile-time
`RAD_HACAPK_BATCH_ENTRIES` boundary. Default CPU builds still call the scalar
entry function in the same order. Opt-in builds request active entries in
chunks of at most 256, scatter to their original positions, and then execute
the existing ACA subtraction/pivot logic. No global callback setter was added.
`check_entry_batch.py` builds and executes both C paths on 100;
`results/aca-entry-contract-20261012.json` records 36 cases per path covering
both directions, masks, permutation and chunk boundaries. This tests extraction
only: the CUDA callback definition, complete native HACApK build and compression
parity remain pending. The GPU build must not enable the macro before that
definition and its physical-entry normalization/error handling are connected.

## CUDA ACA connection: numerical stop

The module now defines the C batch callback with a scoped thread-local binding
to the resident analytic oracle, captures exceptions at the C ABI, and rethrows
before accepting factors. `analytic_aca_vector` passes eight masked/permuted
row/column cases with CPU relative difference 9.688e-16. The original
`calc_vec`, `aca` and `acaplus` bodies are shared verbatim through
`cHACApK_aca.inc`; the narrow native C host unit needs no TaskManager or Python.
Host ACA arithmetic uses the existing non-BLAS implementation intentionally.

Compression is **not accepted**. Clean build `8a3885533` compiles and links the
existing ACA+ to CUDA entries, with charge-diagonal normalization and the
standard 1e-3 threshold multiplier. The two-by-two offdiagonal test returns
rank one and relative error 0.0594196. Independent scalar C execution of the
same shared ACA+ reproduces that rank/error (`check_aca_cpu.c` consumes the
public toy Gram values from the validation MAT). This isolates the failure
from CUDA entry parity; it does not certify compression or justify changing
the mathematical route. See `results/aca-cuda-connection-20261012.json`.
Large assembly and solve tests stop here. Next check actual admissible-leaf
selection and small-block ACA+ behavior, and complete the shared CPU native
build gate. No new whole-time comparison was performed with CI queued.

## Cyclic-search correction

The original reference searches stopped before the final cyclic row/column.
For the two-by-two fixture that skipped the only remaining active candidate,
leaving zero reference residuals and accepting rank one too early. Bounded
full-cycle scans now include that candidate; retry limits, masks, ACA updates
and tolerances are preserved in the shared CPU/CUDA body. Twelve scalar and
twelve batched C regression cases cover square/rectangular full-rank and
rank-one matrices (`check_entry_batch.py --compression`).

Clean CUDA build `af2f12ad6` returns rank two for the physical tiny fixture,
CPU relative difference 1.403e-16; independent scalar C gives 1.235e-16.
`validate_resident.py --aca-block-size 66` also checks a separated 66-by-66
constant-tet block against independent pybind physical entries. The actual
CPU builder reports two low-rank and eight dense leaves for 132 charges.
Physical block error is 2.046e-9, but rank reaches the diagnostic cap of 32.
That is not convergence to requested 1e-10. The full builder must retain its
existing exact-dense promotion when a cap is exhausted below both dimensions.
`results/aca-cyclic-search-fix-20261012.json` distinguishes these gates.

In the uniform one-charge-per-element path, leaf_size=64 yields threshold 65;
both dimensions must reach it and support-box gap admissibility must pass.
The earlier forced 2-by-2 ACA test is normally a dense leaf, but still exposed
a shared algorithm defect worth fixing. Variable-DOF leaf rules require their
own acceptance. Full adaptive assembly, full CPU native build, HDiv-MMM solves,
true residuals and whole-run performance remain pending. This supersedes the
cyclic-search failure, without deleting its original evidence.
