# Experimental AMS CUDA MEX

This experiment applies the existing real FP64 Radia AMS preconditioner on a
CUDA GPU. It is an assembled-operator bridge, not NGSolve's upstream
matrix-free GPU backend. It is not a released MATLAB solver or a claim of
production AMS GPU support.

Python/NGSolve owns the HCurl space, numbering, discrete gradient, assembly
and boundary conditions. The shared `HypreBasedAMS`/`CompactAMG` C++ code
owns interpolation, PMIS coarsening, Galerkin products and coarse
factorization. A validation-only Python extension copies that existing
numeric hierarchy and invokes its actual CPU `Mult` for probe comparisons.
There is no second MATLAB CPU solver implementation.

CUDA MEX owns the FP64 AMS cycle, l1-Jacobi sweeps, residual products,
restriction, prolongation, auxiliary V-cycles and PCG vector operations.
Both real AMS cycle types 1 and 7 are supported. Mixed precision and
non-CompactAMG auxiliary solvers fail explicitly. Coarse sparse solves use
NGSolve's existing ordered sparse LDL factors on the GPU; row products
run within one warp while triangular dependencies stay ordered. Existing dense coarse inverses and the
explicit ten-sweep coarse treatment are preserved when supplied by setup.
Coarse factorization and hierarchy creation remain CPU setup work. The
host computes PCG scalar control from cuBLAS GPU dot-product results; no
CPU vector or coarse solve is used as a fallback.

The gateway takes a versioned struct containing GPU arrays, a GPU RHS,
tolerance, maximum iteration count and a logical apply-only flag. It has
no persistent handles. Every call validates arrays and allocates temporary
workspace; timings include that overhead. Reusing the input hierarchy on
the GPU therefore does not mean reusing all gateway workspace. This is a
correctness-first prototype and can be slower than CPU AMS.

Build on the designated Windows build host in an isolated checkout with the configured MATLAB
Engine Python environment, MATLAB R2026a and its compatible VS2022 toolchain:

```powershell
./Build.ps1 -AmsCudaOnly -RequireNativeProvenance
python validation_test/ams_cuda/run.py --out "$env:TEMP/ams-cuda-owned-run" --maxh 0.5 0.2 0.12
```

The CUDA binary includes sm75, sm80 and compute80 PTX. A different OS needs
its own MEX build. This initial build script targets the confirmed Windows
build environment. A100 performance is unverified until separately measured
on an idle compute GPU; no A100 result is inferred from RTX measurements.

Acceptance compares random preconditioner applications with shared C++ AMS,
and the resulting solution with an independent NGSolve `sparsecholesky`
solve of the same manufactured HCurl curl-curl plus mass problem. The
FP64 true residual must be at most 1e-6; this small run requests 1e-9.
Results include iterations, first-call and warmed times, input upload,
transfer-inclusive time, GPU memory observations, versions and hashes.
New CPU reference timings use serial NGSolve without TaskManager. Historical
JSON retains the threading configuration used when it was measured. Invalid CSR, CPU input, NaN,
invalid smoother, nonconvergence, zero RHS and clear/reload checks run on
the actual gateway. Only sanitized JSON evidence belongs in Git; generated
HDF5 inputs, raw logs and binaries stay outside tracked source.

The first gateway call is measured after MATLAB and `gpuDevice`
initialization; MATLAB startup is reported separately. Available-memory
deltas and live-workspace samples are observations, not allocator peaks.
CPU CG stops on its preconditioned norm and counts the initial residual
check; GPU PCG stops on the Euclidean norm and counts updates. Both undergo
the independent true-residual gate, but iteration counts are not directly
interchangeable.

## RTX 5000 result, 2026-10-11

The launch-reduction candidate fuses CSR Jacobi products and updates with
separate old/new buffers, fuses prolongation and addition, and combines PCG
solution/residual updates. GPU-to-GPU copies are queued on the same default
stream as kernels and cuBLAS. Device-to-device cudaMemcpy is already host
asynchronous according to NVIDIA's API synchronization contract; replacing it
with cudaMemcpyAsync alone is not evidence of fewer synchronization points.
See [NVIDIA's synchronization contract](https://docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html).
PCG dot products now return device scalars; rho and pAp are read together once
per iteration, while next-rho and beta stay on device. This removes the separate
host wait for next-rho. Host scalar convergence and breakdown checks and the final FP64 true-residual
gate remain in place. No CPU solve fallback or CUDA Graph is introduced.
Coarse-factor validation now bounds-checks each row inside the checking kernel,
allowing all validation flags to be checked together before execution. Its
temporary permutation counters live through that check, eliminating intermediate
host checks and frees. Malformed coarse row pointers are rejected on the actual
gateway without permitting an out-of-bounds factor traversal.

Before uploading, the benchmark estimates numeric inputs, AMS/PCG workspace,
duplicate uploads and diagnostic arrays, plus library headroom. It rejects
cases that would consume the available GPU memory minus the larger of 2 GiB
or 20% of total device memory. The estimate is conservative admission control,
not a measured allocator peak. Develop within GPU memory; increase case sizes
gradually. A 10 GB matrix alone is not a 10 GB solver working set.
All Python-to-MATLAB AMS inputs now use versioned plain HDF5, including
artifacts and individual arrays exceeding 2 GB. The package-owned writer is
`radia_gpu.ams_exchange.write_exchange`; the matching MATLAB entry point is
`radia.gpu.readAMSExchange`. This is plain HDF5, not a MAT v7.3 struct file:
do not pass it to MATLAB `load`. FP64 values, zero-based int32 indices, array
shapes and nested auxiliary levels are preserved explicitly in the schema.
Streaming file hashes avoid an additional whole-file memory allocation.
Install the optional Python writer dependencies with `radia-gpu[exchange]`.

CUDA AMS sources, the pybind snapshot adapter and the build helper are owned
by `packages/radia-gpu/{native,tools}`. CPU hierarchy algorithms and NGSolve
FE semantics retain their CPU owner. `Build.ps1 -AmsCudaOnly` is a thin
monorepo dispatcher. The GPU distribution stages the standalone AMS MEX,
`radia.gpu.amsPCG`, the HDF5 reader and exact-build provenance; it does not
ship the validation-only CPU pybind adapter.

`python validation_test/ams_cuda/validate_exchange.py --large --out <owned-dir>`
checks Python-to-MATLAB layout and an actual array/file exceeding 2 GB.
This is exchange acceptance, not GPU capacity or a 10 GB solve certification.

[Saved numerical evidence](rtx5000_20261011.json) covers six small real
curl-curl plus positive-mass cases, with `beta_zero=False`, 96 to 5,326 DOFs,
both AMS cycles and one/two auxiliary levels. Pure curl-curl, complex
systems, A100 scaling and production model integration remain unverified.

At 5,326 DOFs, on the recorded four-thread CPU comparison:

| AMS cycle | CPU warm median | GPU warm median | GPU with transfers |
| --- | ---: | ---: | ---: |
| 1 | 15.80 ms | 41.14 ms | 46.55 ms |
| 7 | 17.22 ms | 34.23 ms | 49.25 ms |

All six true residual and preconditioner-parity checks pass. The prototype
does **not** beat CPU AMS at these sizes. Parallelizing coarse row products
within a warp reduced the earlier ordered single-thread GPU cost, but
per-call validation/allocation, kernel launches, host scalar synchronization
and small triangular dependencies still need measurement and optimization
before claiming fast GPU AMS. A100 must be tested directly on the correct
coordinated host.

## Residual and validation overhead experiment

The next implementation combines CSR products with residual output, combines
the smoother's residual subtraction and update, and starts a zero-initialized
Jacobi sweep directly with `b/l1`. Smoother updates still read a completed
matrix product: updating `x` inside the CSR kernel would introduce row races.
Independent input checks share a device error flag. Coarse row-pointer bounds
are checked before factor validation dereferences them, and all pending
checks complete before any solve. The hierarchy, factors, FP64 precision,
iteration rules and acceptance thresholds are unchanged.

`compare.py` compares two hash-verified MEX builds on the **same saved MAT
inputs**, with before/after/after/before order repeated three times. Each
block discards its first call and measures five warmed solves, giving 30
samples per build and case. It also records transfer-inclusive samples,
CPU-preconditioner probe errors, true residuals and gateway checks in JSON.
Run it only when competing CI and numerical work have finished:

```powershell
python validation_test/ams_cuda/compare.py --before $beforeBuild --after $afterBuild --inputs $inputCycle1 $inputCycle7 --out $comparisonOutput
```

Run `run.py` separately for independent direct-solution checks and CPU timings.
MEX load/clear between blocks and MATLAB's cached GPU allocator can affect
individual samples; retain all samples and compare medians rather than a
single fastest run. These are small acceptance measurements, not GPU scaling
evidence.

[Saved overhead experiment](rtx5000_fusion_20261011.json) records all samples
and six numerical-acceptance cases. The paired measurements use identical
5,326-DOF inputs with 30 warmed samples per build and cycle:

| AMS cycle | Before warm median | After warm median | Reduction | Before/after with transfers |
| --- | ---: | ---: | ---: | ---: |
| 1 | 37.56 ms | 31.17 ms | 17.0% | 43.84 / 37.70 ms |
| 7 | 33.45 ms | 27.49 ms | 17.8% | 39.73 / 34.11 ms |

This is a **shared-host relative comparison**: unrelated CPU CI environment
preparation/build was active, and the measurement process used BelowNormal
priority. The JSON records the observed CI steps. It is not an isolated-host
benchmark or evidence of beating CPU AMS; do not combine its samples with
the earlier CPU table to claim a CPU speedup. Preliminary switching-check
samples are excluded. A100 was neither started nor used.

All six numerical cases pass: maximum preconditioner probe relative error
is 3.45e-14, maximum FP64 true relative residual is 9.36e-10, and maximum
relative solution difference from independent CPU sparsecholesky is 4.76e-9.
All eleven gateway checks pass. The numerical-only run's timing samples
are excluded because it overlapped CPU CI as well.

## Idle-host repeat, 2026-10-11

[Saved idle-host evidence](rtx5000_idle_20261011.json) repeats the same six
small cases after repository CI finished. No CI was active or queued at the
start, midpoint or end. Two process/load probes five seconds apart preceded
measurement; total CPU utilization was 5% before and 9% afterward, with 2%
GPU utilization at both boundaries. Existing idle desktop/MATLAB sessions
were retained, so this is operational idle rather than exclusive reservation.
A100 was not started or used.

CPU uses Python/pybind11 + NGSolve with four threads; GPU uses the unchanged
FP64 mexcuda binary. These are the five-sample warmed solve medians from the
same acceptance run, excluding hierarchy setup for both routes. GPU gateway
validation/allocation is included; the last column also includes transfers.

| AMS cycle, 5,326 DOFs | CPU | GPU resident input | GPU with transfers |
| --- | ---: | ---: | ---: |
| 1 | 20.44 ms | 32.21 ms | 42.41 ms |
| 7 | 17.76 ms | 27.40 ms | 33.36 ms |

The GPU is still 1.58x / 1.54x slower than CPU for these cases. All six CPU
and GPU true residuals are below 1e-9; maximum CPU-preconditioner probe error
is 3.29e-14 and maximum solution difference from independent sparsecholesky
is 4.76e-9. All eleven gateway checks pass. Stopping-norm and iteration-count
semantics remain as described above.

A separate identical-input ABBA comparison retains 30 warmed samples per
build/cycle: cycle 1 improves from 39.01 to 32.52 ms (16.6%), and cycle 7
from 34.77 to 28.60 ms (17.7%). These medians belong to that separate
comparison, not the five-sample CPU/GPU table.

Zero-RHS calls skip PCG iterations but retain gateway validation/allocation,
cuBLAS setup and the RHS norm. Their paired medians drop from about 7.6 ms
to 3.2 ms. Apply-only calls, including one AMS cycle and their own setup,
take about 3.5--3.6 ms after the change. These wall-time diagnostics are not
exact kernel/allocation counters and use different workspace paths. Do not
subtract zero-RHS time from apply-only time. Gateway reuse alone therefore
cannot explain or eliminate the full CPU/GPU gap; iterative launch/sync and
coarse-solve costs need further profiling before selecting an optimization.

Source of the sparse factor layout: NGSolve's public `SparseCholeskyTM`
accessors in `sparsecholesky.hpp` for 6.2.2607. The GPU kernel evaluates the
standard ordered unit-triangular/diagonal solves from those factors; it
does not recompute or replace NGSolve's factorization.

### Package/HDF5 and compute GPU acceptance (2026-10-11)

AMS now belongs to `packages/radia-gpu`, including its mexcuda builder,
MATLAB wrapper and checked HDF5 reader. All snapshots use plain HDF5
(`radia.ams-hdf5.v1`). An actual 2,147,490,184-byte exchange passes
Python/MATLAB shape, dtype and full-value checks. This validates exchange,
rather than a large GPU solve.

[Compute evidence](compute_gpu_20261011.json) records eight CPU/GPU cases,
versions, build/source hashes, idle observations and raw timing samples.
The already-running user-provisioned A100-SXM4-40GB compute host was used;
the agent did not start a machine or access the old A100 SSH alias.
CPU uses Python/pybind11 + NGSolve; GPU uses mexcuda. Maximum CPU AMS
solution difference is 2.57e-9 and maximum GPU true residual is 9.36e-10.
GPU memory admission checks pass for all eight cases.

This separate same-input ABBA comparison retains 30 warmed samples per
variant/case. Its baseline is the previous optimized binary, rather than
the earlier baseline in the historical RTX measurements above.

| DOFs | Cycle | Before | After | Time reduction |
| ---: | ---: | ---: | ---: | ---: |
| 5,326 | 1 | 25.67 ms | 24.72 ms | 3.68% |
| 5,326 | 7 | 22.28 ms | 21.45 ms | 3.76% |
| 12,217 | 1 | 24.02 ms | 22.76 ms | 5.25% |
| 12,217 | 7 | 21.17 ms | 19.98 ms | 5.64% |

In the separate CPU/GPU acceptance run at 12,217 DOFs, CPU takes
17.09/17.11 ms and GPU takes 22.68/19.92 ms for cycles 1/7. CPU remains
faster at these sizes. Both routes satisfy independent true-residual
gates, but use differing stopping norms. These results do not establish
performance or allocator peak memory for a 10 GB problem.
Package CI passed at code commit `170717fa7`, including installed-wheel
AMS acceptance. Exchange evidence remains in
[package_hdf5_20261011.json](package_hdf5_20261011.json).

## AMS setup acceleration

`setup_benchmark.py` compares the same hierarchy and solver before/after
explicit OpenMP setup (`setup_threads=1` versus the requested count), using
ABBA order. It requires identical hierarchy arrays, checks CPU probes and
true residuals, and verifies GPU results against the CPU reference. Setup
phase counters identify strength, coarsening, interpolation, transpose,
Galerkin, l1/allocation and coarse-factor costs. Coarsening and coarse
factorization retain their serial order.

The transfer-inclusive total is hierarchy setup + snapshot export + HDF5
write/read + GPU upload/solve/download. Validation-only CPU probes and MATLAB
startup are reported separately. Resident solve timings show the repeated
application cost when a caller keeps the hierarchy on the GPU. This is a
reusable FEM/AMS kernel study for future control applications, not a CPU/GPU
winner comparison or certification of a particular Simulink model.

Run on 100 for numerical acceptance and diagnostics:

```powershell
python validation_test/ams_cuda/setup_benchmark.py --out "$env:TEMP/ams-setup-run" --maxh 0.08 --setup-threads 4 --rounds 2
```

Formal performance acceptance requires an idle GPU-capable mdx. A CPU-only
run uses `--cpu-only` and leaves GPU/transfer timing pending. Never start A100
to satisfy this benchmark. Each run writes `result.json`, including failures.

## MATLAB-owned native setup

The experimental candidate also exposes `radia.gpu.prepareAMS(model)` within
MATLAB/mexcuda. It reuses the retired MEX gateway's shared `HypreBasedAMS`
constructor and NGSolve sparsecholesky; no Python adapter or prebuilt Python
solver hierarchy is required. Supply a versioned numeric model containing
`A`, `G` (sorted int32 CSR), `coordinates` (vertices by three), an int32 free-DOF
mask, and cycle 1 or 7. NGSolve still owns FE orientation and assembly semantics.
The default uploads the resulting hierarchy for `radia.gpu.amsPCG`; `Upload=false`
returns host data for explicit transfer management. TaskManager regions are
compiled out; indirect upstream dependencies are permitted by the clarified
policy. Python DLL dependencies do not imply interpreter initialization.

Run `native_setup_validation.py --out <owned-output-directory> --maxh 0.12 0.08`
for independent CPU accuracy and serial/OpenMP setup comparison. HDF5 carries
validation inputs; GPU users can supply their numeric model directly. JSON
includes cold/warm, setup plus upload/solve/download, memory, versions, input
rejections and source/build identity. Build-host timings remain diagnostic.
## Nonlinear GPU-resident correctness study

`nonlinear_resident.py --out <job-owned-output> --maxh .35`
runs an independently authored SI unit-cube magnetostatic model with
`H=(1+4*|B|^2)B`, positive unit gauge mass, homogeneous tangential Dirichlet
conditions, and a reproducible manufactured discrete excitation. The explicit
Netgen alternative generates a `.vol`, which is reloaded by NGSolve with a
six-boundary contract. This is a nonlinear constitutive/solver test, not a
measured steel curve, open-boundary device, or eddy-current model.

NGSolve owns mapped HCurl shape functions and the independent CPU nonlinear
assembly. MATLAB evaluates material state, energy, residual and consistent SPD
tangent on GPU. The frozen AMS hierarchy and nonlinear state remain GPU arrays;
only scalar convergence controls cross the boundary during reuse. Rebuilding
AMS explicitly downloads tangent CSR values for the existing native host setup.
The Python-generated FE maps are independent validation inputs; the production
AMS API continues to accept numeric operators without a Python adapter.

Three routes run in order 1,2,3,3,2,1: rebuild each Newton iteration, freeze AMS,
and adaptive reuse. Adaptive reuse rebuilds on the next iteration if PCG exceeds
`max(30, 2*iterations_at_last_rebuild)`. There is no backend fallback. Each linear
solve checks a true residual at 1e-10, each nonlinear run at 1e-9, and independent
CPU solution, B-field RMS and nonlinear residual parity must be <=1e-7.

The current implementation assembles sparse CSR operators from six mapped curl
vectors and 36 CSR slots per affine tetrahedron. It supports lowest-order HCurl
on affine tetrahedra; other element orders or geometry fail explicitly. NGSolve
owns orientation and mapped shape functions. No global dense matrix is formed.
`radia.gpu.AMSContext` owns a checked native handle, persistent PCG buffers and
cached CUDA graph. Hierarchy replacement preserves its PCG workspace. Nonlinear
residual, tangent and element-energy buffers also remain on the device. Native
buffer allocation counters must remain zero during each resident linear solve;
returned MATLAB arrays and CUDA graph internals are outside this counter.

CPU reference assembly and solves use caller-owned NGSolve TaskManager regions;
MATLAB native setup uses OpenMP with serial sparsecholesky coarse factorization.
Both backends enforce the same true linear residual. Near floating-point energy
roundoff, a line-search step must reduce the actual residual by ten percent or
pass the final nonlinear gate; tolerances are not relaxed.

Run `nonlinear_resident.py --out <job-owned-output> --maxh .04 --threads 8
--setup-threads 8 --contract` for sparse scaling and lifecycle checks. Total times
include nonlinear updates, line search, AMS setup, PCG and final download; GPU
total also includes upload. Shared FE preparation and MATLAB startup are recorded
separately. Build-host numbers remain diagnostic timing, not formal performance
acceptance. The historical dense study below is superseded for scaling.

The accepted [nonlinear evidence](results/nonlinear-resident-20261013.json)
contains 260- and 698-DOF cases on the same native binary, plus cycles 1 and 7
linear setup parity. The larger case reproduced a nonfinite gradient coarse
factor before the fix. Its known constant nullspace now receives an explicit
potential gauge per connected zero-row-sum component before sparsecholesky;
no diagonal perturbation or numerical fallback is used.


## Sparse nonlinear scaling evidence

[Accepted sparse study](results/nonlinear-sparse-resident-20261013.json) includes
260, 5,326, 86,735 and 251,896 DOFs on the same native binary, with independent
Python CPU assembly/AMS and MATLAB GPU assembly/AMS. The largest case has
205,728 affine tetrahedra and 4,000,702 CSR nonzeros. Both CPU assembly/solve and
GPU host setup use eight threads; GPU computation uses a Quadro RTX 5000.

Largest-case total times (median of two samples per route, order 1,2,3,3,2,1):

| AMS policy | Hierarchy builds | Python CPU | MATLAB GPU | CPU/GPU |
| --- | ---: | ---: | ---: | ---: |
| Every Newton update | 8 | 14.205 s | 8.992 s | 1.58x |
| Frozen | 1 | 21.212 s | 6.505 s | 3.26x |
| Adaptive reuse | 2 | 14.088 s | 4.642 s | 3.03x |

Adaptive GPU reuse reduces total time by 1.94x versus rebuilding every update;
freezing is slower than adaptive reuse despite avoiding another setup. Each run
requires eight Newton updates. Maximum independent solution difference is
1.52e-13, B-field relative RMS difference 2.46e-15, and independent nonlinear
residual 2.78e-14. Every linear solve passes the 1e-10 true-residual gate.

The retained PCG workspace is 10,075,864 bytes. No additional native Buffer
allocation occurs during resident solves. Maximum sampled device allocation
delta is about 348 MiB; this is not a whole-process peak measurement. Lifecycle
checks pass after hierarchy replacement and caller-array release, reject stale
and incorrectly typed handles, and verify zero-RHS reset and context cleanup.
Shared FE preparation (48.47 s) and MATLAB startup (11.09 s) are separately
reported rather than hidden in solver timings. These are manufactured-model
build-host diagnostics with two samples per policy, not physical-device or
formal timing acceptance. A100 was not used.
