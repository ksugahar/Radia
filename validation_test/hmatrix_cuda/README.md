## 2026-10-18 bounded scale research candidate

The native MATLAB route has explicit experimental options:

```matlab
solution = radia.gpu.solveHDivElements(checkedMeshPath, 2, [1 .2 .1], ...
    2.5, 1e-11, 1000, 'host-cache', 'mass-cholesky', 'sparse');
```

`mass-cholesky` builds a dense FP64 mass factor on the GPU with MATLAB-bundled
cuSOLVER and applies its two triangular solves with cuBLAS. The original material
operator and final true-residual check are unchanged. Non-SPD factorization fails
without switching to diagonal or CPU execution. The default remains the existing
mass-diagonal preconditioner. Factorization has quadratic storage; this is not a
sparse Cholesky implementation or an AMS/IC acceptance result.

`sparse` stages both the FE mass and oriented charge map as MATLAB sparse arrays;
`sparse-mass` stages only the mass sparsely. The default remains `dense`. The same
per-element accumulation, monomial transformation and curved zero floor are used.
`checkHDivSparseMass.m` checks exact storage and whole-solve equivalence;
`checkHDivMassPreconditioner.m` checks repeated solves, zero RHS and non-SPD failure
isolation. Checked-handle reload tests remain required.

The host-potential cache uses bounded source tiles. Mapped quadrature positions and
charge weights are reused, while affine analytical moments, touching curved Duffy
integration, non-touching tensor quadrature and ordered symmetric averaging remain
unchanged. These optimizations do not turn the dense high-order Gram path into
adaptive HACApK assembly.

The configured bounds are 16,384 DoF, 4,096 tetrahedra and 8,192 charges, with
preallocation guards and the existing device-memory reserve. A configured bound is
not a validated geometry or capacity claim. `prepare_hdiv_scale.py` generates an
independent public unit sphere and records dense and tiled-cache memory lower
bounds; actual cuSOLVER workspace/runtime and sampled host-memory limits are
recorded separately. `validate_scaled_hdiv.py` certifies the complete solution by
independent CPU operator actions, avoiding a dense DoF-by-DoF oracle. Optional
`--abba-cycles 2` rebuilds both workflows for every timed sample with fixed CPU
threads and records transfer-inclusive GPU API totals, controller time, residuals,
iterations and memory. First benchmark iterations follow parity warm-ups and are
not process cold starts. The in-process CPU memory sampler has unquantified timing
overhead, and the retained CPU certification operator is included in its baseline.

These are bounded linear, uniform-isotropic public-sphere research checks. They do
not certify nonlinear BH, topology optimization, general interfaces/constraints,
images, arbitrary curved geometry, or a released same-wheel distribution. Raw and
failed-run evidence is retained privately; release acceptance remains separate.

### Scaling beyond the bounded dense route

Sparse FE staging does not remove the dense charge Gram or the optional dense
mass-Cholesky factor. At 100,000 DoF a single FP64 DoF-square factor requires
80,000,000,000 bytes, before workspace or the Gram operator. The scale preflight
reports separate `mass_cholesky_array_gate` and `dense_gram_array_gate` values;
passing the tiled potential-cache gate does not certify solver capacity. These
are lower-bound checks, not complete allocation admission or numerical acceptance.

The next architectural work is an explicit sparse fixed SPD mass preconditioner
and high-order/curved entry evaluation connected to existing HACApK adaptive
leaves. Both remain unimplemented and unvalidated in this research candidate.
Preserve oriented FE charge maps, curved reference frames, touching Duffy rules,
normalization, symmetric leaves and the original-equation true residual. Do not
increase dimension caps, change integration, or replace curved elements with
point sources to obtain a larger benchmark. Compressed storage and convergence
must be measured on actual meshes before claiming 100k capability or improved
speedup at that size.

The staged `curved_entry_batch` diagnostic evaluates requested one-based charge
indices in bounded 256-entry GPU batches with resident geometry and quadrature.
It reuses the direct curved integral and upper-triangle symmetric averaging;
it does not allocate a dense Gram. `checkHDivIndexedEntries.m` prepares a bounded
full-operator, invalid-index, empty-request and batch-boundary comparison. This
new boundary passed a clean native build and a 160-charge curved full-operator
comparison (relative Frobenius difference 1.1016030385600124e-16), empty requests,
invalid indices, the 256-entry batch boundary and handle cleanup. The initial
HDF5 save failed missing required labels; a corrected validation script passed
with the same binary, and both records were retained. This is GPU-to-GPU operator
parity only, not independent CPU or full-solver acceptance. It is not yet
connected to HACApK callbacks, clustering, normalization or compressed leaves.
It does not implement the sparse mass preconditioner or enable 100k solves.

The bounded research candidate adds `create_hdiv_hacapk(config,
epsilon, leafSize, eta, maxRank)` with explicit parameters and a 4,096-charge
bound. It uses the same existing C ACA dispatch, diagonal normalization, rank-cap
dense promotion, native leaf packing and resident sparse `B^T G B` application.
P2 support bounds use the convex hull of quadratic Bernstein controls, including
edge controls `2*midpoint - (vertex_i+vertex_j)/2`, with roundoff padding; vertices
alone do not bound curved supports. The entry oracle is direct and rejects the
cache/profile flags rather than silently treating them as supported. No default
assembly route or leaf size changes. A clean full native build and
`checkHDivHacapk.m` passed affine and P2 curved operator/material GPU route
comparisons. The curved original-direct equation residual was 1.30114e-12;
its ten leaves were all dense. `checkHDivHacapkSeparated.m` exercised one
low-rank-format leaf on two translated affine tetrahedra, with full-operator
relative Frobenius difference 5.35540e-18. Actual ranks and compression ratios
have not yet been recorded, so this is not evidence of a storage or speed gain.
Independent CPU checks, high-order curved low-rank accuracy, conservative-bound
regression and leaf-size/performance sweeps remain required. Large leaves are a
hypothesis to measure, not a universal MMM rule. See
`results/curved-hacapk-research-20261018.json`; no release acceptance is claimed.

The bounded wrapper/statistics research candidate accepts an
explicit `hacapk` mode as follows; no HACApK parameters are silently defaulted:

```matlab
options=struct('epsilon',1e-10,'leaf_size',64,'eta',2,'max_rank',64);
solution=radia.gpu.solveHDivElements(checkedMeshPath,2,[1 .2 .1], ...
    2.5,1e-11,1000,'hacapk','mass-cholesky','sparse',false,options);
```

`validate_scaled_hdiv.py --assembly-mode hacapk` requires the corresponding
four `--hacapk-*` options and still certifies the result with the independent
CPU original operator. Statistics distinguish actual rank, stored FP64
coefficient count, dense storage for the same upper leaves, and low-rank leaves
whose factors actually save coefficients. Metadata, vectors, sparse maps,
workspace, host and runtime memory are separate from coefficient storage.

The 288-DoF P2 sphere passed the independent CPU comparison: coefficient error
6.48756e-14 and CPU original-operator residual 1.30139e-12. The source/wheel
adoption chain for that CPU reference remains provisional. Its HACApK assembly
has ten dense leaves and the unoptimized direct entry oracle takes about 57 s;
this route is not a speed improvement over the accepted cached dense candidate.
Only one first and one warm diagnostic were run, not ABBA or process cold start.

The separate affine leaf-size-8 case reveals why format counts are insufficient:
its one low-rank-format leaf is rank 10 on a 10-by-10 block, using 200 factor
coefficients instead of 100 dense coefficients. Total stored coefficients are
1,828 versus 1,728 for the same dense upper leaves, an increase of about 5.8%.
No compression benefit is accepted. Next optimize repeated mapped source-moment
evaluation in the indexed oracle before attempting expensive scale sweeps.
See `results/hacapk-cpu-research-20261018.json` for numerical and storage records.

The explicit research `hacapk-cached` candidate reuses mapped host
positions, monomial weights and the accepted fused source-moment kernels. It
keeps at most two source-potential tiles with a default 256 MiB cache bound; it
never allocates a dense Gram in the oracle. Ordered entries are evaluated before
the original symmetric averaging. Defaults, quadrature, curved reference frames,
touching Duffy rules and ACA settings are unchanged. It is an explicit alternative
to `hacapk`, not an automatic fallback. `checkHDivCachedEntries.m` passed
complete-operator, forced eviction, partial
moment-group boundary, invalid-budget and liveness checks. Its diagnostic
`curved_entry_batch_cached` command optionally takes a smaller cache byte budget
and returns cache counters. The explicit whole-solver budget option is described
below; omitting it retains 256 MiB.

The compiled cache candidate passed affine and curved full-operator comparisons
(relative errors 1.18e-16 and 1.10e-16), including 400 forced tile loads for the
curved case. The independent CPU 288-DoF whole solve passed with coefficient
error 5.86e-14 and original-operator true residual 1.30e-12. Fresh two-cycle ABBA
gave CPU-16/GPU API medians 0.434248/0.0496885 s (8.7394x), including setup,
assembly, transfers, solve and recovery. Engine/controller transport is excluded;
the controller-observed GPU median was 0.429313 s. These timing boundaries must
not be confused. Exact accounted GPU owned peak was 2,185,476 bytes, excluding
runtime/global peaks. This is a small bounded research case with ten dense
leaves, not compression, large-scale, process-cold or release acceptance. The
CPU immutable-wheel adoption chain remains provisional. The old uncached 57 s
record remains a different-session diagnostic, not a paired old/new ratio.
See `results/hacapk-indexed-cache-research-20261018.json`.

At 2,412 DoF / 1,176 charges, the same cache binary passed independent CPU
comparison (coefficient error 1.32e-13, true residual 7.93e-12). Fresh two-cycle
ABBA measured CPU-16/GPU API medians 4.90803/0.618260 s (7.93846x), with
117,017,988 bytes exact accounted GPU owned peak. Leaf size 64 produced 251
dense leaves. A separate leaf-size-256 probe also passed numerically but produced
17 dense leaves and 812,624 stored coefficients versus 719,148; its single warm
time did not improve. This is not evidence for universally larger leaves.
The next 6,594-DoF input requires six source tiles under the default two-slot
cache. The compiled research extension permits explicit `cache_bytes` in
`hacapk-cached` options, at most 1 GiB, retaining the default 256 MiB budget,
the 2 GiB total oracle workspace check and device reserve. The validation CLI
spells this `--hacapk-cache-mib`. A requested budget is a bound, not an automatic
allocation or a capacity/accuracy guarantee. Recomputation and scaling remain
untested for tile-evicting whole solves. See
`results/hacapk-cache-scale-2412-research-20261018.json`.

The default-budget 6,594-DoF run reached its 240 s wall limit without a numerical
result. The explicit 1 GiB budget then passed the full operator/eviction gates,
independent CPU regression at 288 DoF, and whole-solve comparison at 6,594 DoF
(2,724 charges). Fresh two-cycle ABBA measured CPU-16/GPU API medians
12.9006/2.31717 s (5.56741x), coefficient error 6.16e-12 and independent original
residual 6.94e-12. Exact accounted GPU owned peak was 710,470,180 bytes, including
the source-potential cache. This is not runtime/global peak memory.
Only two of 1,484 leaves actually save storage: coefficient reduction versus the
same dense leaves is about 0.045%, so substantial H-matrix compression remains
unproven. The old timeout is neither numerical disagreement nor a paired speedup
measurement. Defaults, integrals, reserve and 2 GiB oracle workspace checks are
unchanged. The 4,096-charge full-HACApK cap still rejects the existing 10,344-DoF
input with 4,208 charges; no capacity or release claim follows. See
`results/hacapk-cache-budget-6594-research-20261018.json`.

The compiled capacity extension bounds curved
full assembly at 4,608 charges. Above 4,096 it additionally requires leaf size
at least 64 and rank cap at most 256. Before the existing C allocations it checks
signed-index products and a conservative 384 MiB factor-coefficient bound.
An explicit cache budget may be at most 2 GiB, but the total oracle-owned
workspace must still fit 2 GiB including geometry and metadata, and device
reserve checks still apply. Defaults and the analytic 4,096-charge cap remain.
The public 10,344-DoF input needs about 1.71 GB potential storage and an 856 MB
mass factor. Its initial 2 GiB request exposed the signed-int scalar parser;
byte budgets now have a separate bounded parser, leaving C index limits intact.
The corrected binary passed the full small-operator/cache-eviction checks,
extended-capacity rejection checks, and independent CPU material solve.
Fresh two-cycle ABBA measured CPU-16/GPU API medians 21.2148/4.97606 s (4.26337x),
coefficient error 5.23e-11 and original-operator true residual 6.27e-12.
Accounted GPU owned peak was 1,712,869,764 bytes. The 32 storage-saving low-rank
leaves reduce coefficients by only 0.753% versus the same dense leaves.
The previously accepted dense4Mi path is still faster on this public sphere;
these are separate-session/artifact results, not a paired path benchmark.
HACApK assembly/transfer dominates its whole time, so this does not establish a
general hierarchy performance advantage, 100k capacity, or release readiness.
See `results/hacapk-capacity-10344-research-20261018.json`.

The compiled entry-batch research candidate adds
`build_hmatrix_cuda.py --entry-batch-size 1024` (default remains 256). The same
bounded macro controls C active row/column extraction, row-major dense fill,
the C ABI error boundary and native oracle capacity. It does not change ACA
pivot/update order, entry contractions or CUDA thread-block sizes. Scalar CPU
builds retain their existing route. C scalar/batched mask/permutation/boundary
contracts passed 60 cases each at both configured capacities; shared ACA+
compression regression passed its 12 cases at 1024. Full native C scalar/batch
synthetic assembly also passed at the unchanged default capacity.
The compiled candidate passed full small-operator and independent CPU 288/10k
whole-solve checks. At 10,344 DoF, fresh two-cycle ABBA measured CPU-16/GPU API
medians 21.1954/3.61354 s (5.86556x), coefficient error 5.23e-11 and original
true residual 6.27e-12. Entry batches fell from 39,734 to 14,364 while entry
count, leaf ranks and coefficient storage stayed unchanged. The previous
256-entry artifact/session is a diagnostic control, not a paired path timing.
Owned setup peak rose by 36,864 bytes to 1,712,906,628 bytes; no memory-reduction
claim follows. The separate dense4Mi path remains faster on this public sphere.
See `results/hacapk-batch1024-10344-research-20261018.json`. This is explicit
research configuration, not same-wheel distribution or release acceptance.

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

100 timing is diagnostic, not formal idle-historical_cpu_worker performance acceptance. Exact
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
process startup or historical_cpu_worker performance. A later timing attempt was rejected by the
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

## Exact rank-budget promotion bridge

`analytic_aca_leaf` retains the existing ACA+ diagnostic inputs and CPU
rank-budget promotion condition. A saturated rank below both dimensions
is discarded and the physical dense block is evaluated on CUDA in bounded
256-entry batches. The MATLAB result uses column-major dense storage; native
resident dense leaves still use row-major storage. Non-saturated results retain
the factor response. This is a leaf bridge, not adaptive whole-matrix assembly.
Generated dense payloads above 2 GiB are rejected before allocation; HDF5
remains pending. The larger check compares promotion against direct CUDA entries at 1e-12
and against independent CPU physical entries at the existing 1e-8 block
limit. A first attempt requiring 1e-12 CPU physical parity failed at
1.726e-9; the direct-entry comparison isolates promotion from integration
roundoff. This does not certify CPU physical parity at 1e-12. `results/dense-promotion-20261012.json` records direct CUDA difference zero
and CPU physical difference 1.726e-9. Full adaptive assembly remains pending.

## MATLAB host scheduling boundary

`hmatrix_serial_host.c` implements the existing HACApK C scheduling API on
the owning MATLAB thread. It replaces, rather than links alongside, the
TaskManager implementation in a future full native assembly target. CUDA
entry bindings are thread-local, so serial control preserves ownership without
requiring worker binding or oversubscribing CI. `check_entry_batch.py
--serial-host` checks empty work, 257 ordered iterations, one job and thread
identity. `results/serial-host-20261012.json` records that C contract only.
The current narrow CUDA target still does not build the full cluster/leaf
implementation; PCA/LAPACK integration and full physical acceptance remain
pending. No speedup is inferred from serial host control.

## Same-artifact A100 subgate

The user authorized tests on the already running NGSolve A100 host.
`validate_same_artifact.py` verifies every staged input/module/wrapper hash,
uses an owned MATLAB Engine session and releases it on success or failure.
Python transports CPU validation inputs only; it is not a GPU user adapter.
The 100-built `15e2bef2a` binary passed the same numerical/lifecycle suite on
A100-SXM4-40GB, including exact dense promotion and CPU physical parity.
The optional profile uses Windows system CPU/memory samples and GPU utilization
before/after; these samples do not exclude transient external load or measure
process/global peaks. The small 275-dimensional, three-apply ABBA diagnostic
records 4.292 ms rebuild versus 2.481 ms resident (1.73x), with 573312 owned
device bytes in both modes. See `results/a100-same-artifact-20261012.json`.
This is not full HACApK assembly, an HDiv-MMM solve, true-residual acceptance
or full-solver performance. The initial profile attempt failed on absent
psutil and is retained separately; no host package was installed.

## MATLAB-native RT0 material path

`radia.gpu.solveHDivRT0` loads a checked affine-tetra mesh through NGSolve,
constructs RT0 mapped divergence/normal charges, mass and applied-field RHS
inside the mexcuda host bridge, then constructs the constant-charge Gram on
CUDA. No Python FE snapshot is required by this MATLAB entry point. The
current bounded path rejects curved/non-tetra meshes and more than 16 cells;
it uses exact dense Gram, not adaptive whole HACApK assembly. The projection
module now contains sm75/sm80/PTX and links supported nglib mesh access.

The FP64 material PCG keeps mass, vectors and the demagnetization operator on
GPU. It uses an explicitly selected inv_chi mass-diagonal preconditioner;
CPU reference uses NGSolve mass-Riesz sparsecholesky. cuBLAS scalar reductions
synchronize host control; no synchronization-free solver claim is made.
Convergence is checked by recomputing the original equation residual, with
restart on recursive/true-residual disagreement and loud failure at limits.

One and eight tetra physical subgates passed on 100 and the same hashes on
the authorized A100. The eight-tetra A100 coefficient/mean-magnetization
relative differences are 1.126e-11/3.663e-12; independent true residual is
9.729e-12. Accounted owned peak is 31852 bytes. The recorded first solve
0.291 s follows a separate FE parity check and is not warmed performance.
See `results/hdiv-native-rt0-a100-20261012.json`. Experimental MATLAB entry
points remain excluded from the wheel until full release acceptance.

The first live-refined CPU reference exposed four zero mass diagonals in
28 DOFs (24 free). Saving/reloading the same mesh produces 24 DOFs with
positive minimum mass eigenvalue 0.0287241. Validation now reloads its mesh
artifact, matching the native mesh-file entry point; no alternate direct
solver or weakened SPD check was introduced. The initial A100 load failed
on the isolated dependency bin path and was corrected without changing
global installations. Failure evidence is retained. General HDiv-MMM and
adaptive HACApK validation remain incomplete.

## Full native host assembly dependency gate

`check_full_host.py --mkl-root <isolated-venv>/Library --out <new-job>`
builds the existing complete C cluster/leaf/fill implementation, including
MKL PCA/LAPACK dependencies, with `hmatrix_serial_host.c` instead of
TaskManager scheduling. The original production ACA dispatch (`param[60]=1`)
is preserved. Both BBOX and PCA strategies on an independently authored
132-point smooth synthetic kernel pass scalar and batch dispatch against
direct matrix multiplication (maximum relative difference 3.353e-13).
`results/full-host-assembly-20261013.json` records source and executable hashes.
This closes the host dependency/assembly subgate; the batch callback in this
test evaluates the synthetic CPU kernel, not CUDA. The following section
records the subsequent CUDA connection. No timing was measured in this gate.

## Complete constant-charge CUDA hierarchy path

`Build.ps1 -HMatrixCudaOnly -HMatrixMklRoot <isolated-venv>/Library
-RequireNativeProvenance` links the complete existing C cluster/leaf/ACA
implementation with supported MKL and the serial MATLAB host scheduler.
The new `create_analytic_hacapk` command accepts constant affine tetra/triangle
charges, explicitly bounded to 4096 sources. It preserves BBOX support radii,
production ACA dispatch, diagonal normalization and rank-cap dense promotion.
CUDA computes both compressed-entry requests and exact dense leaves. Dense
leaves use contiguous batches of at most 256 entries; CPU default builds retain
the scalar order. Scoped callbacks capture C ABI failures and owners reject
failed builds before creating a resident handle. Symmetric fill, geometry-radius
and cluster-strategy state is reset on exit.

The 132-tetra direct-CPU parity subgate passes at 1.440e-11 relative error,
with 14272 evaluated entries in 113 batches, 7 stored upper leaves and 99688
owned resident device bytes (292324 accounted setup peak). At this tolerance,
all selected low-rank leaves exhausted the rank budget and were promoted to
dense: this is full hierarchy assembly, **not** successful compression.
`solveHDivRT0(...,'hacapk')` connects that builder to native NGSolve FE setup
and the resident material solver. The eight-tetra physical path retains CPU
coefficient difference 1.126e-11 and independent true residual 9.729e-12.

`measure_hdiv_assembly.py` compares the previous scalar-dense CUDA assembly
with batched CUDA assembly using the same physical mesh and CPU reference.
Its guarded ABBA timer includes FE setup, assembly, transfers, solve, result
recovery and destruction. It is an idle-100 diagnostic, not formal historical_cpu_worker
performance acceptance. First payloads are distinguished from cold processes;
available-memory samples are not global peaks. General high-order/curved FE,
images, constraints, nonlinear physics, large compressed-leaf validation and
release acceptance remain pending.

`results/full-cuda-assembly-20261013.json` records the latest checked build
(`984875630`, `job-fop0ay52`) and same-hash A100 numerical subgates. The
idle-100 eight-tetra ABBA whole-solve medians are 1.30363 s before batching
and 0.0480764 s after batching (27.12x). CPU coefficient difference remains
1.126e-11, independent true residual 9.729e-12, 21 iterations, and accounted
owned peak 28996 bytes in both modes. These are CUDA-before/CUDA-after
diagnostics on this bounded physical case, not CPU-vs-GPU or general-solver
performance acceptance. A100 validation here uses independent CPU FE fixtures;
native FE initialization on A100 was not reaccepted for this latest artifact.

## Compression tolerance diagnostic on 256 constant charges

`validate_resident.py --aca-block-size 128 --full-assembly --full-max-rank 64
--full-epsilon 1e-6` keeps the production defaults and the existing independent
CPU parity limit (1e-8) unchanged. One low-rank upper leaf survives; matvec
relative difference is 5.815e-10 and owned resident bytes decrease from 336392
to 217704 (35.3%). No timing improvement is established. The direct CPU
128-by-128 offdiagonal block has optimal rank-64 relative Frobenius tail
4.472e-10, above the stricter 1e-10 assembly tolerance. Thus dense promotion
in the strict case is not evidence of a CUDA-only failure. The origin of the
small singular-value floor still needs isolation. A single matvec comparison
is not an operator norm certificate. See `results/epsilon-sensitivity-20261013.json`.

`results/full-operator-20261013.json` strengthens this diagnostic by comparing
all 256 basis-vector GPU applications with the complete independent CPU direct
matrix: relative Frobenius error 9.960e-10 (unchanged limit 1e-8), symmetry
error 2.536e-17. This certifies the bounded constant-charge operator comparison,
not general HDiv-MMM physical or release acceptance. No timing is measured.

## Native RT0 physical scaling to 64 tetrahedra

The bounded native FE setup now admits up to 128 affine tetrahedra, 1024 DOFs
and 1024 charges, rejecting over-budget or unsupported geometry before the
corresponding dense allocation. A fresh `ProjectionCudaOnly` build from
`7bf728b70` passes the twice-refined public tetra fixture with the existing
full-HACApK artifact. Native mass/RHS/mean operators agree with independent
NGSolve assembly; material coefficients differ by 9.680e-12, mean magnetization
by 4.611e-13, and independent true residual is 8.674e-12 (32 iterations).
Accounted owned setup/solver peak is 166916 bytes. Initial payload whole time
is diagnostic only, not a warm comparison. See `results/hdiv-64tet-20261013.json`.
General high-order/curved/constraint/nonlinear and release gates remain pending;
the new projection binary has not yet been reaccepted on A100.

## High-order and curved FE staging boundary

`inspect_hdiv_elements` in the experimental projection mexcuda module stages
NGSolve-owned mapped HDiv mass and P2 reference-node samples. Orders 1..3 are
admitted for affine tetrahedra, orders 1..2 for curved tetrahedra, with explicit
128-cell/1024-DOF budgets. It starts no TaskManager and requires no Python
fixture for MATLAB use. `validate_element_setup.py` independently assembles
the NGSolve mass in a caller-owned Python TaskManager and compares the native
mapped nodes on an affine public tetrahedron and an original unit sphere.
This diagnostic command is not a charge-map constructor or a GPU material solve.

Remaining integration must preserve these existing mathematical contracts:

- `_charge_basis` uses local L2/SurfaceL2 mass inverses and an orientation-aware
  change of monomial basis for **each** element, not one reused element basis.
- `_charge_basis_curved` uses reference-frame monomials, mapped P2 nodes and
  the existing structural-zero floor. Its supported curved orders are 1 and 2.
- `PhiInnerHOHostVec` uses analytic moment contractions, including the distinct
  reference-moment branch for higher degree. RadFld scalar kernels alone do
  not replace this host-block contract.
- `QuadBlockHOTet` uses the existing singular host-vector Duffy path on touching
  curved pairs; only non-touching curved pairs use the existing direct tensor
  rule. Symmetrization, image transforms and near/far selection must survive.
- The resulting charge map and Gram must feed the resident material equation;
  compare coefficients and physical quantities with independent CPU assembly
  and check the true residual before any whole-workflow timing acceptance.

FE staging parity alone leaves general HDiv and release acceptance pending.
The NGSolve A100 is offline and is not started or contacted for this gate.

### Bounded curved whole solve (2026-10-13)

`results/curved-hdiv-solve-20261013.json` records native NGSolve charge-map setup,
existing curved singular integration, resident material PCG and physical recovery
for the original public 13-cell P2 sphere. HDiv orders 1 and 2 pass independent
CPU coefficient, mean magnetization and true-residual gates. The earlier volume
normalization failure is retained in durable evidence. This dense near-rule
subgate does not certify general adaptive/high-order HDiv, constraints, images,
nonlinear materials, performance, packaging or release. A100 was not used.

The latest same-artifact cases additionally cover the public affine tetrahedron
at HDiv orders 1 and 2 through `solveHDivElements`. The wrapper preserves the
existing CPU capability limit of order 2. Native charge-map staging at order 3
is not whole-solver acceptance. Affine analytic reference moments and curved
touching Duffy rules feed the same resident material solve. Latest curved order 2
coefficient disagreement is 5.02e-10 and independent true residual 7.95e-12.
`validate_hdiv_artifact.py` evaluates the same native setup and built artifacts
on an independently staged host; Python supplies independent validation fixtures
and Engine orchestration only, with no GPU-user Python adapter.

### Same-artifact A100 evaluation

`results/hdiv-elements-a100-20261013.json` records all four native FE whole
solves on the user-started NGSolve A100. Warm medians (two samples) are
12.73 ms / 34.34 ms for affine orders 1 / 2 and 1.752 s / 86.30 s for curved
orders 1 / 2. Independent true residuals remain below 7.96e-12. These are
A100 diagnostics with exact owned-buffer counters and available-memory samples,
not CPU speedup, global memory peaks or idle-historical_cpu_worker acceptance. Initial idle/DLL
failures were retained privately. All 40 recovered job files were verified
before deleting the isolated remote job; the host remains running. General
adaptive high-order assembly, images, constraints, nonlinear material, broad
physics, IH Simulink and release/distribution gates remain pending.

### Bounded host-vector cache numerical acceptance

`host-vector-cache-stage-20261013.json` tracks an opt-in `host-cache` assembly
mode. It reuses source potential at identical mapped target-host quadrature
points across target monomials. Charge weights, inner integral math, ordered
outer sums and symmetric assembly stay unchanged. Reused hosts must have identical
geometry and topology; inconsistent input fails. Temporary workspace increases
setup owned memory and is budget checked. Default `direct` stays selected.
Clean native source `57a3dda15`, HMatrix artifact `job-6c403wdx`, now passes
direct/cache whole-solve parity for affine order 2 and P2 curved orders 1 and 2.
`results/cache-curved2-20261013.json` records identical direct/cache coefficients
on both 100 and the already-running NGSolve A100. Independent Python CPU
coefficient error is 5.02e-10, physical mean error 1.49e-12 and original-equation
true residual 7.95e-12 for curved order 2 (288 DOFs, 160 charges).

A100 transfer-inclusive whole-solve warm diagnostics are 86.4778011 s direct
and 0.2696088 s cached, with one warm sample per mode, sequential rather than
ABBA. These establish neither formal idle-historical_cpu_worker performance acceptance nor CPU
speedup. Exact accounted setup/combined peak rises from 450288 to 2286980 bytes;
cache does not reduce memory. Available-memory samples exclude CUDA/cuBLAS,
host and global peaks. First payloads share initialized modules, with owned
Engine startup reported separately. All 581 remote non-environment files were
hash-recovered before guarded job deletion; environment freeze and DLL hashes
are retained privately. The A100 host remains running.

Default `direct` remains unchanged. Affine order 1 cache, broader operator and
error/lifecycle coverage, general adaptive FE semantics, formal timing and exact
release-wheel distribution acceptance remain pending. The wheel currently
excludes these experimental HDiv wrappers; staged A100 acceptance is not wheel
acceptance. AMS and FEM matrix-build work belong to the separate AMS session.

For a release candidate containing bounded HDiv, explicitly set
`RADIA_HDIV_HMATRIX_BUILD` and `RADIA_HDIV_PROJECTION_BUILD` to clean accepted
build directories before building the wheel. Both receipts must name the exact
wheel source commit and match every source and binary hash. Staging includes
the two MEX files, the three standalone MATLAB wrappers and `hdiv-build.json`;
the wheel verifier rejects a partial or mismatched payload. User-installed
compatible native dependencies remain required and are not bundled. This
packaging path has contract coverage but still needs fresh same-commit native
builds and actual installed-wheel acceptance on 100/A100. The earlier staged
artifacts deliberately cannot satisfy the new release-commit check.

### Formal A100 same-host comparison (2026-10-13)

The user designated NGSolve A100 as the official CPU/GPU measurement host.
`measure_hdiv_cpu_gpu.py` verifies the staged GPU payload and independent
Python/pybind11 reference hashes, then performs two CPU/GPU/GPU/CPU cycles
under before/after idle guards. Python CPU setup, assembly, solve and mean
recovery are compared with MATLAB GPU setup, transfer, solve and recovery.
Each sample checks independent coefficient/mean agreement and original-system
true residual. Python does not execute GPU numerical kernels.

For the independently authored public 13-tetrahedron P2 sphere, HDiv order 2,
linear isotropic material and opt-in host-cache, warmed medians are CPU
5.16058735 s and GPU 0.27303460 s (18.90x). Maximum true residual across ten
first/warm samples is 7.96e-12. Exact accounted GPU owned peak is 2286980 bytes;
available host/device samples do not measure host/global/CUDA runtime peaks.
Engine startup and first payload are separate; this is not a complete
process-cold CPU/GPU comparison. The loaded CPU native hash is recorded, but
its native build-source receipt has not been independently established.
See `results/formal-a100-curved2-20261013.json`. This bounded measurement
does not certify other geometries, nonlinear materials, general adaptive
HDiv, an installed wheel or release completion. Prior diagnostics retain
their original classification.


### Optional assembly event diagnostics

The tenth `solveHDivElements` argument is an optional logical
`profileAssembly` (default `false`). Cached assembly modes accept it; direct
assembly rejects it. The native `profile_assembly` option requires a scalar
logical and a module implementing the corresponding stats fields.
`result.assembly_profile` records mapping, tetrahedral/triangular/scalar
potentials, contraction and finalization stream-event spans. Events use the
owned stream and existing synchronization boundaries. They include launch
gaps and instrumentation overhead: these are diagnostic spans, not pure
kernel active times or a new transfer-inclusive benchmark. Runtime event
allocations are outside exact owned-buffer counts. Disabled stage times are
zero. `checkHDivFusedMoments(mesh,output,true)` checks profiled operator parity,
finite durations and handle cleanup.

The bounded 10,344-DoF linear public sphere profile attributes about 96.5%
of the measured stream spans to tetrahedral source potentials. Its independent
CPU coefficient and true-residual checks pass; no additional speedup or release
acceptance follows from instrumentation. See
`results/assembly-profile-20261018.json`.

### Separated public curved spheres

`prepare_hdiv_scale.py --out <new-directory> --maxh .6 --sphere-separation 8`
independently generates two unit-radius spheres centred at (-4,0,0) and
(4,0,0) metres. Both use `magnetic_body` and `exterior`, P2 geometry and
order-two HDiv. The saved/reloaded mesh and label contract are checked before
solving. This is a linear isotropic research case, with applied H=(1,.2,.1)
A/m and inverse susceptibility 2.5, unrelated to an engineering device.

The accepted 5,184-DoF/2,496-charge case uses explicit HACApK epsilon 1e-10,
leaf size 128, eta 2 and rank cap 128. One rank-37 block reduces total stored
coefficients from 3,222,112 to 1,759,664 versus the same upper-leaf dense
partition (45.39%). This is not a global memory reduction: the potential cache
still dominates setup, with accounted owned peak 527,796,612 bytes.
Independent CPU solution error is about 5.03e-13 and original true residual
4.20e-12. `validate_scaled_hdiv.py --operator-probes 3` also checks seeded
random B-transpose G B actions against the independent CPU operator; these
checks are excluded from performance timing and do not constitute a full
operator-norm bound. The seed and HDF5 input hash are recorded.

ABBA was stopped by the host idle guard; its partial samples are rejected for
performance acceptance. No speedup follows from this case yet. First/warm
numerical diagnostic times are not process-cold measurements. See
`results/hacapk-separated-spheres-20261018.json` for numerical acceptance,
actual storage counters, owned/sample memory limits and pending release gates.

A separate fresh two-cycle ABBA completed on the same input and native hashes:
CPU16/BLAS1 median 10.0512 s versus GPU API median 1.38544 s (7.2549x).
This includes FE setup, assembly, transfers, solve, recovery and destruction;
Engine/interpreter startup, independent validation, HDF5 and controller
transport are excluded. The CPU workflow uses mass-Riesz sparsecholesky,
while the GPU uses the explicitly selected mass-Cholesky preconditioner.
All samples pass independent solution and original true-residual gates.
The interrupted prior samples are not included. These are warmed API workflow
measurements, not process-cold or common-wheel release acceptance, and there
is no paired dense-GPU performance comparison on this two-sphere input.
See `results/hacapk-separated-spheres-abba-20261018.json`.

The subsequent **same-module GPU ABBA** compares `hacapk-cached` against
`host-cache-moments` on exactly this two-sphere input. Warm API medians are
1.38446 s and 1.04413 s respectively: dense is about 1.33x faster here.
Actual storage reduction does not imply a whole-workflow speedup. Both modes
pass the independent CPU solution and original-equation residual checks on
every sample. The cache policies differ explicitly: HACApK has a 2 GiB cap
(about 524 MB of potentials on this case), while dense retains its existing
256 MiB source-tile policy and 4 Mi independent-work batching. This comparison
therefore does not isolate compression as the sole cause of time or peak
differences. `validate_scaled_hdiv.py --gpu-pair-cycles 2 --abba-cycles 0`
selects this paired route with explicit `hacapk-cached` options. See
`results/hacapk-dense-separated-spheres-20261018.json`; dense remains a valid
control rather than being replaced by a slower compressed path.

`results/hacapk-final-entry-sync-20261018.json` records a later numerical
candidate that copies the cached-entry sticky error flag alongside the final
tile output, sharing its existing stream wait. No entries reach ACA before
that flag is checked. RTX full affine/curved operator and forced-eviction
checks pass, followed by independent CPU solution, original residual and
seeded action checks on the two-sphere case on the same compiled artifact.
This removes a redundant wait without changing arithmetic or tolerances;
whole-workflow speed improvement has **not** been established. The initial
bounded local mesh-check failure and the retry with an explicitly single-thread
NGSolve canonical checker are retained in private operational evidence.

The full-cache host-workspace candidate bypasses tile maps and intermediate
batch copies when every source column fits in one cache tile. It reuses
49,152 bytes of host vector storage at entry capacity 1,024 (allocator overhead
and other host temporaries excluded). This is separate from device-owned
buffers and is not a host/global peak. Multiple-tile and forced-eviction
paths remain available. Assembly stats expose
`assembly_cache_full_tile_batches` and
`assembly_cache_host_scratch_bytes` to verify actual use.

RTX full-operator/eviction checks and independent CPU solution, original
residual and random-action checks pass. A same-module GPU ABBA on the two
spheres gives HACApK 1.26805 s versus dense 1.06987 s; dense remains faster.
All 4,417 HACApK batches used the full-cache path, with unchanged entry,
leaf, rank and coefficient counts. Prior artifact measurements are unpaired
controls, so no causal old/new speedup is claimed. See
`results/hacapk-full-cache-host-reuse-20261018.json`. An additional analytic
CUDA regression did not execute because the isolated local CPU reference
module was absent; it remains an explicit pending gate.

### Deferred dense-leaf boundary (C stage only)

`dense-leaf-device-fill.md` describes the ownership, layout and memory changes
required to remove dense-leaf CPU round trips. The opt-in C boundary is tested
with `check_full_host.py --deferred-dense`; ordinary builds do not define
`RAD_HACAPK_DENSE_LEAF_SINK`. Eight synthetic BBOX/PCA, rank-cap and symmetric
fill conditions preserve the original leaves, permutation and low-rank
factor bytes. After test-only CPU materialization, the maximum direct-action
relative difference is 2.245e-13. Destruction before materialization also
passes. Default scalar/batched C regression remains accepted.

This is **not connected to GPU filling**, and is not a native GPU build,
transfer reduction or performance acceptance. See
`results/deferred-dense-boundary-20261018.json` for the bounded contract and
remaining gates. Existing packing rejects an unmaterialized dense payload;
no supported command publishes an incomplete matrix.


### GPU direct dense-leaf materialization (bounded research acceptance)

The subsequent explicit `hacapk-device-dense` mode completes the C-stage
boundary above. It retains the full source-potential cache until normalized
dense spans have been written directly into the resident coefficient buffer.
Only low-rank coefficients cross the host upload boundary; incomplete contexts
are never registered. Insufficient full-cache capacity is rejected. Concurrent
oracle plus Resident owned buffers are checked against 2 GiB and GPU reserve.

RTX all-column affine/curved operator comparisons agree exactly with the
same-module cached route; failed construction, stale handles and live-count
cleanup pass. On the independently authored separated P2 spheres (5,184 DoF),
original-equation true residual is 4.205e-12 and three independent CPU action
errors are at most 2.426e-12. Same-module transfer-inclusive GPU ABBA medians
are 1.091684 s direct versus 1.278653 s host round trip (1.1713x). Fresh CPU/GPU
ABBA medians are 10.008421 s / 1.093932 s (9.1490x), CPU16/BLAS1.

226 dense leaves / 1,667,312 coefficients are directly materialized; host
coefficient upload is 738,816 bytes for low-rank factors. Overlap increases
accounted setup owned peak to 543,239,548 bytes; no global memory reduction
is claimed. These are warmed API totals including construction, transfers,
solve, recovery and destruction, excluding Engine/HDF5/controller; not process
cold or general/release acceptance. Full measurements and limitations are in
`results/hacapk-device-dense-20261018.json`.


The same native artifact also passes the public 10,344 DoF P2 sphere
(4,208 charges). Independent CPU coefficient relative error is 5.234e-11,
original-equation true residual 6.268e-12, and three CPU action errors are
at most 1.860e-11. Same-module GPU ABBA medians are 2.334304 s direct versus
3.289161 s cached (1.4091x). Independent CPU/GPU ABBA medians are
21.209392 s / 2.339421 s (9.0661x), with the same warmed API boundaries.

3,179 dense leaves / 8,685,668 coefficients are directly materialized, leaving
1,584,848 bytes of coefficient uploads. Concurrent accounted owned peak is
1,791,115,980 bytes, below the 2 GiB cap; this is not a global peak or a
memory-reduction result. Stored compression is only 0.75335 percent against
the same-leaf dense representation for this sphere; 100k readiness is not
established. See `results/hacapk-device-dense-10344-20261018.json`.


### Sparse SPD mass polynomial and frozen H-matrix reuse

`mass-polynomial` explicitly selects eight fixed damped-Jacobi inverse-action
steps, with damping from a normalized row-sum spectral bound. The resulting
symmetric polynomial is SPD; it uses the existing sparse mass matvec and
linear device workspace rather than a dense factor. It is not IC, AMG or an
exact inverse. `HMatrixOperator` accepts explicit `hdiv_hacapk_device_dense`
configuration and `solveMaterial` calls with explicit mass/RHS; geometry and
charge map stay frozen. Solver workspaces are not cached.

Initial tests found a pageable H2D/default-stream upload completion boundary
missing before nonblocking stream consumption. The fix synchronizes that
upload stream before publishing its data; see NVIDIA's
[CUDA 12.8.1 synchronization contract](https://docs.nvidia.com/cuda/archive/12.8.1/cuda-runtime-api/api-sync-behavior.html).
Original Cholesky/sparse-breakdown failures are retained. No tolerance change
or CPU fallback was used. Corrected RTX affine and curved reuse/positive,
negative, zero RHS/invalid mass/lifecycle gates pass.

On the same public 10,344 DoF sphere, independent CPU coefficient relative
error is 4.508e-10 and true original residual 8.687e-12. CPU/GPU ABBA medians
are 21.151502 / 2.061307 seconds (10.2612x), with the warmed transfer-inclusive
API boundaries above. Sparse solver-owned buffers are 7,670,312 bytes versus
864,077,036 for the dense Cholesky control on the same frozen operator. The
full setup peak remains 1,791,115,980 bytes, dominated by potential caching.
Reuse solve-only samples are diagnostic and exclude construction; they are
not the whole-workflow benchmark. Different-session prior artifacts do not
provide a causal preconditioner speedup. This remains a **linear** bounded
research case. GPU nonlinear BH/Newton is pending. See
`results/hdiv-sparse-spd-reuse-20261018.json`.


### Experimental nonlinear BH material and bounded solve

`nonlinear-rtx-20261018.json` records a synthetic monotone SI BH table on the public 13-tetrahedron P2 sphere (288 HDiv DoFs). NGSolve owns IntegrationRuleSpace numbering, mapped shapes, and weights. GPU sparse contractions implement the inverse-law coenergy gradient and consistent tangent; one geometry H-matrix is reused through energy-Newton iterations. Independent CPU material energy/gradient/tangent agree below 1e-15; weak, strong and saturation-region final coefficients pass the declared 1e-8 comparison at tightened nonlinear tolerance 1e-10. The strong-case returned GPU coefficients also pass an independently assembled CPU nonlinear residual.

This is bounded experimental numerical acceptance, not large nonlinear or release acceptance. Single-call RTX times are diagnostics, not speedup measurements. Saturation takes 26 GPU Newton iterations versus 11 CPU iterations; investigate before performance claims. Gradients/tangents and coefficient vectors still cross the host boundary. MATLAB array bounds and native owned counters are not global CUDA/host peak memory. Loose-tolerance comparison failures and build failures are retained privately. `checkNonlinearSolve.m` is the reusable whole-solve comparison entry point; independent fixtures come from the CPU route, never GPU snapshots.


### Nonlinear 2412-DoF numerical gate

`nonlinear-2412-20261018.json` records independent CPU coefficient agreement on the public P2 sphere at H=(0,0,100000) A/m, using the explicit all-high-order integration profile. Development and compute GPU routes both reach coefficient relative error about 1.81e-9, mean error 2.67e-13, and nonlinear residual below 5e-15 in 44 Newton updates. The returned development-GPU vector also has independently reconstructed CPU nonlinear residual 3.46e-11. Inner linear original-equation residuals satisfy the declared tolerance. Operator handles are released. These are numerical gates and single-call diagnostics, not formal performance or release acceptance.

A failed comparison with CPU production defaults is retained: the default far_quad=3 / ho_far_factor=2 operator differs from the GPU all-high-order route and produced nonlinear coefficient error 2e-4. Matching the explicit integration profile resolved the discrepancy without relaxing tolerances. Do not call this parity with CPU default split integration. Material and RHS decomposition errors were separately checked.

The frontend now offers explicit material and geometry-cache byte budgets, each at most 1GiB, with a conservative combined bound at most 2GiB and a separate 20-percent/minimum-2GiB free-device reserve. The accepted 2412-DoF call requests 1GiB material and 256MiB cache, with computed combined bound 1327978352 bytes. BH tables occupy 38416 owned bytes and are uploaded once per solve. Material gradients/tangents still cross the host boundary. Whole timing releases material GPU owners before stopping; counters and array bounds do not establish global memory peaks. Native 65536-point/4000000-shape-entry bounds remain, so larger material maps require the chunked implementation described in `nonlinear-chunking-plan.md`.


### Chunked nonlinear material and 6594-DoF example

`nonlinear-chunks-6594-20261018.json` records the public P2 unit sphere with 112161 material points. The explicit native element-range command preserves NGSolve global sample IDs, orientation and weights, with the existing 65536-point/4000000-shape-entry bound per range. Full-map rejection remains tested. Chunk sizes 64/128 agree with independent CPU energy/gradient/tangent within 1.2e-15; small chunk sizes1/2/7/13 and invalid ranges are also checked. GPU partials are summed before gathering the complete gradient/tangent.

Whole nonlinear solves on both GPU roles agree with independent CPU coefficients within 1.41e-9, and independently reconstructed residuals are about 2.42e-11. In this H=(0,0,100000) A/m case both use two Newton corrections; this does not predict iteration counts for other meshes/loads (the 2412-DoF case took44). The explicit all-high-order CPU integration profile is required. The ideal nonlinear sphere mean error is0.01064percent and the declared81-point core magnetization RMS error is0.22501percent, below the declared1percent physical gate. This does not certify arbitrary geometry or the complete near-surface field.

The compute-GPU single-call diagnostic is7.0461s, including about3.8306s of native FE/map preparation and0.9178s of H-matrix creation. Host phase timers include synchronous API boundaries; asynchronous upload may carry into warmstart. These are not CUDA kernel timings or a formal same-host CPU/GPU speedup. Material maps are uploaded after the deferred geometry cache is released. The controlled-array bound is1988212592bytes; CUDA-library workspace and host/global peaks are excluded. Initial/repeated lifecycle checks end with zero native operator handles.

Energy fields in these experimental Newton records are the integrated coenergy merit with the common mu0 factor omitted, in A^2*m. Multiply by mu0 for joules. They are not induction-heating power or heat. `validate_nonlinear_sphere.py` is the independent CPU FE mean/core-stencil physical check for the declared public mesh. Native material chunking is accepted only for the recorded bounded examples; general nonlinear, 10k, arbitrary constraints and common-wheel release gates remain pending.

### Chunked nonlinear BH at 10,344 DoF (2026-10-18)

[nonlinear-chunks-10344-20261018.json](results/nonlinear-chunks-10344-20261018.json)
records the same public P2 unit sphere at 10,344 DoF, with 176,645 material
integration points. Chunk sizes 64 and 128 match independent CPU material
energy, gradient and tangent to 1.07e-15 relative. The complete monotone BH
Newton solve uses an explicit fixed SPD mass-polynomial preconditioner and
reuses the geometry H-matrix through its two accepted Newton steps. Both GPU
roles match independent CPU coefficients within 1.48e-9; the independent CPU
reconstruction of the returned nonlinear residual is below 2.87e-11.
The ideal-sphere mean error is 0.00863%, and the declared 81-point internal
magnetization vector RMS error is 0.121%, below the unchanged 1% physical gate.
The CPU reference explicitly uses all-high-order quadrature; this is not an
agreement claim for the CPU default near/far split.

Single-call whole API diagnostics are 32.93 s on the development GPU and
12.58 s on the compute GPU. These include FE staging, assembly, transfers,
Newton, recovery and release, but exclude Engine startup and validation.
They are not ABBA, full-process cold measurements or a CPU/GPU speed ratio.
The conservative controlled-array bound is 2,068,412,007 bytes, below 2 GiB;
CUDA workspace, host and global peaks remain unmeasured. Dense Cholesky does
not fit this bound. The cache request is explicit, and the native cache and
combined allocation limits remain enforced. Fixed FE/map preparation and the
1,375/1,739 inner iterations are remaining performance targets. This is a
bounded synthetic BH case, not general nonlinear, IH, topology optimization,
100k DoF or release acceptance.

### Prepared nonlinear reuse and bounded inexact inner solves

[nonlinear-prepared-reuse-20261018.json](results/nonlinear-prepared-reuse-20261018.json)
records the explicit `radia.gpu.PreparedHDivNonlinear` owner. FE setup, mapped
material quadrature, GPU maps and geometry H-matrix remain immutable between
solves; `delete` releases their resources. Repeated zero-start solves,
explicit-initial-coefficient solves and invalid-input recovery pass on both
GPU roles, with final native handle count zero. Independent CPU coefficient
difference remains below 1.48e-9 and independently reconstructed nonlinear
residual below 2.87e-11 on the public 10,344-DoF sphere.

`InexactNewton=true` is optional and remains disabled by default. Inner true
residuals stay below 1e-6; the outer target stays 1e-12. Relaxing early inner
work explicitly reduces warmstart iterations from 1595 to 485 and Newton
inner iterations from about 3116 to 1962 on this case. It is not universally
better: the small saturation example takes 16 Newton steps rather than 11.

Two complete one-shot ABBA cycles on the compute GPU give medians 11.3624 s
(strict) and 9.8330 s (inexact), a 1.156x whole-API ratio. Engine/GPU startup
is excluded; the first strict sample includes the first native invocation.
A separately prepared inexact context takes 7.1214 s to construct, followed
by three zero-start solves of 2.4669/2.4632/2.4644 s. Construction, three solves,
small result checks and release total 14.5261 s. The 2.46-second value therefore
must not be called a complete one-shot nonlinear analysis.

The 2,068,412,007-byte controlled-array bound is retained. Native solver-owned
buffers and available-device-memory samples are recorded; library pools,
host and global peaks are not measured. The older validation-only
`total_through_release_seconds` field stopped before release and is excluded
from whole-time claims. Failed measurement serialization and incorrect
function-path selection are retained privately. This is bounded research
acceptance, not general nonlinear, multi-language backend or release acceptance.

### Historical polynomial CUDA Graph candidate failure

[polynomial-graph-20261018.json](results/polynomial-graph-20261018.json)
retains both successful subgates and an unresolved failure. The explicit
`mass-polynomial-graph` option captures the same memset and eight CSR/update
pairs; original PCG reductions and true-residual certification remain outside
the graph. The original `mass-polynomial` option remains the default.

On identical frozen G/mass/RHS systems, graph/original solutions are identical
for positive, negative and zero RHS scalings; the nonzero cases both take
167 iterations with true residual 8.67e-12. Successful 10,344-DoF nonlinear
samples match independent CPU coefficients within 1.48e-9 and independent
CPU residual reconstruction below 2.87e-11. Actual graph launches and resource
release are checked.

However, repeated nonlinear performance measurement failed at the first graph
solve after the original warmup with an energy line-search error. No ABBA
performance result is accepted. Three bounded diagnostic repetitions did not
reproduce the failure; that does not establish a fix. The small saturation
example also took more Newton steps and ran slower with the graph candidate.
The failure cause remains unresolved, with no tolerance relaxation, rejected
step acceptance or CPU fallback. Driver graph allocation is not included in
native owned-buffer counts. Do not present this candidate as a speedup or
release-ready implementation; prior accepted non-graph records remain valid.


### Residual-certified Newton trial regression

[nonlinear-residual-first-20261018.json](results/nonlinear-residual-first-20261018.json)
records the subsequent diagnosis and bounded fix. Deterministic perturbations
of the independent CPU root exposed full Newton trials with original residual
about 1e-14, rejected by absolute-energy Armijo comparisons differing by only
1--9 floating-point ULPs. Both original and graph routes were affected.

A trial meeting the unchanged requested original nonlinear residual now precedes
energy globalization; the next outer iteration independently recomputes that
residual. All nonconverged trials retain the original Armijo condition. No
energy allowance, relaxed residual target or CPU fallback is introduced.
Five near-root perturbations per route pass on the compute GPU, together with
whole-solve CPU comparison, invalid-input recovery and handle release. Independent
Python/NGSolve reconstruction is below 2.864e-11. The earlier failed evidence
is preserved. This regression alone is not a performance or release claim;
the graph route remains opt-in. Graph-driver, library and global memory peaks
are not included in controlled buffer bounds.

Fresh compute-GPU two-cycle ABBA completed after the fix: original/graph whole
medians are 9.72653/9.65676 seconds (1.0072x). The approximately 0.7 percent
difference does not establish a meaningful speedup. Three prepared zero-start
nonlinear solves take median 2.48497/2.43500 seconds, excluding initial FE and
H-matrix construction. Construction plus three solves plus release takes
14.59695/14.49196 seconds. These are moderate-field public-sphere results,
not CPU/GPU speedup, process-cold or release acceptance.


### Direct CSC material staging

[material-csc-staging-20261019.json](results/material-csc-staging-20261019.json)
replaces one ordered-tree allocation per material evaluation nonzero with
per-column append/sort and direct MATLAB CSC construction. Duplicate local DOFs
or output entries fail explicitly. NGSolve mapped shapes, weights, global point
IDs and FE semantics are unchanged. Old/new sparse structure and values are
bitwise identical at 288, 2,412, 6,594 and 10,344 DoF; independent CPU material,
whole nonlinear solve, invalid interval and release gates pass.

For the same moderate synthetic-BH public sphere at 10,344 DoF, fresh compute-GPU
two-cycle ABBA measures complete nonlinear APIs at median **9.74162 to 5.43383
seconds (1.7928x)**. FE/material preparation falls from 5.25558 to 1.03583 seconds.
This compares old and new GPU routes, not CPU speedup. Both retain two Newton
updates; independent CPU original nonlinear residual is below 2.864e-11.
The 2,068,412,007-byte controlled device-array bound is unchanged. Host staging,
library/driver and global memory peaks are unmeasured, so no quantified host
memory reduction is claimed. First calls are separated; complete process cold,
general geometry, installed common wheel and release acceptance remain pending.


### Explicit fixed Chebyshev mass inverse action

[chebyshev-preconditioner-20261019.json](results/chebyshev-preconditioner-20261019.json)
and [SPD derivation](chebyshev-mass-plan.md) describe the opt-in
`mass-chebyshev` eight-step inverse action. Existing defaults are unchanged.
The spectral upper bound is conservative; the lower interval endpoint is a
smoothing target. Original PCG and nonlinear residual checks remain mandatory.

On identical frozen G/mass/RHS, positive and negative loads take 111 instead of
167 PCG steps, with original true residual 8.36e-12; zero load and invalid-mass
recovery/release pass. Independent CPU whole nonlinear coefficients agree
within 1.48e-9, and independent CPU reconstructed nonlinear residual is below
2.864e-11. The 288-DoF strong-saturation case needs 11 Newton updates in both
routes, while inner PCG totals fall from 6704 to 5118.

For the 10,344-DoF moderate-field public sphere, both routes need **two Newton
updates**, not thousands: the thousands are summed inner PCG iterations.
Fresh compute-GPU two-cycle ABBA whole medians are **5.51465 to 5.12746 seconds
(1.0755x)**, including preparation, transfers, nonlinear solve and destruction.
After retaining FE/material/H-matrix preparation, three zero-start nonlinear
solves have medians **2.47946 to 2.08163 seconds**; those exclude preparation.
Construction plus three solves plus release totals 10.48438/9.27957 seconds.
The extra direction buffer is 82,752 bytes at this size, included in native
solver accounting (7,753,352 bytes). The existing controlled combined bound
remains 2,068,412,007 bytes; host/library/driver/global peaks are unmeasured.
This is a bounded GPU/GPU comparison, not CPU speedup, process cold, 100k-DoF,
RadFld projection, topology-optimization or common-wheel release acceptance.


### Residual-only final recertification

[final-residual-only-20261019.json](results/final-residual-only-20261019.json)
validates explicit `FinalResidualOnly=true` (default false). After an accepted
trial satisfies the requested original residual, the next outer check recomputes
the material gradient, energy and demagnetizing operator without constructing an
unused tangent. The trial residual is never reused as convergence proof. If the
new residual fails, a tangent is constructed and Newton continues. Targets,
Armijo conditions and native kernels are unchanged.

For 10,344 DoF, two Newton updates need two instead of three outer tangents;
the strong-saturation 288-DoF case retains eleven updates and needs eleven
instead of twelve tangents. Both pass CPU comparison and lifecycle gates.
Fresh compute-GPU two-cycle ABBA whole medians are **5.07906/4.91936 seconds**
(1.0325x); prepared nonlinear solve medians are **2.07351/1.93075 seconds**,
excluding preparation. The small whole-time improvement does not establish a
universal benefit. Independent CPU original residual is below 2.864e-11.
There is no extra retained buffer and the controlled memory bound is unchanged;
omitting a final tangent does not prove a lower whole-run peak. Library, host
and global peaks, 100k-DoF scaling and release acceptance remain unmeasured.


### Strong synthetic excitation at 10,344 DoF

[strong-nonlinear-10344-20261019.json](results/strong-nonlinear-10344-20261019.json)
raises the applied field to 1e6 A/m on the same public mesh and synthetic table.
The declared finite-barrier continuation above the final magnetization knot
is active (maximum M about 991740 A/m, cap about 990099 A/m). This is a numerical
stress case, not measured steel or a within-table-only saturation case.

The independent CPU reference converges with fourteen Newton updates, while
the explicit Chebyshev/inexact-Newton GPU route converges with seven updates
and 13326 inner PCG iterations. CPU coefficient difference is 3.202e-10;
independent CPU reconstruction of the GPU original nonlinear residual is below
9.365e-12. Mean magnetization differs from the declared ideal-sphere/barrier-law
root by 2.198e-6 relative. Repeated solve, invalid-initial-state recovery and
handle release pass. Different iteration controls prevent interpreting the
Newton-count difference as hardware speedup.

The compute GPU was unreachable before verified staging. Its strong-case
timing remains pending; the single local GPU call is diagnostic, with no
CPU/GPU performance ratio. No host was started or restarted and no size/memory
limits were raised. The earlier moderate-field compute-GPU timings remain
separate evidence.


### Large-mesh capacity probe: no large solve

[large-capacity-preflight-20261019.json](results/large-capacity-preflight-20261019.json)
and [preflight_hdiv_scale.py](preflight_hdiv_scale.py) use NGSolve-owned uniform
refinement, canonical checked labels/P2 mapping and FE-space metadata. No Gram,
material evaluation matrix or GPU buffer is allocated. The accepted 10,344-DoF
mesh passes the interpreted current capacity checks. One refinement produces
**78,456 DoF, 4,120 tetrahedra, 25,072 charge basis functions and 1,413,160 material
points**. Both source and refined meshes pass the canonical mesh check.

The refined case exceeds current element/DOF/source caps. Its existing material
array formula gives a **4.204 GB conservative bound**, above the current 1 GiB
material budget. A completely dense upper Gram would have 2.515 GB of FP64
coefficients; actual compressed storage at that mesh is unknown. Above-2GiB
exchanged payloads require HDF5, which does not itself establish GPU capacity.

The accepted coarse sphere has 32 low-rank leaves among 3211 stored leaves and
saves about 0.75 percent of coefficients against the same leaf structure made
fully dense. Do not extrapolate that compression ratio to a refined mesh.
Sparse FE staging also still encounters guards derived from hypothetical dense
mass/charge maps, so failing them is not proof of a sparse-memory requirement.
Before lifting caps, replace those guards by verified sparse/leaf allocation
bounds, bound or stream material maps, validate C index products and account
for geometry-cache tiling. This probe establishes neither 100k-DoF solve support
nor a time estimate, physical FE accuracy or release acceptance.


Sparse FE allocation admission (2026-10-19): [result](results/sparse-staging-guard-20261019.json) records a clean native build and exact old/new mass, charge map, RHS, mean and geometry parity at 288, 2,412, 6,594 and 10,344 DoF, across sparse mass with dense or sparse charge maps. Invalid-order recovery and zero live handles passed. Sparse routes now use overflow-checked local FE contribution bounds; dense allocations retain their previous limits. Element/DOF/charge caps are unchanged. This is FE staging regression evidence, with no new solve residual, speedup, larger-capacity or release claim.

The same sparse-admission binary also passed the bounded 10,344-DoF nonlinear solve and lifecycle regression: [independent result](results/sparse-staging-nonlinear-20261019.json). Two Newton updates, coefficient relative difference 1.48e-9, and independently reconstructed original nonlinear residual at most 2.87e-11 were confirmed. Synthetic-sphere mean error is 8.63e-5. Times are single-call RTX diagnostics, with no new performance or larger-capacity claim; A100 remains pending.

Material chunk admission now checks its monotonically increasing partial array allowance before retaining each staged host chunk. [Small-case regression](results/material-chunk-admission-20261019.json) confirms the dedicated low-budget error, zero live handles and successful subsequent strong-saturation solve (288 DoF, 11 Newton updates, CPU coefficient difference 4.75e-14). The current native chunk still needs bounded staging; no global peak reduction, larger-capacity or performance claim is made.

[Immutable material staging context](results/material-staging-context-20261019.json) reuses an NGSolve-owned loaded mesh, HDiv space and IntegrationRuleSpace across chunk requests. Old path-based and new context-based outputs are exactly equal at 288/2,412/6,594/10,344 DoF. Wrong-type, invalid-range, destroyed and reloaded handles are rejected; contexts are explicitly owned and counted. Chunked nonlinear preparation now releases its context after staging, including the early-budget error path. A 288-DoF strong-saturation solve/repeat/recovery passed CPU parity; full 10k integration and transfer-inclusive performance remain pending. File changes require context recreation; this is no measured speedup claim.

The immutable material context is also [checked in the full 10,344-DoF nonlinear pipeline](results/material-context-nonlinear-20261019.json): two Newton updates, frozen independent CPU coefficient difference 1.48e-9, native original nonlinear residual below 1.25e-14, repeat/invalid-initial recovery and zero live handles. The controlled combined array bound remains 2,068,412,007 bytes. These RTX samples are ordered diagnostics; they establish no speedup. A100 transfer-inclusive ABBA remains pending.
