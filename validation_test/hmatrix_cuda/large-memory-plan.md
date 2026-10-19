# Large HDiv memory work: admission before numerical allocation

Status: source investigation and implementation plan only. No cap change,
large operator, large nonlinear solve, speedup or release acceptance.
The canonical counts are in [capacity preflight](results/large-capacity-preflight-20261019.json).
The currently accepted nonlinear size remains 10,344 DoF.

## Measured inputs and distinct bounds

The one-refinement public P2 sphere has 78,456 DoF, 25,072 charges,
4,120 tetrahedra, 1,413,160 material points and 127,184,400 potential
mapped shape entries. No full material evaluation map was allocated by
that preflight. The material tangent contribution bound is 3,708,000.

| Quantity | Bytes | Interpretation |
| --- | ---: | --- |
| Current material array allowance | 4,204,297,600 | Conservative formula, not measured peak |
| Current resident allowance | 10,124,791,808 | Conservative square bound, not compressed storage |
| Current solver allowance | 581,902,336 | Conservative formula, not measured peak |
| Full upper-triangle FP64 charge coefficients | 2,514,521,024 | Only if all represented coefficients were dense |
| Potential mapped shape values alone | 1,017,475,200 | 8 times the potential entry count; excludes all indices/workspace |

The last row is NOT a proposed implementation's total memory. Replacing
24-byte sparse allowances by 8-byte value storage alone would ignore local
indices, accumulation structure, tangent workspace and library allocations.
The coarse sphere saved only 0.75335% of same-leaf dense coefficients;
that observation cannot predict compression on the refined mesh.

## Leaf-frame admission insertion point

`src/ext/HACApK/cHACApK_cpp_impl.c`, function
`HACApK_build_hmatrix_wrapper`, already separates cluster construction,
`cHACApK_count_lntmx`, leaf allocation/generation/sort, and numerical
`cHACApK_fill_leafmtx_hyp`. The candidate insertion point is after the
sorted frame and before any ACA numerical fill. This is not implemented.

Before allocating a frame, retain a checked signed-product/leaf-count bound.
After generating it, a GPU-build-only admission hook could inspect dimensions,
starts, admissibility and rank caps without evaluating any entry. Default CPU
behavior must remain unchanged. A rejected hook must unwind every cluster,
frame and control allocation; it cannot return a partially accepted plan.

For the stored upper/diagonal leaves, compute overflow-checked sums for dense
coefficients, bounded ACA factors, maximum serial ACA workspace and metadata.
A cap-exhausted ACA leaf must budget its existing exact-dense promotion, not
assume compression succeeds. A conservative per-admissible-leaf value bound
must include the larger of dense storage and bounded rank-factor storage.
`cHACApK_base.c::fill_one_leaf_block` skips strictly lower leaves when symmetric
fill is enabled, but their frame metadata still exists. Diagonal/mirror rules
and existing pivot/normalization arithmetic stay unchanged.

The C temporary factors, packed host factors, resident CUDA coefficients,
sparse B/B-transpose, device dense-fill workspace and retained oracle/cache
have overlapping lifetimes. Admission must describe those phases and prove
all sums and byte products before allocation. Existing `checkedProduct` and
`Resident` reserve checks remain mandatory, including its separate 2-GiB
combined deferred-setup constraint. A tighter leaf bound alone does not
justify removing that constraint or increasing any source/DoF cap.

## Material representation candidate

NGSolve continues to supply mapped shapes, quadrature weights, orientation
and element DOFs. An element-local value buffer plus checked DOF maps could
avoid repeated generic sparse-matrix metadata and SpGEMM temporaries.
The kernel must implement the SAME inverse-BH energy, gradient and symmetric
tangent, including finite barrier continuation. It must not interpret FE
coefficients as nodal values or replace the mapped quadrature.

First construct symbolic element-to-CSR scatter once using NGSolve-owned DOFs.
For deterministic reduction, accumulate local contributions into bounded
workspaces and reduce each CSR entry in a declared fixed element order.
Do not use unchecked global atomic accumulation or alter the bilinear form
merely to obtain a different convergence history. Keep FP64 external state,
original-equation residual certification, error-state isolation and explicit
release. Streaming chunks from host on every Newton evaluation is a separate
transfer-cost tradeoff, not an assumed optimization.

## Acceptance sequence

1. Keep current element/DoF/charge/material caps unchanged while implementing
   and testing symbolic admission and allocation/index overflow rejection.
2. Prove exact staged NGSolve data identity and independently compare material
   energy/gradient/tangent on the existing 288/2.4k/6.6k/10k cases. Exercise
   failure cleanup and accepted-state retention before any cap increase.
3. Require CPU coefficient agreement and original nonlinear/linear residuals;
   retain previous failed evidence and actual same-artifact identity.
4. Only after total host/device phase bounds and memory reserve fit, expose a
   bounded larger case. Use HDF5 for exchanged payloads exceeding 2 GiB; HDF5
   capacity does not imply GPU capacity. Never split to evade that rule.
5. Measure fresh same-host fixed-thread CPU/GPU transfer-inclusive whole time,
   first and warmed samples and owned versus sampled/peak-memory limits.
   Neither the current RTX geometry diagnostic nor the 10k result predicts
   a 100k runtime. Existing A100 strong-case and new-context ABBA remain pending.

## Implemented coefficient primitive (not frame admission)

`packages/radia-simulink/native/hmatrix_leaf_bounds.h` now provides a transactional,
overflow-checked accumulator for validated leaf extents. It counts all frame
leaves, skips lower values, and separately sums upper dense, factor-cap and
maximum-of-dense-or-factor storage coefficients. The maximum serial candidate
factor size is separate. Overflow and invalid extent rejection preserve the
previous accumulated bounds. The [native contract result](results/leaf-coefficient-bounds-20261019.json)
checks an analytical block partition, rank-cap dense promotion, invalid extents
and byte overflow with supported MSVC 14.44.

At the primitive-only checkpoint this helper was not connected to the C leaf-frame builder or CUDA allocation
admission; the bounded connection is recorded below. It neither budgets all ACA work arrays nor replaces existing memory
limits. No production MEX was rebuilt for it and no capacity cap changed.
The next connection requires a GPU-only frame hook with tested unwind paths,
followed by exact same physical operator/solver and lifecycle gates.

## Implemented opt-in C rejection boundary (GPU consumer pending)

`RAD_HACAPK_FRAME_ADMISSION` now enables a borrowed-frame C callback after
uniform-charge leaf generation/sort and before numerical fill. A zero return
frees the unfilled leaf frame, tree and temporary coordinates, returns status
-2 and publishes no matrix ownership. The outer caller still frees its control
arrays. Ordinary CPU builds do not define this macro.

The [full-host contract](results/frame-admission-boundary-20261019.json)
rejected three builds before any entry evaluation, then accepted two builds
and compared their synthetic operator to direct evaluation (3.36e-13 relative).
Default scalar and batch CPU routes independently passed the same BBOX/PCA
operator regression. This is not allocator/global-peak certification.
At that C-boundary-only checkpoint, the production CUDA callback and coefficient helper were pending and the mexcuda builder did not enable the macro. See the subsequent bounded connection below. No cap changed.

## Connected bounded mexcuda path

The full-host mexcuda build now enables the frame callback and uses the checked
coefficient accumulator in the module-local entry binding. Errors are retained
inside C++ and rethrown by the owner after C unwinds; exceptions do not cross C.
The coefficient-only allowance is 2 GiB. Previous source/index/square limits,
oracle/Resident phase bounds and device reserve checks remain active.

[RTX numerical evidence](results/cuda-frame-admission-20261019.json) covers
288 and 10,344 DoF full nonlinear solves, CPU coefficient agreement, native
original residuals, repeat/error recovery and zero handles. On the 10k case,
the actual frame has 6,325 leaves, upper dense count 8,951,208 and conservative
upper storage cap 9,059,268 coefficients; actual stored coefficients are
8,883,774. These are coefficient counts, not total memory or peak reduction.

No cap was increased. At this first connected-build checkpoint, deliberate CUDA TLS rejection was still pending; it is tested in the following section. Expanded-size admission, actual host/global peak and fresh A100 performance remain pending. Prior plans above describe work still
required before replacing square-derived admission or enabling larger models.

## Real module rejection and recovery

The native analytic HACApK configuration now accepts optional
`frame_coefficient_budget_bytes`, an integer double in 1..2 GiB with unchanged
2-GiB default. This lowers the coefficient allowance only; it does not bound
normalization/oracle/setup memory or lift any dimension/phase reserve limit.

[Real-module rejection evidence](results/frame-budget-rejection-20261019.json)
uses a one-byte coefficient allowance to reach the TLS callback and unwind the
actual C frame three times, with precise callback error and zero published
handles, then successfully assembles/applies the same independent CPU-referenced
operator. Five invalid budget values are separately rejected. Ten resident
cases, three material solves, reload/stale-handle checks, and a 288-DoF strong
nonlinear solve/repeat/recovery also pass. This closes the deliberate module
rejection test gap; it does not certify larger capacity, global memory peak,
performance or the new binary's full 10k/A100 acceptance.


Borrowed factor packing (2026-10-19): clean native source 455daee04 passes the bounded RTX resident suite, real frame-budget rejection/retry, and 288-DoF nonlinear CPU coefficient comparison. Native factors are borrowed only during synchronous packing into one reserved owned Plan buffer. This removes an intermediate host factor copy; GPU storage is unchanged. No measured whole-time or peak-memory reduction, 10k acceptance of this binary, or larger-capacity claim follows. See results/borrowed-factor-pack-20261019.json.


Subsequent full-solve gate: the same 455daee04 binary passes 288-DoF strong and 10,344-DoF moderate nonlinear frozen-CPU coefficient comparison, native residual recertification, prepared reuse and failure/recovery. See results/borrowed-pack-full-solve-20261019.json. This extends numerical acceptance only; earlier timing and memory-saving limitations remain.

## Host lifetime composition still required

The current borrowed-pack change does not make the host and device phases
mutually exclusive. Admission must account for these simultaneous owners:

| Phase | Live storage to include |
| --- | --- |
| Before frame callback | C cluster tree, control arrays, coordinate arrays, `lodfc`, leaf pointer array including sentinel and every leaf descriptor |
| Serial numerical fill | Existing retained native factors, current output factors or dense promotion, ACA masks/rank workspace/reference vectors, callback staging and resident oracle/cache |
| `fromHacapk` packing | Native factors still owned by `NativeOwner`, borrowed descriptors, final owned Plan coefficients, direction/upload/span arrays, row contribution vectors and permutation/scale |
| Resident construction | Plan still alive, final device coefficients, directions, maps, input/tmp/partial/state buffers, concurrent oracle/cache and library workspace |

`NativeOwner` is destroyed when `assembleFullEntries` returns. Borrowed views
therefore remove one copy but do not remove native-factor/Plan overlap.
`Resident(const Plan&)` likewise cannot treat the Plan as already freed.
For deferred dense leaves the host Plan stores no dense values; device
coefficient capacity still includes their exact physical dense spans.

The tested `serialAcaAuxiliaryBytes` helper bounds only masks, rank workspace
and the larger of ACA row bounds or ACA+ reference vectors. Output factors
are excluded to avoid hiding their ownership and promotion overlap. The
frame pointer count is `nlf+1`, not `nlf`; descriptor storage covers all C
leaves even when lower-triangle coefficient values are skipped. Account for
container capacity, allocator overhead and library limits separately before
calling any composition a host/global peak bound. No larger admission or
performance acceptance follows from this source audit.

### Coefficient and serial ACA coexistence guard (source candidate)

The frame callback now stages an additional 2-GiB working-payload ceiling.
The fill bound is retained upper-leaf storage plus the largest serial
rank-cap factor scratch plus the largest ACA/ACA+ auxiliary scratch. The
packing bound is twice retained upper-leaf storage because native factors
and the owned Plan coexist. The maximum of these two bounds is checked
before numerical leaf fill. Deferred dense values are conservatively
included even when their host copies are absent. The existing configurable
single-copy coefficient allowance remains a separate check.

The pure C++ contract checks exact hand-calculated examples, dense promotion,
the 2-GiB boundary and overflow in both scratch and overlapping copies.
This callback change is not yet CUDA-built or numerically accepted. It
does not cover C tree/control metadata, container capacity, callback staging,
oracle/cache, allocator overhead or the Plan/device/library coexistence.
It does not authorize higher dimension caps or establish a total peak bound.

## Requested next performance targets (not achieved evidence)

The next research targets are 100k-DoF nonlinear HDiv-MMM within 600 seconds,
including FE preparation, H-matrix construction, transfers, nonlinear solve
and recovery; the same target for polynomial order greater than two; and
integrated RadFld evaluation at approximately 1,000 target points in a few
seconds including construction/transfers/recovery. These are targets, not
extrapolated benchmark results. Record the exact geometry, order, material
law, excitation, tolerances, nonlinear/inner iterations and timing boundary.

Large 100k-DoF CPU completion runs are not required by the user. Preserve
small independent CPU agreement, original-equation nonlinear residual,
physical consistency and discretization-convergence checks. Do not replace
the original residual with a compressed or preconditioned residual alone.
Current full-solve order acceptance is one/two; higher-order staging does
not establish a higher-order nonlinear solver.

The memory target is below the 40-GiB device capacity. Maintain the existing
20-percent reserve: plan at most 32 GiB of simultaneously live allocations,
with earlier failure if measured free capacity is lower. Controlled-buffer
accounting, device-wide samples, allocator/library exclusions and host
memory limits remain distinct. HDF5 is required for exchanged payloads
above 2 GiB. Neither this target nor the CPU comparison waiver lifts current
geometry, index, dimension or admission limits.
