# FEM/AMS tuning

The AMS session returns its implementation focus to FEM/AMS following the
2026-10-13 user instruction. The completed FMM candidate remains separate.

Evidence: `results/formal-a100-20261013.json`, four warmed native GPU calls on
the explicitly authorized NGSolve(A100) host. Medians of individual phases:

| Phase | Seconds |
| --- | ---: |
| Native solve total | 0.515404 |
| AMS linear solves | 0.355292 |
| Hierarchy construction/update | 0.089305 |
| Initial FE construction | 0.062191 |
| NGSolve host basis/topology staging | 0.060475 |
| Upload and GPU element assembly/reduction | 0.001126 |
| Nonlinear residual/tangent evaluations | 0.001975 |

Phase medians need not sum to the median total. This is the existing bounded
5,326-DOF case, not a fresh measurement or a large-mesh scaling result.

First investigate resident PCG host synchronization and sparse coarse solves:
`ams_resident.cuh` downloads curvature controls and a residual norm each
iteration, while `sparse_coarse` preserves factor dependencies with serial
row traversal within one warp. Neither source inspection alone establishes
their individual wall-time contribution. Measure them before selecting a
replacement. Graph-cached preconditioning and workspace reuse already exist;
do not count these as new improvements.

Next investigate mesh/basis/topology reuse across solves and AMS hierarchy
construction. NGSolve must continue to own FE semantics. Any cache must have
checked ownership/lifetime and explicit invalidation for geometry, labels,
space, coefficients and device changes; do not add an implicit global cache.

Acceptance for a new candidate requires same-input baseline/candidate timings,
independent Python CPU operator/solution checks, the original true-residual
gates, nonlinear convergence and handle/error/repeated-run checks. Report
cold construction, repeated RHS/parameter updates, transfers and owned memory
separately. Coordinate 100 builds and A100 timing with radia-gpu. This note
records the original prioritization; the candidate measurements follow.

## Implemented candidate and A100 acceptance

Candidate native source `cb1dd5add` removes the PCG curvature download: an
FP64 device kernel consumes rho and pAp and marks an invalid curvature with a
nonfinite residual, which the immediately following host norm rejects. The
regular iteration now has one host synchronization for its residual norm;
initial/final and convergence-triggered true-residual checks remain. This is
not a claim of completely asynchronous or communication-free PCG.

For sparse coarse factors with 1..512 used rows, construction solves identity
columns on GPU using the original validated NGSolve sparse LDLt factor and
dependency order. The cached inverse action is reused in preconditioner graph
replays. Larger factors keep the original GPU sparse factor kernel. There is
no CPU fallback or PARDISO. The cache changes floating-point summation, so
independent residual/solution acceptance is required. It adds at most
8,396,800 bytes for four cached coarse spaces in one hierarchy. During a
hierarchy replacement old and new caches may coexist; admission checks use
actual free device memory with the existing 2 GiB/20% reserve.

`radia.gpu.HCurlSystem` retains one linear FE matrix, RHS, native AMS hierarchy,
PCG workspace and graph. A solve changes only the explicit current multiplier:

```matlab
system = radia.gpu.HCurlSystem(volPath,boundaryLabels,SetupThreads=1);
cleanup = onCleanup(@()delete(system));
[x1,info1] = system.solve(CurrentScale=1,Tolerance=1e-12,ReturnOnGPU=false);
[x2,info2] = system.solve(CurrentScale=2,Tolerance=1e-12,ReturnOnGPU=false);
```

Geometry, labels, mass and reluctivity are fixed and private. Construct a new
object to change them; deleting the object invalidates it. This API is linear
only. `solveHCurl` retains its existing nonlinear adaptive tangent/hierarchy
updates; a new nonlinear mesh/session cache is not claimed.

Fresh A100 testing used the same checked 3,937-tet/5,326-DOF cube and original
accepted baseline MEX, two baseline/tuned/tuned/baseline cycles, one CPU/host
setup thread and unchanged 1e-12 linear / 1e-9 nonlinear tolerances. Every
sample passed an idle CPU/GPU/free-memory guard. All API timings include
required transfers and host output recovery, excluding fixture loading,
MATLAB/GPU initialization and result-file serialization. The baseline also
uses the new MATLAB reuse wrapper, so its native kernel comparison is fair.

| Median API time | Baseline | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| Nonlinear FEM/AMS including construction | 509.40 ms | 286.41 ms | 1.78 |
| Fresh linear FEM/AMS including construction | 137.34 ms | 117.13 ms | 1.17 |
| Prepared linear system construction | 110.32 ms | 105.48 ms | 1.05 |
| Reused linear solve and output recovery | 31.27 ms | 9.68 ms | 3.23 |

The candidate reused solve is 12.10 times faster than its own fresh linear
pipeline, after paying the separate 105.48 ms preparation cost. Reuse samples
change current scales between 2, 0.5 and 1; the first graph-building solve and
zero-RHS shortcut are excluded from this warm median. All six calls, including
zero RHS and first solve, remain in the JSON. The nonlinear solve phase alone
improved from 348.73 to 107.17 ms, while hierarchy setup increased from 84.85
to 91.39 ms. No separate ablation attributes these gains to only one change.

Independent Python NGSolve checks found maximum nonlinear true relative
residual 4.36e-12 and solution relative difference 1.60e-14. All reused linear
solutions passed independent sparsecholesky solution and operator checks;
maximum true relative residual was 8.41e-13. Additional 100/RTX5000 tests
checked FE matrix/mass/gradient/RHS and volume-weighted vector curls. The
existing resident contract on a 251,896-DOF fixture passed stale/wrong-type
handles, failed-update preservation, source ownership, zero RHS and workspace
reuse, plus new NaN RHS, nonpositive-curvature rejection and recovery checks.
Those 100 measurements are not formal timing evidence.

Added coarse action caches were 2,536,300 bytes for the reusable linear
hierarchy and 2,986,180 bytes for the final nonlinear hierarchy. Persistent
PCG workspace was 213,064 bytes; warmed solves allocated no native buffers.
Device-free endpoints and 50 ms MATLAB process RSS/private samples are saved.
MATLAB's allocator retains/reuses memory across variants, so endpoint deltas
are not exact live allocation or global peak measurements. The previous
`max_phase_owned_bytes_upper_bound` excludes these new caches and is not a
total pipeline bound; use the separately reported cache/workspace bytes and
the documented sampling limitations. No total-memory reduction is claimed.

`results/ams-tuning-a100-20261013.json` contains all samples and identities.
This candidate has not replaced the radia-gpu release owner's fixed wheel.
Broader materials, high-order/curved elements and large-scale A100 timing
remain outside this bounded acceptance.


## Installed wheel acceptance

The CI-built `radia-gpu 0.1.4.dev0` wheel from source
`a88e158410b330279a5cb21354c264108db9cec7` passed installed acceptance on
RTX 5000 and A100 with identical wheel, manifest and loaded AMS MEX hashes.
See `results/installed-wheel014-cross-host-20261013.json` and the adjacent
per-device records. The wheel remains a development candidate.

Eight idle A100 ABBA contexts measured whole API calls including setup and
transfers: nonlinear median 576.462 to 322.031 ms (1.79x), fresh linear
142.117 to 125.756 ms, prepared setup 122.598 to 117.501 ms, and warmed
matrix/hierarchy reuse 32.746 to 12.239 ms (2.68x). Initialization and first
API calls are recorded separately. RTX timing remains diagnostic.

Independent Python CPU checks passed. Cross-device coefficient differences
were at most 3.91e-14. Cache/workspace memory and sampled process memory were
rechecked; the sampling limitations above still apply. No global peak or
large winding scaling certification is claimed.

The A100 module identity probe required an explicit process-scoped
PowerShell execution policy. The validation controller was corrected; the
wheel was unchanged. The first failed probe and all successful samples were
recovered with hashes before job-owned remote cleanup.
