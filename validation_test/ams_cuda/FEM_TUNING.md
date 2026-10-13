# FEM/AMS tuning priorities

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
records prioritization only; no new FEM speedup is certified here.
