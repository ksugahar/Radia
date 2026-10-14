# Induction-heating solver development handoff

## Scope and publication boundary

This document covers reusable induction-heating solver development only. Project-specific identities, geometry, CAD, meshes, material datasets, operating conditions, outputs, comparisons, and reciprocal reviews belong in private project storage. Anonymizing project data does not make it publishable.

Public changes must use neutral electromagnetic domain names and independent validation against analytical solutions, public sources, or publishable original cases. Do not copy private case fixtures or numerical references into tests, documentation, packages, or CI.

## Solver baseline

Scalar P1 single-layer/double-layer HACApK entry-on-demand construction and bounded weak genus-1 support are implemented on main. See `src/radia/bem_sibc_solver.py`, `src/radia/bem_loop_work.py`, and the [workpiece contract](workpiece_holes.md). Numerical evidence and its limitations are documented in the [validation workflow](../../validation_test/induction_heating/README.md). Release acceptance must bind the current source and native binary; do not reapply completed work from historical branches.

Contracts to preserve:

- Explicit regular integration degree and singular integration order.
- Bounded shared pair-block cache, including a cache-disabled route.
- Sparse mass/operator assembly and reusable sparse factorization.
- Automatic dense-entry construction through 7,000 total mesh nodes; on-demand construction above that threshold. Count all nodes, including interior nodes.
- Explicit construction overrides with errors for unsupported combinations; no implicit numerical fallback.
- Extended weak genus-1 support limited to 20,000 faces and 10,000 total-node degrees of freedom, requiring valid native on-demand P1 single-layer and double-layer handles and the actual P0 FMM route.
- Existing strong, dense, direct-P0, P2, and application-artifact limits remain separate contracts.

The limits apply to the selected operator route, not every IH model or backend.

## Cohomology and independent FEM validation

A genus-1 workpiece requires one circulating-current cohomology unknown in
addition to the single-valued scalar potential. The production route detects
topology, constructs and certifies an interior carrier, projects its Neumann
trace into the P1 space, and closes the circulation using distributed electric
and magnetic work. It does not fill the bore or silently discard circulation.
Supported automatic carriers surround the z axis; arbitrary holes and multiple
handles are outside this contract.

The [synthetic bored-workpiece study](../../validation_test/induction_heating/bored_workpiece/README.md)
compares this BEM route with independent air-domain HCurl FEM using linear SIBC.
It includes three mesh levels, carrier sensitivity, frozen-circulation
diagnostics, true residuals and power balances. The physical comparison gate
is 2%; it is separate from solver convergence and compression error.
Both routes share the SIBC assumption, so agreement does not certify a
skin-resolved volume model, nonlinear materials, or general geometry.
Historical saved results remain historical evidence; candidate acceptance must
record the current executed source/native identities and fresh numerical JSON.

## Interpretation and remaining work

Entry-on-demand trades construction time for lower memory use. Do not describe it as a demonstrated speedup or claim linear memory scaling. Working-set and committed-memory measurements are distinct. Native allocations and shape-dependent near blocks/factor fill limit generalization.

Shared pair evaluation and batched ACA evaluation are future optimization candidates, not unfinished blockers in the completed baseline. Large strong formulations and expanded application limits require separately defined scope and validation.

## Resume and acceptance

1. Confirm the execution host, checkout, branch, remote main, working-tree ownership, and CI state. Preserve other sessions' staged and unstaged changes.
2. Use a clean integration worktree for reviewed solver changes. Verify source SHA, native sidecar/binary hashes, runtime, and input provenance.
3. Run focused deterministic contracts and independent numerical checks on public inputs. Save numerical validation JSON under `validation_test/`; retain operational logs outside tracked source.
4. Run builds and numerical validation on the designated build/test host in isolated environments. Formal performance acceptance requires an otherwise idle timing host under current repository policy. Historical timings retain their original provenance.
5. Preserve physical and residual acceptance limits. Compare vector fields using field norms; do not tune tolerances to observed results.
6. Verify current CI before claiming completion. Reconnection and live import-path checks are required for affected editable MCP clients.

Python/NGSolve owns CPU execution; MATLAB mexcuda owns GPU execution. No implicit CPU/Python fallback is permitted. MATLAB/GPU and application-block changes belong to their independently released distribution and acceptance workflows.
