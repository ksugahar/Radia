# Integration ownership audit — 2026-10-11

Scope: the CPU/GPU distribution boundary, not a complete Radia API inventory.
The user selected complete ownership transfer within the existing monorepo.

| API family | Classification | Owner and action |
| --- | --- | --- |
| CPU electromagnetic formulations and NGSolve assembly | method-keep | Keep in `radia`, callable through Python/pybind11. Preserve numerical contracts. |
| CUDA CSR PCG and MATLAB GPU entry points | user-intent | Move canonical sources, build and validation ownership to `packages/radia-gpu`. No CPU solver fallback. |
| Simulink blocks and lifecycle adapters | user-intent | Move to `packages/radia-gpu`; keep domain solver logic in its CPU or GPU owner. |
| LTspice circuit integration and command-line entry points | user-intent | Move implementation, dependencies and entry points out of `radia` into `radia-gpu`. External LTspice remains user-installed. |
| Cross-package MCP tools and operational knowledge | user-intent | Keep in `radia-mcp`; optional GPU backend must not prevent discovery on machines without the GPU distribution. |
| Independent Optuna MATLAB APIs | user-intent | Preserve `radia-optuna` ownership and release boundary. Do not silently absorb into `radia-gpu`. |
| Legacy CPU Radia MEX build route | deprecated-drop | Already unsupported and removed from the supported build dispatcher; do not restore during migration. |
| LTspice compatibility package shells in the CPU distribution | deprecated-drop | Migrate consumers before removing from the CPU package; do not leave integration aliases in the CPU wheel. |

## Migration and acceptance

1. Establish package-owned canonical sources; update build paths and consumers without duplicate copies. Preserve MATLAB namespaces where ownership changes alone do not require an API rename.
2. Remove integration dependencies, CLI exports and package data from CPU distribution metadata. Separate reusable CPU domain methods from integration adapters.
3. Update MCP routing with explicit optional-backend errors and preserve cross-package knowledge. Update CI selectors to distribution ownership.
4. Verify CPU wheel and source-distribution boundaries, GPU wheel inventory and provenance, and MCP discovery without the GPU package installed.
5. Rebuild the exact final commit on 100 and repeat installed CUDA, Simulink and LTspice acceptance before PyPI publication. A100 acceptance remains pending while unreachable; never deploy GPU to 102.

## Current evidence and limits

PR #339 remains a draft. Its pre-migration installed acceptance passed, but that does not certify the new ownership boundary. The existing candidate still stages MATLAB sources from the repository root and retains LTspice Python implementation in the CPU source tree. Publication is held until migration and the corresponding checks pass. No claim is made that AMS, ICCG or HDiv-MMM GPU support is certified by the Poisson test.
