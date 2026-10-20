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
5. Rebuild the exact final commit on 100 and repeat installed CUDA, Simulink and LTspice acceptance before PyPI publication. Develop radia-gpu on 100; final accuracy acceptance compares against Python/NGSolve on gpu evaluation worker; formal GPU timing uses the gpu evaluation worker A100 when online. CPU-only reference execution does not establish CUDA acceptance. A100 acceptance remains pending while unreachable; never deploy GPU to 102.

## Current evidence and limits

PR #339 remains a draft. Its pre-migration installed acceptance passed, but that does not certify the new ownership boundary. The pre-transfer candidate staged MATLAB sources from the repository root and retained LTspice Python implementation in the CPU source tree. The accepted GPU surface now stages package-owned sources. Publication is held until migration and the corresponding checks pass. No claim is made that AMS, ICCG or HDiv-MMM GPU support is certified by the Poisson test.

The accepted CUDA PCG and LTspice MATLAB adapters, CUDA build helper, Python LTspice implementation and tests have moved to `packages/radia-gpu`. The LTspice MCP transport moved to `radia-mcp` with a lazy GPU backend. The subsequent source migration also relocates the remaining MATLAB tree, Simulink adapters and circuit/MAT-file exporters to their integration owner, and independent Optuna MATLAB sources to `packages/radia-optuna`. Reusable thermal assembly and state-space mathematics remain CPU domain code. Path repairs, ownership classifications, final artifact inspection and installed acceptance remain pending; this is not a completed distribution-boundary acceptance.

Transfer commit `54a577e7f` passed 397 Python/package tests, 16 focused integration contracts, installed CUDA/Simulink/LTspice acceptance and real MCP tools/list plus topology comparison. Its Poisson CPU parity difference was 3.52e-10 and GPU true relative residual 9.99e-11. Public migration evidence is in `validation_test/radia_simulink/results/transfer-54a577e7f.json`; raw operational evidence and the exact wheel are retained privately. These checks certify this transfer stage, not the remaining source migration or publication.
