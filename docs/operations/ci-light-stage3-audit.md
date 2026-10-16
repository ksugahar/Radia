# Stage 3 lightweight runner audit

The migrated `Radia / fast-contracts` job selects exactly
`[self-hosted, Windows, X64, ci-light]`. It intentionally omits `100` while
legacy local jobs still select generic `100`. Check names are unchanged.
No runner registration, account, service, label or scheduling action is included.

## Routing contract

`tools/ci_runner_routing.json` classifies all 45 jobs in the current workflow
files by owner and target. Each existing generic local job is an explicit
`legacy_100` exception, not a default for new jobs. New, missing, renamed,
expression-selected or unclassified jobs fail closed. `--mode final` rejects
all legacy routes and requires `100` plus exactly one declared capability on
every local job. GitHub-hosted jobs remain separately classified.

The lightweight audit examines actual run commands, actions, declared
requirements, the real worst-case impact selector, local Python imports,
intermediate pytest conftests, and explicitly reviewed dynamic-loader sources.
Its SHA-256 binds executable source, test selector, dynamic sources and job
steps. An executable change needs a repeated F2 review and deliberate digest
update; it is not automatically approved by regenerating the digest.
This is a review contract, not a security sandbox or complete Python call-graph
proof. Dynamic Python cannot generally be proven safe by syntax inspection.

## F2: actual execution

| Step | Actual requirement / operation |
|---|---|
| Checkout | Git; checkout action only |
| Isolated environment | Python venv; pytest, PyYAML, NumPy, SciPy, setuptools, wheel |
| Routing validation and mutations | Read workflows; parse Python; copy temporary source; small Python subprocess |
| compileall | Python bytecode compilation; no native compiler |
| TaskManager audit | AST/source inspection, not TaskManager execution |
| Version | Read Python and TOML; no Radia package import |
| Fast tier | Nine compact contracts plus eligible source-only impact contracts |
| Showcase | Read saved notebook/result records; no solver execution |
| Ref context and artifacts | PowerShell JSON/file IO and upload action |

Two inherited hazards were found and corrected:

1. The parent pytest conftest previously probed NGSolve, Netgen and Cubit even
   in the source lane. `RADIA_SOURCE_ONLY_CONTRACTS=1` now avoids those probes
   and raises on subsequent native/GPU imports. Required tests cannot disappear
   as dependency skips. Non-source contexts retain their original behavior.
2. `test_lamination_api.py` contains real native public-import and MCP-adapter
   parity checks. The whole file now belongs to `solver-numerics`, which the
   existing native workflow runs after building the package facade. Source and
   test paths were added to both native workflow triggers and their checked
   impact projections. Tests prove the actual router selects this lane and
   that removing those paths loses the selection.

The dynamic surface-ROM contract deliberately supplies NGSolve and model test
modules; it executes only the residual-check helper. The heat audit imports
its module but does not enter deferred NGSolve solver functions. These two
specific imports are declared and bound by the audit digest; other MATLAB,
CUDA or native imports are not exempted. Synthetic batch guard tests execute
Python exit codes, not CMake/MSVC builds. Package tests construct dummy wheel
archives; they do not load their dummy native payloads.

The initial worst-case run exposed the lamination import: 447 passed, one
native-dependent failure and one optional MCP-adapter skip. After moving that
file to its required native lane, the same worst-case selector completed
424 tests with zero skips. Routing mutations cover wrong/multiple capabilities,
unknown jobs, direct MATLAB/CUDA/native commands, hidden entrypoint imports and
native compiler subprocess changes, declared heavy dependencies, final mode,
conftest probe prevention, and retained native coverage. No native solver/build
or GPU run is needed for this audit.

## Later migration / policy delta

The MATLAB capability runner must stay in an authenticated interactive account,
started by `run.cmd`; do not service-host it. Verify Engine startup, calculation
and shutdown in that actual runner context. Update the plan sections 3 and 5
(stage 4) accordingly. The stage-3 `100` omission is a temporary safety
exception; restoring the label belongs to the management executor only after
all local jobs use capabilities. GPU exclusivity, memory/CPU admission and
uncontended timing windows remain later operational gates, not claims made by
this source-only PR.
