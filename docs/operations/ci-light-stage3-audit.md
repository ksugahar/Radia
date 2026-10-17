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

## Integration with the strengthened pre-push checks

The stage-3 branch merges main's pre-push change rather than replacing its
lower CI policy test or impact rules. The additional eligible tests load
`tools/pre_push_check.py` and `tools/local_push_checks.py` through explicit
file loaders. Both helpers use standard-library source/Git inspection; their
tests mock remote dispatch and subprocess execution, so no SSH, live hook,
package build or numerical runtime is invoked by this lane. Both helper files
are now part of the reviewed execution closure. The latest worst-case selector
contains 447 cases, including the added pre-push contracts; the earlier 424-case
record describes the pre-integration selector, not a skipped subset.

## Current-main conflict resolution (2026-10-17)

The integration retains main's CPU/GPU ownership, isolated GPU source and
HDF5 migration, and preserves stage 3's capability-isolated lightweight route.
Only the two policy documents conflicted. Their shared body remains identical
and below the policy size limit. No runner registration or runtime is changed.

A fresh isolated environment containing exactly the declared lightweight
requirements passed the worst-case fast selector: 447 passed, zero skipped.
The routing, policy, workflow-owner and release-ref contracts passed 68 tests.
All 45 workflow jobs still pass the stage-3 routing audit. The reviewed
execution-closure digest remains unchanged; no newly merged solver or HDF5
fixture test was silently added to the lightweight execution closure.

The follow-up merge includes main's A100 routing policy (#440). GPU host
names in the CI paragraph and its contract refer to A100, distinct from the
NGSolve host. Updated policy/owner/ref contracts passed 49 tests.
The 45-job execution-closure audit is unchanged, so the completed 447-case
lightweight execution remains applicable without repeating solver-free tests.


## Follow-up validation-source gate closure (2026-10-17)

The validation-source reuse candidate adds 39 deterministic metadata/trust
contracts to fast-contracts. Its dynamically loaded
`tools/validation_source_gate.py` is explicitly part of the audited closure.
A fresh Python 3.12.10 venv containing only pytest, PyYAML, NumPy, SciPy,
setuptools and wheel executed **486 tests, zero skips**, in 13.50 seconds
(14.28 seconds through the tier runner). No NGSolve, MATLAB Engine, Cubit,
CUDA toolkit, h5py or MCP dependency was installed in that audit environment.
The 51-file execution closure digest is
`bd004bc5dd2dcd23a3a1a6bf3e679eb17e6d25a681b071a419bbdb1f12119e5b`.
This audit does not certify a live validation cache hit, native artifacts or
physics results; those remain separate acceptance gates.
