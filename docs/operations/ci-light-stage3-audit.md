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

## Validation-reuse review amendments (2026-10-17)

The reviewed changes inspect untracked inputs repository-wide, bind compiler
flags and numerical environment variables, and exercise the actual receipt
search. The validation job now fetches complete main ancestry. In the same
isolated lightweight environment, the worst-case fast selector executed
**505 tests, zero skips**, in 29.84 seconds (30.98 seconds through the runner).
Only the declared lightweight dependencies listed above were installed; no
solver, MATLAB Engine or GPU compilation was used. The 51-file closure digest
after these execution changes is
`8dee341aab927e5fc050b1025d53fa98ee7f3b44ac59fe9a9bdf8c363431bc2f`.
All 45 workflow jobs pass the routing audit. This repeats the execution audit,
rather than updating a stale digest alone. Live main receipt hit/miss and
current native artifact checks remain separate required evidence.

## Bound native build metadata and job-context audit (2026-10-17)

The follow-up binds four exact-name native modules and copied build outputs
to current clean-source sidecars with a declared build directory. Only known
CMake/Ninja metadata there and two reviewed runner-temp context paths are
normalized; unknown native and runtime inputs still force execution.
The isolated lightweight worst-case selector passed **515 tests, zero skips**,
in 23.71 seconds (24.58 seconds through the runner). Its 51-file digest is
`07d599ee89666da020d2d6fefc4583358ac0ef101245a41357537c10d9fa0f0e`. No solver or MATLAB dependency was installed for this audit.
Real native-build lanes and main hit/miss evidence remain required separately.

## Installed runtime-copy and rejection diagnostics audit (2026-10-17)

The explicit MKL runtime-copy check binds its installed package version and
DLL hash. Rejected untracked paths are now retained in verification JSON, and
compiler-output suffix exemptions are limited to the declared build directory.
The isolated lightweight selector passed **520 tests, zero skips**, in 23.77
seconds (24.72 seconds through the runner). The 51-file digest is
`eae589c6e263d297e8533035bbde26e0950a67ec0093a0dca10c52447f4b15fb`. A read-only observation of the actual native CI build found
no persistent extra per-target response files or manifests requiring a broader
exemption; its complete observed file list is retained as private operations
evidence. Real main receipt hit/miss acceptance remains a separate gate.

## Distribution rename and admission changes (2026-10-20)

The lightweight closure changed with the distribution rename (wheel and
TaskManager audit tooling paths) and the explicit build-admission tests. In the
isolated lightweight environment (only the declared dependencies), the
fast-contracts job steps passed: source compilation, the caller-owned
TaskManager audit, the routing contract (19 tests), the version check, and the
worst-case fast selector (137 tests, zero skips, 15.3 s). No solver, MATLAB
Engine or GPU compilation was used. The 51-file closure digest is now
`e0f823ec0563ecdfaf3305abca6b6994c206d8585fe1968f1ce037d6a9d38e23`.
All 45 workflow jobs pass the routing audit.


## HDiv BH contract and Newton budget repair (2026-10-21)

The changed `tests/test_hdiv_bh_table_saturation_contract.py` and its dynamic
source `validation_test/feec/bh_saturation_audit.py` were re-audited in a fresh
isolated Python 3.12.10 venv containing only pytest, PyYAML, NumPy, SciPy,
setuptools, wheel and their lightweight dependencies. The curve test executes
selected scalar production helpers and SciPy interpolation; it does not import
NGSolve or a native solver. No MATLAB Engine, Cubit, GPU or solver dependency
was installed in this environment.

All F2 steps passed: routing (45 jobs), routing mutations (19 tests), source
compilation, caller-owned TaskManager audit, CPU version (1 test), the actual
fast selector with `--since origin/main` (342 tests, zero skips, 24.31 seconds;
25.28 seconds through the tier runner), and showcase readers (6 tests).
The Newton inner-solve failure contract also passed 16 tests, including the
explicit compensated-budget forwarding regression. This is source-contract
coverage, not native numerical or release acceptance.

The 51-file execution closure digest is
`519820bd4a6447418bdbd59ce3fec8ec393ad63d76a4cbad56fa65ef520ee345`.
The single audit used two physical cores on one NUMA node, one library thread,
BelowNormal priority, an 8 GiB job-commit cap and a 240-second wall cap, with no
GPU and no retry. It completed in 81.14 seconds with a peak job commit of
132526080 bytes. Exact topology, affinity, dependency versions and raw outputs
are retained privately; no runner or shared environment was changed.
Re-audit after the execution-role change. The closure files changed were
tests/test_public_runtime_privacy.py, tests/test_radia_promotion_gate.py,
tools/ci_distribution_manifest.json, tools/policy_lint.py, tools/pre_push_check.py
and tools/release_acceptance.py. In the isolated lightweight environment, the
fast-contracts job steps passed: routing audit, routing contract (19 tests),
source compilation, caller-owned TaskManager audit, version check, fast
selector (135 tests) and showcase readers (6 tests). The 51-file closure digest
is now `c59f7d870b417d07d7de003cdd27bdc1c3bec3e53a06646203029acbcf437e9e`.

## License metadata contract (2026-10-21)

Added source-only SPDX expression, packaging metadata and external DLL boundary
checks. Policy lint, repository contracts (481 tests), caller-owned TaskManager
audit and fast contracts (142 tests) passed without native builds. The metadata
check uses setuptools configuration and metadata/file selectors; it does not
execute package setup scripts or compile solvers. Updated the light closure
digest after these checks: `9bfbe0f4f595e581cc5d3ae297ad0a21ebd0cf46316def46dad1aeed22e335cf`.

After registering the new root license files (LICENSE-APACHE-2.0, NOTICE,
THIRD_PARTY_NOTICES.md) in the distribution manifest, the isolated fast-contracts
steps passed again (routing contract 19 tests, TaskManager audit, version check,
fast selector 142 tests, showcase readers 6 tests). The 52-file closure digest is
now `2d6f0f35bf7de375a25bd0e7dd86fe24c26180af43cbb4c0e84403f74938e94e`.
After merging the Newton budget fix with current main, the combined 51-file
closure digest is `0d933121845010f273cd88353a810aa7253f472f539f56783a19471ae4e7f26d`; the isolated fast-contracts steps were re-run on the
merged tree.

After merging current main (resident AMG build switch, Apache-2.0 relicensing)
into the HDiv runtime candidate with the setuptools-free distribution contract,
the shared precheck passed (light fast-contracts steps in an isolated venv,
policy checks, repository contracts). The 52-file closure digest is now
`2b7175ea2d8a8b081d2419d4bd549de6a8d04766abd285f78919b1f5f989df2d`.

## Heavy CI throughput routing (2026-10-21)

Reviewed the lightweight execution closure after introducing identical-tree
receipt verification and full-closure AMS PR input selection. Both helpers use only
standard-library modules. Network verification executes in hosted routing jobs;
source-only contracts use injected API fixtures. The AMS selector parses files
without importing solver modules, loading native libraries or starting builds.
The closure now contains 54 files; its digest is
`03eaeb5881ac7cd81d7e400c811218236a1a9c924c845c77154b39c18d92aea0`.
The isolated ci_precheck verifies routing, routing mutations, compilation,
TaskManager ownership, version consistency, fast contracts and showcase readers,
followed by policy and repository contracts. Runner count and heavy labels are
unchanged. See [heavy CI execution policy](ci-throughput.md) for release gates
and the full numerical audit schedule.

After merging current main (CI throughput routing) into the NGSolve 6.2.2608
migration, the shared precheck passed on the merged tree. The 54-file closure
digest is now `12f26b1613e69e623b4314a92af07c7d2ad10e14a91a8c8b2f2bf9148f7e8c71`.

## GPU capability routing (2026-10-21)

Re-audited the 54-file light execution closure after adding `gpu-rtx` and
`gpu-a100` routing and the compute-only distribution ownership declarations.
The routing mutations check missing GPU capabilities and reject PR-triggered
A100 validation or compiler commands. No MATLAB, CUDA, native solver or network
execution enters light contracts; workflow commands are inspected as strings.
The artifact helper is standard-library code tested separately using local
wheel fixtures and an injected empty Actions response. Light requirements and
execution steps remain unchanged. The reviewed closure digest is
`2b4544a4d0a685032d07845b89e869e378c84944ad8b1bcf7d63fcd7c755fef2`.

Re-audited after the GPU nightly amendment and rebase: the routing contracts
inspect selected-SHA workflow strings, separate dispatch concurrency and
failure-only reporting without adding executable GPU or network dependencies.
The closure remains source-only; the current digest is
`7c7b29639ad0c648e648be3dec79ef43a02445054a65f3fcdbe305da3311edfa`.

## Linux source-build route contracts

The Linux CPU build entry point is included explicitly in
`audited_dynamic_sources` because `test_linux_build.py` loads it dynamically.
The new impact contract loads definitions and package metadata; pipeline tests
mock compiler lookup and subprocess execution. The cross-platform provenance
check uses only temporary Git repositories and the existing PowerShell source
identity helper. No native build, NGSolve import, solver, scheduler submission
or remote access executes in this lightweight test.

On the Windows Python 3.12 source-check environment, the eight new tests passed
in 3.43 seconds; the affected fast-contract selection passed 192 tests in
19.18 seconds. These are observed runtimes, not acceptance thresholds. The
executable closure digest was recalculated with `audit_digest` and the declared
dynamic source included. Runner labels, denied imports, dependencies, job steps
and the active profile budget remain unchanged. Reciprocal CI review still
applies; source contracts do not qualify Linux or Windows native artifacts.

## Runner-owned scratch (2026-10-23)

The fast-contracts job now initializes its output and standard temporary
variables from RUNNER_TEMP in a shell step, before checkout and child processes.
No fixed host output/history path remains in the workflows. Artifacts preserve
the result files before runner cleanup; this changes storage, not gate criteria.
Runner context is resolved at step execution, not in unsupported job-level env.

After retaining the Linux source-build route above, the F2 executable source
closure contains 56 files. Relative to the integrated main, only the light job
execution declaration changes; requirements, denied imports, runner capability,
profile budget and source-only import guard remain unchanged. Recomputing with
the prior job declaration reproduces the prior seal. The new audit_digest is
`5ff1ee6e101a637220e585e2de887bf925839516c66710348d08e86d0e61f6d2`.
Projection allocation checks belong to hosted policy contracts and do not add
MATLAB, a GPU, a compiler or a native build to the light execution closure.


## Historical trace provenance (2026-10-24)

The trace-reuse contract now verifies the historical Git blob and marks current
runner applicability stale; it uses standard-library Git inspection, preserves
the historical numerical gates, and adds no solver or GPU imports. Recomputing
with the prior executable source and unchanged declaration reproduces the old
seal; the reviewed closure digest is `df15b00bfccb0056284d27f6a5259e8401102aa3ec5542299cca272f4c4bf4a4`.
The shared precheck repeats F2 in its isolated lightweight environment.
## Numerical plan consumers (2026-10-24)

The tier runner can consume an exact-source numerical plan. Its shared verifier
and planner are standard-library-only; they inspect immutable Git blobs and do
not import production numerical modules. The sealed declaration JSON is included
explicitly because dynamic loaders are reviewed inputs, not inferred exemptions.
The source contracts exercise real CLI processes and temporary Git repositories,
without solver execution, native compilation, MATLAB or GPU use. Source-only
import rejection, runner labels, requirements and the active light budget remain
unchanged. Recompute the executable-closure-v2 seal with `audit_digest` after
integrating these sources; source contracts do not qualify native artifacts.


## Independent equation-product retirement

The equation product's workflows and routing entries are removed; its installed
package consumer remains separate. The publication-policy wrapper now loads
reviewed external criteria through a generic checker, explicitly included in
`audited_dynamic_sources`. No native build, solver or Office execution is added.
The prior seal was reproduced from the checked base blobs before computing the
new executable-closure seal: `373305577e2d0f2d77743741f09bd2b45d1d49e416fecfccf33f9af83a96f5e2`.
Runner labels, requirements, source-only import rejection and the active profile
budget are unchanged. The shared precheck executes the routing mutations and
light source checks; reciprocal review confirms the revised closure.
## Explicit BLAS provider build contracts

The NVPL source-build route adds source-only admission tests for architecture,
provider mismatch, required development files, LP64 library hashes, and shared
provider-header ownership/event routing. They use temporary byte fixtures and
read the build helper without invoking its native build entry. No library is
loaded and no compiler, numerical solve or GPU call is made. Existing Linux
pipeline tests continue to mock toolchain and package metadata; their actual
subprocess work is temporary Git and the PowerShell provenance helper.
The six new provider contracts ran in 0.59 s; the combined provider/Linux
contracts initially ran 13 tests in 4.28 s. These are source checks only.
The executable seal is regenerated through ci_runner_routing.audit_digest;
independent F2 review and native provider acceptance remain separate gates.

## HDiv Python ownership build preparation

The HDiv ownership change adds two impact-selected source contracts. The
build contract loads `tools/build_python_gpu.py` but replaces command execution,
SDK metadata and compiler inputs with temporary fixtures. It never calls the
builder CLI or a compiler. The workflow contract loads `radia/hdiv_gpu.py`,
replaces native construction with test doubles, and exercises lifecycle,
explicit backend rejection and immutable mesh snapshots. The real label-checker integration cases live separately in
`solver-numerics`: package initialization loads Netgen even though the checker
does not launch Cubit. The light workflow contract mocks the subprocess; it does
not execute the checker or import Netgen.
NGSolve and Radia VIM imports in the workflow are deferred into FE/material
methods that these contracts do not enter. The source-only native-import guard,
capability routing, denied requirements and elapsed-time budget remain intact.
The repeated F2 check executes the worst-case fast tier and the routing mutation
suite; native compilation, numerical GPU execution and performance are separate
unverified gates.

## Python GPU distribution and portable checked-mesh receipts

The previous executable-closure seal was reproduced from the merged base
before inspecting this extension. The explicit dynamic closure now includes
the Python GPU payload verifier/stager and the build-role receipt producer.
The payload contracts use synthetic non-executable binary bytes and replace
compiler execution, SDK metadata and native construction. The producer imports
only standard-library modules; its real checker subprocess is not entered by
the light tier. Receipt workflow tests mock FE/native entry and package metadata.
No compiler, checker, NGSolve, GPU or MATLAB operation is added to fast-contracts.

The base fast tier executed 200 tests successfully in 8.18 seconds
(tier elapsed 9.42 seconds). All 45 routing mutation tests passed in 111.95
seconds. The first attempt used an incomplete test environment lacking
setuptools and wheel, yielding eight failures; the corrected run used the
existing complete isolated test environment without source or oracle changes.
The denied-import set, source-only guard, runner labels, requirements and active
profile budget are unchanged. The resulting executable closure seal is
`f276d0f3d1651c0ab1d9f9c9ba094e3173e924f31aaf409e3d49114dc4d57c77`.

Separately, a real installed checker produced plain and compressed input
receipts. A fresh consumer environment without that checker accepted both
and rejected 18 content/provenance/options/report mutations. This CPU transport
check does not certify GPU runtime, numerical results or the eventual wheel.

The subsequent Windows compile check exposed SDK macro re-evaluation of
backslash-containing CMake values. The builder now renders each CMake path
argument with forward slashes while preserving separate subprocess argv.
The corresponding source contract uses mocked SDK files and command execution;
it guards compiler, SDK, Python, object, output and CUDA-library path values.
Its closure has no additional paths or denied-import exemptions. Reproducing
the preceding seal and inspecting these two changed executable files yields
`26a991e809caa0a0f96a00d9cd8035648c7517835a00beef93711a5699a52e39`.
Fresh exact-source fast/routing checks and native compilation are recorded
separately; the previous partial native build does not qualify this revision.

The base-profile run above is not the all-impact worst case. The latter must
invoke `tools/run_test_tier.py --profile fast-contracts --since ""`; omitting
`--since` selects only the base profile. The first such all-impact execution
ran 639 tests and exposed one stale Windows path assertion in a mocked Linux
SDK contract. The assertion now compares the same canonical forward-slash
CMake representation as the builder, without changing SDK identity, compiler,
libraries, numerical criteria, denied imports or the profile budget. Fresh
source-bound execution reports distinguish these runs and the unchanged-code
compile evidence from any later exact-release native artifacts.
The final test-path-only closure seal is
`6a4402eac907a661e0222fbd8e8c15b046319def31fb1576162c58fc8da31bac`.

## IH fin physical verdict contract

The tier manifest maps `validation_test/induction_heating/fin_return_principle.py`
and its contract `tests/test_ih_fin_physical_verdict.py` to that contract. The
contract loads the checker by file path, so the checker is declared in
`audited_dynamic_sources`. Its NGSolve, Netgen, h5py and Radia solver imports
are function-local to the solve and mesh routes; the contract enters only
`physical_verdict` and `PHYSICAL_TARGETS` (standard library and NumPy), so
`ngsolve` and `netgen.occ` are declared deferred for that file. The prior seal
`43d9572de6d61b665ed543e605c7865cc41abf71139b9b388cf20b99300bbe35` was
reproduced on the merged base before the change; the closure grows from 74 to
76 files (the contract and the checker). Runner labels, requirements, denied
imports, the source-only guard and the profile budget are unchanged. The
`audit_digest` seal is
`8426c0d457c73593aca30ae3c66003dc1e61edcbaf02c8f79e0355aafdd2873f`.
