# Three-runner CI plan (review draft, 2026-10-16)

Status: approved staged plan; no blanket service authorization. Baseline: `ebb117b80` (includes
the CPU workflow rename in PR #427). Codex implements and validates; Claude
repository management reviews and owns GitHub integration. The user reviews
this plan before runner registration or service configuration. Stage 3 is now
approved and the isolated light runner is registered. Tokens remain
user-owned and are never placed in source, logs, manifests or command history.

## 1. Routing and contracts

Migration exception: light must NOT carry `100` while heavy jobs still select
the generic label. Stage 3 uses `[self-hosted, Windows, X64, ci-light]`.
Restore `100` only after all heavy jobs select explicit capabilities. Contract
tests validate the transitional and final states separately.

Three runner instances share the Windows build host; A100 is a separate
compute role. CPU heavy jobs currently use `100`, RTX jobs add `gpu-rtx`,
and light jobs use `ci-light` without `100`. The final CPU capability is
`ci-native`; A100 uses `gpu-a100` without `100`:

| Capability | Purpose | Proposed job budget: CPU / private memory / wall time |
| --- | --- | --- |
| `ci-native` | Native compiler, CPU numerical validation and MATLAB CPU contracts | 16 threads / 64 GiB / 60 min |
| `gpu-rtx` | RTX GPU wheel build and acceptance | 8 host threads / 32 GiB / 60 min |
| `gpu-a100` | Artifact-only GPU compute validation | 8 host threads / 32 GiB / 90 min |
| `ci-light` | Source/policy/compact contract checks only | 2 threads / 8 GiB / 20 min |

These are initial proposed ceilings, not measured adequate settings. Admission
of each existing job requires its observed peak and duration to fit; otherwise
review a specific exception before migration. Do not reduce numerical gates to
fit. The current host reports 104 logical processors and about 190.7 GiB RAM;
this is not 104 independent physical cores. The proposed total is 26 threads
and 104 GiB, retaining substantial capacity for users and the OS. Keep a global
64 GiB available-memory floor. Re-measure hardware at installation; private
hardware identity and concrete paths remain in the private operator record.

| Current workflow / job | Target | Per-job override |
| --- | --- | --- |
| `build-test.yml / build-test` | native | 16 threads, 64 GiB, propose 90 min; measure before adoption |
| `sparsesolv.yml / ams-regression` | native | 16 threads, 64 GiB, existing 45 min |
| `radia-gpu.yml / cuda-mex` | MATLAB/GPU | 8 threads, 32 GiB, existing 60 min |
| `radia-optuna.yml / build-test` | MATLAB/GPU | 8 threads, 32 GiB, propose 60 min |
| `radia-optuna.yml / installed-wheel-matlab-e2e` | MATLAB/GPU | 8 threads, 32 GiB, existing 45 min |
| `ltspice-data-safety.yml / matlab-data-contracts` | MATLAB/GPU | 8 threads, 32 GiB, existing 5 min |
| `radia.yml / fast-contracts` | light | 2 threads, 8 GiB, existing 20 min |

Optuna `build-test` actually calls mexcuda: its generic name does not make it
lightweight. Existing GitHub-hosted jobs remain GitHub-hosted, including impact
selectors, policy, MCP minimal/selftests, Cubit native packaging,
Windows/Web jobs and publication jobs. They do not need another local runner.
See the complete baseline job inventory below; unknown/new jobs require a
classification instead of inheriting the legacy generic label.

Implementation stage: add a checked routing manifest and contract tests that
parse every workflow job, reusable workflow and local action entry point.
Assert that each local job has exactly one capability, all seven listed jobs
have their declared owner, and no workflow selects only generic `100`. Mark
MATLAB/Engine/CUDA/Simulink requirements explicitly per job; test those
requirements against commands and called wrappers. Negative mutations must
reject Optuna mexcuda on light, Engine on light, missing classification, and
conflicting capability labels. Light jobs must not install MATLAB/CUDA or
invoke native compilation. Do not equate a simple word scan with transitive
dependency verification; unknown script entry points fail closed for review.

## 2. Concurrency and resource enforcement

Each runner executes one job, but all three can execute simultaneously. A
single dedicated MATLAB/GPU instance serializes GPU CI. Every GPU workflow
must select it. Do not start additional generic runners capable of accepting
GPU work. Cross-workflow, manual MCP and specialist use share one host-local
OS GPU lock, with owner/run/PID and bounded acquisition; an abandoned OS lock
releases on process exit. Test two competing owners and process termination.
Read-only telemetry is allowed without GPU ownership; execution is not.

Do not rely on a GitHub concurrency group as a FIFO queue: default groups can
replace pending work. Avoid dropping required checks. The dedicated single
runner plus OS lock is the primary exclusion mechanism; retain existing
workflow concurrency only after confirming its cancellation semantics.
See [GitHub concurrency documentation](https://docs.github.com/en/actions/concepts/workflows-and-actions/concurrency).

CPU environment limits alone are advisory. Future job launcher must assign a
Windows Job Object before resuming its owned process tree, enforce aggregate
private committed-memory limits, restrict CPU affinity to each role's allocated
logical processors and set tool thread limits (`/m`, CMake, BLAS/OpenMP and
NGSolve under its existing region policy). MATLAB uses supported MATLAB
scheduling; mexcuda never starts NGSolve TaskManager. Select physical-core
topology-aware CPU sets during setup; do not guess affinity from core numbers.
Record actual limits and descendants. Per-job timeout, host memory floor and
GPU reserve `max(2 GiB, 20% total device memory)` stop only owned work. GPU
reserve is a precheck/monitored threshold, not a hard VRAM allocation guarantee.

On a limit or lock failure, fail with JSON diagnostics; never kill another
runner, user's MATLAB/CAD, or silently switch backend. Native and GPU jobs may
be co-admitted only within aggregate caps. An oversized native job requests an
exclusive window rather than stealing the light runner's capacity.

## 3. MATLAB authentication in the real runner context

User decision: keep MATLAB/GPU in the authenticated account's interactive
session started with `run.cmd`; do not service-install this runner. Verify
license state in that actual runner context, not SYSTEM, SSH or a separate
manual desktop test. Do not copy credentials between accounts. Authentication
failure blocks this lane; any later service migration needs separate review.

In a diagnostic job routed to the new MATLAB runner, record account identity
privately, runner role, Python executable, MATLAB/Engine/CUDA/compiler
versions, then start a new official Engine, calculate a deterministic matrix
product, verify its value, query the CUDA device and execute a tiny GPU
calculation, and quit the owned Engine in `finally`. Verify its PID ended.
Repeat in that actual interactive runner context after approved reboot/relaunch;
logout and remote disconnection must not be assumed equivalent.
GUI/Simulink acceptance stays separate where required. Use the MATLAB release's
supported Python Engine package:
[installation](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html),
[startup and shutdown](https://www.mathworks.com/help/matlab/matlab_external/start-the-matlab-engine-for-python.html).
Authentication failure is a red gate, never an automatic skip.

## 4. Isolation of working folders and environments

Give each instance a separate installation directory, `_work`, `_temp`, job
environment and output/cache namespace. Concrete absolute paths and account
mapping are private operator configuration. Never copy a configured runner
directory with its registration identity. Caches are role/toolchain keyed and
immutable on restore; do not share writable venvs, MEX build directories,
source checkouts or MATLAB preference/temp folders across runners.

Run artifacts include role/run/job/SHA and hashes; keep production/test-only
MEX receipts isolated. No system-Python installs, mapped-drive dependency or
running editable runtime mutation. Recover and hash-check failed and successful
diagnostics before exact job-owned scratch cleanup; verify absence. A rejected
cleanup is recorded with owner/reason/trigger, not retried through another tool.

## 5. Staged migration and recovery

1. Review this plan and proposed policy text; agree on budgets, account and
   folder layout. Read-only inventory is allowed now; no registration/services.
2. After user approval, Claude management is sole configuration executor.
   User handles short-lived registration tokens directly; no token files or
   echoed commands. Install distinct runner copies and names/identities; grant
   each only its own capability. Keep the current runner until migration gates.
3. Create light runner first, validate source-only routing and resource guard
   negative tests, then change only `fast-contracts` in a reviewed PR. Keep
   required check names stable. Confirm a real required PR check finishes.
4. Create/validate MATLAB runner account and Engine/GPU lifecycle in the actual
   interactive `run.cmd` context (no service), then migrate all four MATLAB jobs together. Require the GPU OS lock
   for manual work too. Validate a small published numerical contract.
5. Convert or replace the old runner as native-only after its current job drains.
   Migrate native jobs, retire generic eligibility, and test all three lanes
   simultaneously with bounded workloads. Verify no shared-path writes.
6. Soak over representative PRs and one approved restart. Compare queue latency,
   failure recovery, hardware pressure and gate completeness; store JSON.

Recovery: pause admission of the affected new lane, let owned running jobs end
or follow their declared stop rules, preserve logs, and revert only the reviewed
`runs-on`/launcher/policy migration commit to the known previous configuration.
Temporarily use the old runner only after confirming it is idle and MATLAB
authentication still holds. Do not roll back published artifacts or other WIP.
Keep instance directories for diagnosis with owner/cleanup trigger; removal and
unregistration require the applicable explicit authorization. Never leave
required checks pending silently on a label with no online eligible runner.

## 6. Performance measurement windows

Adopt an explicit benchmark reservation. Drain and stop admission to all three
local runners, then stop their idle listeners for the approved window
(configuration executor only; never stop an active job). Also check manual
jobs, free RAM and GPU, and prohibit competing specialist work. A paused
workflow alone does not stop an already queued job reaching another runner.
Record start/end listener states, process inventory and CI status.

Only benchmark after all local runner listeners are stopped and exclusive
ownership is confirmed. Fix CPU threads, include transfers/setup and cold/warm
runs, and save hardware/software/hashes in JSON. Unexpected CI/manual overlap
invalidates affected timing. Restart the same listeners afterward and verify
online state and a harmless job per lane. Accuracy checks may run concurrently
under caps, but are never relabelled benchmark timings. Existing contended
measurements remain diagnostics, not accepted performance evidence.

## 7. Policy alignment and push-before-CI checks

Preserve the `100` build-host label, independent distribution boundaries, wheel
tests on the wheel host, MATLAB-account proof and numerical release gates.
Existing hosted workflows remain explicit exceptions. Proposed focused text
for BOTH AGENTS.md and CLAUDE.md (not edited in this draft):

> Build-host CI uses ci-native, gpu-rtx and ci-light capabilities;
> ci-light does not require the 100 label. Artifact-only compute validation uses
> gpu-a100 on the compute host. GPU CI and manual GPU jobs share one
> exclusive host lock. Each job enforces declared CPU, aggregate process-tree
> memory and timeout budgets. Existing GitHub-hosted package/source jobs remain
> explicitly classified. Performance acceptance runs only in reserved windows
> after all local runner listeners drain and stop; restore and verify afterward.

Do not apply this text before review; existing root policy edits belong to
other work and remain untouched. No new runner on the wheel host or retired
evaluation targets. No heavy solver/MATLAB execution in the pre-push gate.

Parallel pre-push implementation candidate expands the existing exact-SHA,
impact-scoped wheel-host preflight from main-only to every branch push,
including PR branches and new branches. Deletes and tags stay separate;
release-tag checks remain release-workflow owned. Validate all stdin records
before any remote work, stop a multi-ref push on first red result, fail closed
on transport/missing helper, and retain an explicit visible emergency bypass.
Use checkout-scoped Git trust rather than wildcard trust. Do not install the
updated shared hook or run SSH preflight before this candidate is reviewed.

Pre-push tests cover main/feature/new branches, exact outgoing SHAs (not HEAD),
multi-ref failure, malformed input and unavailable helper/transport. Extend
CI contract selection for these tests when integrating; three-runner routing
contract implementation belongs to the separately approved migration. Local
preflight is early feedback, not a replacement for GitHub CI/release gates.

Acceptance before registration: Claude review disposition and user's plan
confirmation. Acceptance before routing adoption: actual runner identities,
account/Engine evidence, contract mutations, bounded simultaneous-load test,
recovery test and staged PR CI. Unresolved choices: final account/context proof,
measured per-job exceptions, physical CPU sets and enforceable benchmark
reservation procedure. No implementation may silently assume these resolved.

## Baseline complete job inventory

All 45 workflow jobs at the baseline are classified below. The table is a review snapshot, not executable routing configuration.

| Workflow | Job | Proposed execution |
| --- | --- | --- |
| `build-test.yml` | `build-test` | native |
| `cubit-mesh-export.yml` | `ci-impact` | GitHub-hosted (retain) |
| `cubit-mesh-export.yml` | `build-test` | GitHub-hosted (retain) |
| `cubit-mesh-export.yml` | `publish` | GitHub-hosted (retain) |
| `ltspice-data-safety.yml` | `ci-impact` | GitHub-hosted (retain) |
| `ltspice-data-safety.yml` | `matlab-data-contracts` | MATLAB/GPU |
| `policy-lint.yml` | `policy-checks` | GitHub-hosted (retain) |
| `policy-lint.yml` | `repository-contracts` | GitHub-hosted (retain) |
| `radia-gpu-release.yml` | `qualify` | GitHub-hosted (retain) |
| `radia-gpu-release.yml` | `publish` | GitHub-hosted (retain) |
| `radia-gpu-release.yml` | `github-release` | GitHub-hosted (retain) |
| `radia-gpu.yml` | `ci-impact` | GitHub-hosted (retain) |
| `radia-gpu.yml` | `cuda-mex` | MATLAB/GPU |
| `radia-mcp-matrix.yml` | `ci-impact` | GitHub-hosted (retain) |
| `radia-mcp-matrix.yml` | `scope` | GitHub-hosted (retain) |
| `radia-mcp-matrix.yml` | `minimum-sdk` | GitHub-hosted (retain) |
| `radia-mcp-matrix.yml` | `selftest` | GitHub-hosted (retain) |
| `radia-mcp-matrix.yml` | `release-evidence` | GitHub-hosted (retain) |
| `radia-mcp-matrix.yml` | `publish` | GitHub-hosted (retain) |
| `radia-optuna.yml` | `ci-impact` | GitHub-hosted (retain) |
| `radia-optuna.yml` | `build-test` | MATLAB/GPU |
| `radia-optuna.yml` | `installed-wheel-matlab-e2e` | MATLAB/GPU |
| `radia.yml` | `ci-impact` | GitHub-hosted (retain) |
| `radia.yml` | `fast-contracts` | light |
| `release-radia-optuna.yml` | `publish` | GitHub-hosted (retain) |
| `release.yml` | `upload-binaries` | GitHub-hosted (retain) |
| `release.yml` | `publication-hold` | GitHub-hosted (retain) |
| `release.yml` | `verify-promotion` | GitHub-hosted (retain) |
| `release.yml` | `publish-pypi` | GitHub-hosted (retain) |
| `sparsesolv.yml` | `ci-impact` | GitHub-hosted (retain) |
| `sparsesolv.yml` | `ams-regression` | native |

Eqnedit64 source and hosted release workflows are owned by the independent
`ksugahar/eqnedit64` repository, outside Radia runner ownership.
