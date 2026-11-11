# CI runner roles

The CI pool has three roles. Capability labels select hardware; they are not distribution names.

| Role | Labels | Work |
| --- | --- | --- |
| CPU heavy | `self-hosted, Windows, X64, 100` | CPU builds and heavy contracts; either compatible heavy runner may execute them. |
| Light contracts | `self-hosted, Windows, X64, ci-light` | Source-only contracts with the audited lightweight dependency closure. |
| A100 compute | `self-hosted, Windows, X64, gpu-a100` | Nightly/manual GPU validation only. A100 40 GB; MATLAB R2026a; no compiler or builds. |

The build host compiles CUDA payloads and executes CPU checks. GPU numerical acceptance runs on A100.

The GPU job audit identified these actual consumers:

| Job | Before | After | Evidence in steps |
| --- | --- | --- | --- |
| radia-simulink.yml:cuda-mex | heavy + RTX | heavy | Compiles CUDA/MEX and seals the compute bundle. |
| radia-simulink-a100.yml:validate | absent | gpu-a100 | Validates the sealed bundle. |
| `gpu-validation.yml:validate` | absent | `gpu-a100` | installed GPU/Simulink acceptance, larger CSR solves and Optuna GPU oracle driver |
| `gpu-validation.yml:nightly-failure` | absent | `ubuntu-latest` | create/update failure issue |
| `hdiv-a100-validation.yml:validate` | absent | `gpu-a100` | HDiv device tests on the installed exact-SHA radia wheel |

`HDiv A100 validation` runs only through workflow dispatch, at a SHA whose `Radia Native Release` succeeded: that workflow seals its radia wheel as `gpu-radia-wheel`, and this one recovers it, installs it into a run-local environment and runs every test gated on `RADIA_HDIV_NATIVE_TEST_DIR`. A skipped, failed or missing case fails the run. It does not run nightly because Native Release does not build on main pushes, and it is not a required merge check.

`GPU compute validation` runs only nightly or through workflow dispatch. Dispatch at the release ref to validate that exact release SHA. It does not run on PRs and must not be a required merge check. Releases continue to use their existing gates; release management may explicitly request this dispatch and require its successful result for that SHA. No implicit A100 release dependency is introduced.

The build jobs seal wheel SHA-256 values and the source commit in `gpu-artifact.json`, then upload dedicated Actions artifacts for 30 days. The compute job searches successful push runs of the owning workflows at its exact SHA, downloads the two dedicated artifacts, verifies source identity and all wheel hashes before installation, and verifies the installed Simulink payload manifest. Missing, expired, stale or modified artifacts fail; receipt reuse from another SHA is forbidden. Both distribution build workflows must have run successfully for that SHA before dispatch. Scheduled runs select the newest commit on main with both unexpired sealed artifacts from successful push runs, check out validation scripts at that SHA and record it in the run summary and artifact receipts. Dispatch remains strict to its requested SHA. Schedules fail visibly when no complete pair exists; they never start a build.

Compute validation invokes the existing installed GPU/Simulink and Optuna validation drivers, with 131,072 density samples, plus bounded analytical SPD CSR cases of 32,768, 131,072 and 524,288 rows through the existing GPU solver. JSON records include driver/toolkit/MATLAB versions, device and memory, cold/warm transfer-inclusive timing and independently recomputed true residuals (limit 1e-6). These are validation timings, not speedup claims. The compute job does not certify the full AMS/FMM physics suite; that suite runs in the A100 Simulink acceptance workflow. It installs neither a compiler nor NGSolve and has no CPU fallback. Binary architecture/ABI failures remain failures.

Evidence uploads run even after failure. A hosted reporting job uses the same create-or-comment issue pattern as the nightly AMS audit for scheduled failures only, excluding cancellations. Dispatch and schedule use separate concurrency groups. An offline compute runner can remain queued; it does not block PRs. Humans own compute-host availability.
