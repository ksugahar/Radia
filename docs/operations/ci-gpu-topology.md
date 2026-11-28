# CI runner roles

The CI pool has three roles. Capability labels select hardware; they are not distribution names.

| Role | Labels | Work |
| --- | --- | --- |
| CPU heavy | `self-hosted, Windows, X64, 100` | CPU builds and heavy contracts; either compatible heavy runner may execute them. |
| Light contracts | `self-hosted, Windows, X64, ci-light` | Source-only contracts with the audited lightweight dependency closure. |
| Windows GPU validation host | `self-hosted, Windows, X64, gpu-a100` | All Windows GPU builds (Python GPU wheel, CUDA MEX/MATLAB), their CI and release qualification. A100 40 GB; MATLAB R2026a; one GPU job at a time. |

The build host is CPU-only: CPU builds, CPU CI (`100`, `ci-light`), early debugging and CPU release-candidate preparation. It compiles no CUDA, installs no NVIDIA/cuDSS packages and runs no GPU or MATLAB GPU work; the `gpu-rtx` label is retired. The Windows GPU validation host (`gpu-a100`) owns GPU builds, GPU CI and GPU qualification. Linux GPU build debugging runs on the supercomputer outside CI. A release waits for four lane receipts of one SHA: `cpu`, `linux-python-gpu`, `windows-python-gpu` and `matlab`.

The GPU job audit identified these actual consumers:

| Job | Before | After | Evidence in steps |
| --- | --- | --- | --- |
| radia-simulink.yml:cuda-mex | heavy + GPU | `gpu-a100` | Compiles CUDA/MEX and seals the compute bundle on the Windows GPU host. |
| radia-simulink-a100.yml:validate | absent | gpu-a100 | Validates the sealed bundle. |
| `gpu-validation.yml:validate` | absent | `gpu-a100` | installed GPU/Simulink acceptance and larger CSR solves |
| `gpu-validation.yml:nightly-failure` | absent | `ubuntu-latest` | create/update failure issue |
| `hdiv-a100-validation.yml:validate` | absent | `gpu-a100` | HDiv device tests on the installed exact-SHA radia wheel |

`HDiv A100 validation` runs only through workflow dispatch, at a SHA whose `Radia Native Release` succeeded: that workflow seals its radia wheel as `gpu-radia-wheel`, and this one recovers it, installs it into a run-local environment and runs every test gated on `RADIA_HDIV_NATIVE_TEST_DIR`. A skipped, failed or missing case fails the run. It does not run nightly because Native Release does not build on main pushes, and it is not a required merge check.

`GPU compute validation` runs only nightly or through workflow dispatch. Dispatch at the release ref to validate that exact release SHA. It does not run on PRs and must not be a required merge check. Releases continue to use their existing gates; release management may explicitly request this dispatch and require its successful result for that SHA. No implicit A100 release dependency is introduced.

The build jobs seal wheel SHA-256 values and the source commit in `gpu-artifact.json`, then upload dedicated Actions artifacts for 30 days. The compute job searches successful push runs of the owning workflows at its exact SHA, downloads the two dedicated artifacts, verifies source identity and all wheel hashes before installation, and verifies the installed Simulink payload manifest. Missing, expired, stale or modified artifacts fail; receipt reuse from another SHA is forbidden. Both distribution build workflows must have run successfully for that SHA before dispatch. Scheduled runs select the newest commit on main with both unexpired sealed artifacts from successful push runs, check out validation scripts at that SHA and record it in the run summary and artifact receipts. Dispatch remains strict to its requested SHA. Schedules fail visibly when no complete pair exists; they never start a build.

Compute validation invokes the existing installed GPU/Simulink validation driver, plus bounded analytical SPD CSR cases of 32,768, 131,072 and 524,288 rows through the existing GPU solver. JSON records include driver/toolkit/MATLAB versions, device and memory, cold/warm transfer-inclusive timing and independently recomputed true residuals (limit 1e-6). These are validation timings, not speedup claims. The compute job does not certify the full AMS/FMM physics suite; that suite runs in the A100 Simulink acceptance workflow. It installs neither a compiler nor NGSolve and has no CPU fallback. Binary architecture/ABI failures remain failures.

Evidence uploads run even after failure. A hosted reporting job uses the same create-or-comment issue pattern as the nightly AMS audit for scheduled failures only, excluding cancellations. Dispatch and schedule use separate concurrency groups. An offline compute runner can remain queued; it does not block PRs. Humans own compute-host availability.
