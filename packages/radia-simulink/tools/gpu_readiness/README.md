# Manual runner readiness

The manager dispatches `gpu-runner-readiness.yml` once on an admitted idle GPU.
It has no push, PR or schedule trigger. `target` selects only existing `gpu-a100`
or `gpu-rtx`; this does not approve or add shared functional/availability labels.
Before dispatch the manager supplies secret `GPU_READINESS_A100` (or later
`GPU_READINESS_RTX`), a JSON object with string fields `runner_root`,
`expected_sid`, `expected_runner_name`, `matlab_root`, `python`, `expected_gpu`.
These identify the existing runner, authenticated account, R2026a installation,
already Engine-capable interpreter and GPU. Missing configuration fails. The
configuration is a single masked secret: private identities are not dispatch
inputs or ordinary job environment values printed by Actions. No secrets or
repository/runner settings are created or changed by this candidate.
No setup-python/pip, installation, license activation, compilation, benchmark,
service change, process attach or kill is performed.

The probe requires its actual Runner.Worker and Runner.Listener ancestors in the
approved interactive account/session (session zero is rejected), refuses active
GPU compute processes, selects VS2022 with vswhere [17,18), and checks installed
MSVC 14.44.35207. This is read-only discovery, not compiler qualification. CUDA
build early-return paths do not inherit Build.ps1 CPU checks; no build portability
or CUDA recipe parity is certified here.

A fresh MATLAB `-wait -batch` process and a fresh Python Engine process each
check installation/release, MATLAB/Simulink/PCT license test and checkout, and
`gpuDevice` without an index (no explicit device reset). The Engine quits only
its own session. No pool, numerical kernel or performance timing is requested.
License checkout and creating a CUDA context are transient effects of readiness.

Each probe step has a five-minute Actions timeout and the job ten minutes.
Missing items, nonzero exit, mismatched batch/Engine GPU/install identity or
shutdown error fail. Timeout/cancellation never produces passing evidence:
`dispatch.json` remains incomplete unless `summary.json` says finished/passed.
Actions enforces timeout; scripts have no process-kill operation. If MATLAB startup
or shutdown hangs, manager must inspect for owned orphan processes before another
dispatch; a hard cancellation cannot promise graceful MATLAB shutdown.

Artifacts contain only dispatch/summary JSON, including failed/incomplete attempts;
raw MATLAB logs and Engine JSON stay in the run-owned RUNNER_TEMP directory for
manager recovery/cleanup. Do not upload raw files: they contain account/path/license
information. Context paths/SID/runner name and GPU UUID are hashed in summary.
The artifact is readiness evidence, never exact-wheel numerical acceptance.

Concurrency serializes this diagnostic per target only. It does not reserve the
GPU against other workflows or externally admitted owners. The manager must admit
an idle window and drain other GPU work before dispatch. Shared labels and the
all-workflow owner-window drain proposal remain unapproved until this probe passes.

References: [GPUDevice properties](https://www.mathworks.com/help/parallel-computing/parallel.gpu.gpudevice.html),
[MATLAB Engine startup](https://www.mathworks.com/help/matlab/apiref/matlab.engine.start_matlab.html).
