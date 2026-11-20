# Validation-role verified-wheel deployment

A Radia release waits for four lane receipts of one source SHA: `cpu` (built
on the CPU-only build host), `linux-python-gpu` (built and validated on the
supercomputer), `windows-python-gpu` and `matlab` (both built and validated on
the Windows GPU validation host). The `release-quad` skill coordinates them;
build-only evidence does not certify numerical execution. Host-specific skills
and host maps remain private.

Editable development installs exist only on the build host. The validation
role runs verified wheels only: its normal Windows solver runtime holds the
exact `windows-python-gpu` lane wheel with the `[gpu]` extra, which the
CPU-only build host never installs. There is no editable route, source
binding or separate QA interpreter on the validation role. MCP, Cubit and
Simulink keep independent release and runtime ownership.

The reviewed controller resolves private configuration for the suffix
`VALIDATION`:

| Setting prefix | Meaning |
| --- | --- |
| `RADIA_RELEASE_SSH_` | One configured SSH alias |
| `RADIA_RELEASE_HOSTNAME_` | Actual target hostname |
| `RADIA_RELEASE_PYTHON_` | Absolute normal-runtime Python executable |

No PATH, WindowsApps or hostname substitution is accepted. The actual remote
host and interpreter are verified before any staging or installation.

`python tools/release_quad.py phase8 --source-sha SHA --manifest FILE`
deploys one release:

1. The source SHA must be the exact release tag commit, and the four lane
   receipts must verify against it.
2. The published lane must hold exactly one `radia-<version>-*.whl`; its
   SHA-256 is re-checked locally and again after it is staged on the target.
   A mismatch fails before pip runs.
3. pip installs the wheel's `[gpu]` dependencies from binaries only, then
   force-reinstalls the exact wheel without dependencies.
4. A fresh target interpreter refuses an editable `direct_url.json`, checks
   version, import origin under site-packages, every installed RECORD hash,
   that the installed hashes equal the candidate wheel's own RECORD, the
   metadata-pinned NGSolve/Netgen and GPU runtime pins, and `pip check`.
5. The staged wheel copy is removed and its removal verified.

`phase9` re-probes the runtime read-only: non-editable registration, unshadowed
import, the release version, ABI pins and native payload digest.

An import or successful installation alone is not qualification. Retain wheel
and RECORD hashes, native source/ABI provenance, required adopted numerical
cases and explicit CPU/GPU acceptance. Never relabel branch-wheel preparation
as final tag-wheel acceptance. Publication also requires exact-release native
CI and the original scientific gates; preserve experimental labels and
user-owned criteria.

The controller candidate must be independently reviewed and pass affected tests
and the shared precheck before use. A missing lane receipt keeps the release
incomplete even when the installed runtime is healthy.
