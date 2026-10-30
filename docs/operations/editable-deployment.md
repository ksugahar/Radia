# Editable operation and isolated wheel qualification

Radia releases have four cases: build-role editable operation, validation-role
editable operation, CPU compute build and GB200 GPU compute build. The
`release-quad` skill coordinates all four; build-only evidence does not certify
numerical execution. Host-specific skills and host maps remain private.

The two normal Windows runtimes use clean, explicitly selected editable sources
with same-source native artifacts. Each role also has a separate isolated venv
for exact-wheel qualification. A normal editable import is never evidence that
the released wheel works. MCP, Cubit and Simulink keep independent release and
runtime ownership.

The reviewed controller resolves private configuration for each suffix `100`
and `VALIDATION`:

| Setting prefix | Meaning |
| --- | --- |
| `RADIA_RELEASE_SSH_` | One configured SSH alias |
| `RADIA_RELEASE_HOSTNAME_` | Actual target hostname |
| `RADIA_RELEASE_PYTHON_` | Absolute normal-runtime Python executable |
| `RADIA_RELEASE_EDITABLE_REPO_` | Absolute persistent clean source checkout |
| `RADIA_RELEASE_EDITABLE_COMMIT_` | Exact release source commit |
| `RADIA_RELEASE_EDITABLE_PATH_ALIASES_` | Optional JSON pairs for equivalent private path prefixes |
| `RADIA_RELEASE_QA_PYTHON_` | Absolute isolated wheel-QA Python, separate from normal operation |

No PATH, WindowsApps, hostname substitution or implicit source selection is
accepted. Verify the actual remote host, interpreter and clean source before
installation or intent writes. Resolve all selected roles before mutation.

The assigned executor stages verified native artifacts from the exact release
source into each persistent checkout. Native binaries are not rebuilt on the
consumer. The normal install uses the selected interpreter and the source's
`[gpu]` extra. Verify metadata, dependency pins, native hashes and editable
`direct_url.json`, fresh import origin and `pip check` before recording the
Radia source intent. Preserve other packages' intent entries and environments.

Run exact-wheel qualification in each distinct QA venv. Retain wheel and RECORD
hashes, native source/ABI provenance, required adopted numerical cases and
explicit CPU/GPU acceptance. An import or successful installation alone is not
qualification. Never relabel branch-wheel preparation as final tag-wheel
acceptance. Publication also requires exact-release native CI and the original
scientific gates; preserve experimental labels and user-owned criteria.

The controller candidate must be independently reviewed and pass affected tests
and the shared precheck before use. A missing compute-build receipt keeps the
four-case wrapper incomplete even when both normal installs are healthy.
