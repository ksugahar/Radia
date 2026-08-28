---
name: release-eqnedit64
description: "Release and deploy Eqnedit64 from one exact tag, including signed standalone EXE verification, GitHub/PyPI publication checks, and byte-identical synchronization to O:\\Eqnedit64.exe. Use for Eqnedit64 releases, publication, or updating the shared O: copy."
---

# Release Eqnedit64

Finish with three matching outcomes: the tagged GitHub Release, the PyPI wheels,
and `O:\Eqnedit64.exe`. Do not call the release complete while any one is stale.

## Release invariants

- Build the standalone executable on LAB from the exact release tag with
  `tools/eqnedit64/build/build_eqnedt64.bat`. The signing key is deliberately
  non-exportable; the CMake executable produced by GitHub-hosted CI is an
  unsigned test artifact and must not be published or copied to `O:`.
- Require product version equal to the tag version, Authenticode status
  `Valid`, signer subject exactly `CN=ksugahar`, and a recorded SHA-256.
- Attach only that signed `Eqnedit64.exe` and its `SHA256SUMS.txt` to the
  `eqnedit64-v<version>` GitHub Release.
- Require the exact tag's Eqnedit64 CI to pass before publication. The PyPI
  workflow downloads the signed GitHub Release EXE and builds Python-specific
  wheels; wait for Python 3.10--3.13 verification and trusted publication.
- Run full private-font process suites only in the isolated Eqnedit64 CI/VM
  session. Do not bypass `EQNEDIT64_ISOLATED_TEST_SESSION` on the interactive
  LAB desktop.

## Synchronize O:

After the GitHub Release is public, run:

```powershell
pwsh -File .agents/skills/release-eqnedit64/scripts/sync_to_o.ps1 `
  -Tag eqnedit64-v<version>
```

The helper downloads the public release rather than trusting a local build,
checks the published checksum, version, and Authenticode signature, replaces
only `O:\Eqnedit64.exe`, then verifies byte identity at the destination. Use
`-WhatIf` for a read-only preflight. If replacement fails because the shared
EXE is locked, leave the existing file intact and report the lock; do not kill
unrelated or remote user processes.

An explicit release/deploy/update request authorizes the matching `O:` update.
For a request that only asks to inspect or test Eqnedit64, obtain authorization
before changing `O:`.

Report the tag, GitHub Release URL, PyPI version, destination path, SHA-256,
product version, and signer status.
