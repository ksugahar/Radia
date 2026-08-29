---
name: release-eqnedit64
description: "Release and deploy Eqnedit64 from one exact tag, including signed standalone EXE verification, GitHub/PyPI publication checks, and byte-identical synchronization to O:\\Eqnedit64.exe. Use for Eqnedit64 releases, publication, or updating the shared O: copy."
---

# Hand-test and release Eqnedit64

Finish with three matching outcomes: `O:\Eqnedit64.exe`, the tagged GitHub
Release, and the PyPI wheels. Do not call the release complete while any one is
stale.

## O: hand-test invariant

- `O:\Eqnedit64.exe` is Sugahara's canonical hand-test entry point. After a
  candidate is committed, developer-signed, and passes the local background
  tests, update O: **before asking Sugahara to test it**. Do not leave the only
  testable candidate under `C:\temp` or a worktree path.
- Hand-test staging does not require a version bump, a push to `main`, or a
  release tag. It does require a clean committed source, a build stamp matching
  that exact commit, a valid `CN=ksugahar` signature, and byte-identical copy.
- Use `scripts/sync_handtest_to_o.ps1`. It writes
  `O:\Eqnedit64.handtest.json` and preserves the former formal release manifest
  as `O:\Eqnedit64.last-release.json`. While a hand-test candidate is current,
  `O:\Eqnedit64.release.json` must be absent so tag publication cannot mistake
  the candidate for a release artifact.
- A later formal `sync_to_o.ps1` recreates `Eqnedit64.release.json` and removes
  the hand-test marker after all release-only checks pass.

## Release invariants

- Push the release commit to `main` first and require the main Eqnedit64 CI to
  pass. Do not create or push the release tag yet.
- Build the standalone executable on LAB from that exact `origin/main` commit
  with `tools/eqnedit64/build/build_eqnedt64.bat`. The signing key is
  deliberately non-exportable; the CMake executable produced by GitHub-hosted
  CI is an unsigned test artifact and must not be published or copied to `O:`.
- Require product version equal to the tag version, Authenticode status
  `Valid`, signer subject exactly `CN=ksugahar`, and a recorded SHA-256.
- Copy the verified executable to `O:\Eqnedit64.exe` and write the adjacent
  `O:\Eqnedit64.release.json` gate manifest. Only after both pass byte-for-byte
  verification may the exact recorded commit receive the release tag.
- Push `eqnedit64-v<version>` last. The tag CI verifies the application; the
  release workflow reads the already-staged O: executable through its short
  self-hosted job, requires the manifest source SHA to equal the tag SHA, then
  builds Python 3.10--3.13 wheels and publishes GitHub Release and PyPI from
  that one signed binary.
- Run full private-font process suites only in the isolated Eqnedit64 CI/VM
  session. Do not bypass `EQNEDIT64_ISOLATED_TEST_SESSION` on the interactive
  LAB desktop.

## Required order

1. Commit and push the versioned release source to `main`.
2. Wait for the main Eqnedit64 CI and Policy Lint to pass.
3. Confirm tracked files are clean and `HEAD == origin/main`.
4. Run `tools/eqnedit64/build/build_eqnedt64.bat` on LAB.
5. Before creating any release tag, stage the exact signed build to O::

```powershell
pwsh -File .agents/skills/release-eqnedit64/scripts/sync_to_o.ps1 `
  -Tag eqnedit64-v<version> `
  -SourceExe tools/eqnedit64/dist/Eqnedit64.exe `
  -SourceSha (git rev-parse HEAD)
```

6. Verify the helper reports `Updated=True`, the intended version, `Valid`,
   `CN=ksugahar`, and the recorded SHA-256.
7. Create the annotated tag at the manifest's exact source SHA and push it.
8. Wait for tag CI, GitHub Release, and PyPI publication to pass.

The helper refuses an existing remote release tag, a dirty or unpushed source,
a mismatched build stamp, version, signature, or hash. Use `-WhatIf` for a
preflight that does not change O:. If replacement fails because the shared EXE
is locked, leave the existing file intact and report the lock; do not kill
unrelated or remote user processes.

An explicit release/deploy/update request authorizes the matching `O:` update.
For a request that only asks to inspect or test Eqnedit64, obtain authorization
before changing `O:`.

For an authorized hand-test update, run:

```powershell
pwsh -File .agents/skills/release-eqnedit64/scripts/sync_handtest_to_o.ps1 `
  -SourceExe tools/eqnedit64/dist/Eqnedit64.exe `
  -SourceSha (git rev-parse HEAD)
```

Report the O: destination, hand-test manifest, source SHA, SHA-256, product
version, and signer status. A hand-test update is not a public release.

Report the pushed main SHA, O: manifest, tag, GitHub Release URL, PyPI version,
destination path, SHA-256, product version, and signer status.
