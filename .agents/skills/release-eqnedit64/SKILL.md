---
name: release-eqnedit64
description: "Release Eqnedit64 EXE and Web/JS from one exact tag: signed EXE, O: synchronization, GitHub/PyPI publication, and verified JS publication on the Sugahara Laboratory homepage. Also supports O: hand-test updates."
---

# Hand-test and release Eqnedit64

Finish with four matching outcomes: `O:\Eqnedit64.exe`, the tagged GitHub
Release, the PyPI wheels, and the Web/JS edition published on the Sugahara
Laboratory homepage. All four must correspond to the same release tag.
GitHub/PyPI success alone does not complete an Eqnedit64 release.

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
- The later formal release synchronization recreates
  `Eqnedit64.release.json` and removes the hand-test marker after all
  release-only checks pass. Its exact command appears only in the required
  order below, after the main-branch gate.

## Release invariants

- Before merging the release candidate to `main`, require one Claude Code
  Fable review of the candidate commit and its specification diff. Record the
  result in the PR or handover. Address findings and repeat the affected tests
  and O: hand test before formal publication; do not create a tag from an
  unreviewed candidate.
- Record both the reviewed SHA and the finding-resolution SHA. A follow-up
  commit limited to the recorded findings and their tests completes the one
  review gate after affected CI and hand testing pass. Any unrelated model or
  specification change creates a new candidate and requires another review.
- After the hand-test and Fable gates, push the approved release commit to
  `main` and require the main Eqnedit64 CI to pass. Do not create or push the
  release tag yet.
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
  release workflow then downloads the staged executable from the
  `eqnedit64-staging` release on a GitHub-hosted runner, requires the manifest
  source SHA to equal the tag SHA, and builds Python 3.10--3.13 wheels and
  publishes GitHub Release and PyPI from that one signed binary. No job reads
  O:, and none runs on a self-hosted runner: a runner service is NETWORK
  SERVICE, which can reach neither a per-user mapped drive nor the workgroup
  share. The signature is what carries the guarantee across that transport -
  the key is non-exportable, so a matching `CN=ksugahar` signature, SHA-256,
  product version, and source SHA together prove the binary was built on LAB
  from the tagged commit, and an unsigned CI build cannot impersonate it.
- On the interactive LAB desktop, the release does exactly two things:
  **compile and sign**. `build_eqnedt64.bat` invokes the compiler and
  `signtool`; it starts no window, renders no equation, and registers no font,
  so it stays clear of the `fontdrvhost.exe` failure that leaves the session's
  Office fonts inkless. Everything that *runs* the built binary - the GUI,
  rendering, `--self-test`, the font-session endurance suites, and
  `accept_release.ps1` - belongs to the isolated CI or VM session. Do not
  bypass `EQNEDIT64_ISOLATED_TEST_SESSION` on the interactive LAB desktop.
  Compiling on LAB is not optional: the signing key is non-exportable, so the
  shipped executable can only be produced there.

## Required order

1. Commit and push the versioned release candidate branch.
2. Wait for its Eqnedit64 CI and Policy Lint to pass, build the exact signed
   candidate, and stage it with `sync_handtest_to_o.ps1`.
3. Complete Sugahara's O: hand test.
4. Run one Claude Code Fable review on that candidate and its specification
   diff. Record and resolve the findings, then refresh affected tests and O:
   hand testing. Do not merge to `main` before this gate is recorded.
5. Merge the approved release commit to `main` and push it.
6. Wait for the main Eqnedit64 CI and Policy Lint to pass.
7. Run the build, staging, and tag-push phase as one command. It checks the source against
   `origin/main`, refuses a tag that already exists, refuses a tree whose
   declared version disagrees with the tag, compiles and signs, verifies the
   product version and `CN=ksugahar` signature, backs up whatever O: currently
   holds - a hand-test candidate staged there is someone's work in progress -
   stages to O: and to the `eqnedit64-staging` release, and only then creates
   and pushes the annotated tag:

```powershell
pwsh -File .agents/skills/release-eqnedit64/scripts/release_eqnedit64.ps1 `
  -Tag eqnedit64-v<version>
```

   `-WhatIf` runs every check and the build and changes nothing; `-SkipBuild`
   reuses an existing `dist/Eqnedit64.exe`, which the identity checks still
   have to accept. Run the individual steps below only when diagnosing a
   failure of this one.

8. Wait for tag CI, GitHub Release, and PyPI publication to pass.
9. Publish and verify the Web/JS edition from that exact tag using the homepage
   procedure below. Record all four destinations before reporting completion.

## Web/JS homepage publication

Radia owns both editions. The canonical Web sources are
`tools/eqnedit64/web/equation-editor.fragment.html` and `equation-editor.js`.
The publication destination is
`https://www.ele.kindai.ac.jp/laboratory/sugahara/elemag/equation-editor.php`.

1. Use a clean checkout of the released tag as `RADIA_REPOSITORY`. Do not use
   a moving main branch or the shared backup-branch checkout. The homepage
   workspace is `W:\20_各委員\01_ホームページ\00_sitedata`; read its local
   `AGENTS.md` and `equation_handover.md` before publishing.
2. Update the EXE link in `05_3Dで学ぶ電磁気学/pages/equation-editor.md`
   to the tag's verified GitHub Release asset. Keep the visible label versionless
   (`Windows版をダウンロード`), the GitHub source link, and the Windows
   screenshot. The learning portal links to this independent page.
3. In the homepage workspace, run
   `site_builder/tools/ci_equation_editor.py --radia-repository <tag-checkout> --profile browser --verify-download`.
   This builds only the editor page/assets and checks the browser contract.
   Run the existing Office integration gate in the approved isolated session
   when shared clipboard behavior changed. Do not run Wolfram, Mathematica,
   NGSolve, or the undergraduate/advanced teaching suites for this publication.
4. Use `site_builder/tools/upload_winscp.ps1` sequentially for additive upload
   of `elemag/equation-editor.php`, `elemag/equation-editor.js`, and changed
   editor-specific images only. Pass `-Files` as an array from the same pwsh
   process. Use the currently verified server certificate fingerprint; never
   disable certificate checking. Update portal/sitemap only when their content
   changed. Keep deployment credentials in the homepage deployment setup.
5. Fetch the public page and JS, check the direct EXE link and public/source
   SHA-256 equality, and run `run_equation_editor_qa.ps1` against the public
   base URL with `-SkipPowerPointIntegration`. Record the tag SHA, published
   URLs, JS hash, browser result, and uploaded files in the homepage handover.
   If the public JS already matches the tag, verify it and record that fact;
   an unchanged JS file need not be uploaded again.

If homepage publication or verification fails, report the release as partially
published with the failing destination. Preserve successful EXE publication and
resume the homepage phase from the same tag. Do not tag a new version merely
to retry an upload. `release_eqnedit64.ps1` dispatches the binary release;
its return value is not proof of end-to-end publication.

The individual steps, for diagnosis: `tools/eqnedit64/build/build_eqnedt64.bat`
compiles and signs; then, before any tag exists,

```powershell
pwsh -File .agents/skills/release-eqnedit64/scripts/sync_to_o.ps1 `
  -Tag eqnedit64-v<version> `
  -SourceExe tools/eqnedit64/dist/Eqnedit64.exe `
  -SourceSha (git rev-parse HEAD)
```

must report `Updated=True`, the intended version, `Valid`, `CN=ksugahar`, the
recorded SHA-256, and the two `StagedAssets`. Create the annotated tag last, at
the manifest's exact source SHA, and push it.

The helper refuses an existing remote release tag, a dirty or unpushed source,
a mismatched build stamp, version, signature, or hash. Use `-WhatIf` for a
preflight that does not change O:. If replacement fails because the shared EXE
is locked, leave the existing file intact and report the lock; do not kill
unrelated or remote user processes.

An explicit release/deploy/update request authorizes the matching `O:` update.
For a request that only asks to inspect or test Eqnedit64, obtain authorization
before changing `O:`.

## Who may run this

Eqnedit64 is the exception to the usual boundary where an agent stops at the
local commit and hands the push, the tag and the publication to another agent.
Here one authorized request carries the whole release: build and sign, stage to
O: and to `eqnedit64-staging`, push the tag, and see the GitHub Release and the
PyPI wheels appear, then publish and verify the matching Web edition on the
laboratory homepage. A formal Eqnedit64 release request includes this homepage
phase (user direction, 2026-09-08). An O: hand-test-only request does not.
Splitting it was worse than keeping it together - the
sequence is order-critical and every step verifies the one before it, so the
agent that produced the signed binary is the one that should carry it to the
end and report what actually published. Sugahara authorized this on
2026-09-05, after Claude Code took 3.0.15 from source to PyPI that way.

This does not widen anything else. The release still needs an explicit request,
the hand-test and review gates still come first, and the boundary for other
Radia work is unchanged.

For an authorized hand-test update, run:

```powershell
pwsh -File .agents/skills/release-eqnedit64/scripts/sync_handtest_to_o.ps1 `
  -SourceExe tools/eqnedit64/dist/Eqnedit64.exe `
  -SourceSha (git rev-parse HEAD)
```

Report the O: destination, hand-test manifest, source SHA, SHA-256, product
version, and signer status. A hand-test update is not a public release.

Report the pushed main SHA, O: manifest, tag, GitHub Release URL, PyPI version,
destination path, SHA-256, product version, signer status, homepage URL,
published JS SHA-256, and public browser verification result.
