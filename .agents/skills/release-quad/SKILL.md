---
name: release-quad
description: Four-machine Radia release gate. Use when the user asks for release-quad, release_quad, post-release deploy, GitHub Release publication, or the Definition Of Done for Radia releases. Coordinates PyPI and Simulink release candidates across LAB, 100号機, mdx1, and mdx2 via tools/release_quad.py.
---

# release-quad

## Canonical Entry Point

Use only `tools/release_quad.py`; release work must go through QUAD.

```powershell
python tools/release_quad.py preflight
python tools/release_quad.py phase0
python tools/release_quad.py phase8 --target lab,100
python tools/release_quad.py phase8e
python tools/release_quad.py phase9
python tools/release_quad.py simulink-candidate --package <zip> --target all
python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all
python tools/release_quad.py optuna-done --wheel <path>
python tools/release_quad.py all
python tools/release_quad.py done --simulink-package <zip>
python tools/release_quad.py restore-editable
```

## Machine Policy

| Machine | Install tier | Release command path |
|---|---|---|
| LAB | NAS editable | `phase8 --target lab` |
| 100号機 | NAS editable over SSH | `phase8 --target 100` |
| mdx1 | Radia PyPI wheel consumer over `ssh mdx1`, no MCP/Cubit deployment | `phase8e` |
| mdx2 | Radia PyPI wheel consumer over `ssh mdx2`, no MCP/Cubit deployment | `phase8e` |

Deploy `cubit-mesh-export` and the Cubit plugin/toolbar only to LAB and
100号機. The mdx lane neither installs the package nor runs Cubit commands,
even if Cubit is detected there. Do not uninstall existing Cubit packages as
part of this routing change. Isolated Cubit-independent checker CI is separate
from deployment and remains allowed.

`radia-optuna` is an independent release lane. It does not run `phase8` or
install Radia/Cubit. `optuna-candidate` downloads the exact wheel artifact from
one successful `main` push CI run and runs that wheel's installed-wheel
MATLAB/Simulink test on LAB, 100号機, mdx1, and mdx2. `optuna-done` requires
the wheel SHA256, source commit, package version, CI run, and all four target
results to agree.

`phase9` is the hard gate: LAB / 100号機 / mdx1 / mdx2 must agree on
versions, compatibility constants, and tracked file hashes. Both mdx hosts report
`radia-mcp`, `cubit-mesh-export`, and Cubit compatibility fields as `N/A`;
that is intentional and is excluded from drift
comparison.

Both mdx runners use the shared `mdx` label and a unique `mdx1` or `mdx2`
label. Before accepting a replacement runner, verify Python 3.12, Git,
PowerShell 7, MSVC x64 via `vswhere`, Windows SDK, CMake, MATLAB, and the
MATLAB Engine import. A successful fast-contract job does not certify MEX
build readiness: also run the Optuna native build lane. Keep pip build
dependencies in the workflow's run-local venv. Remove retired runner
registrations after their replacements are verified online.

MATLAB-capable mdx runner services use `.\Administrator`, the current
MATLAB-authenticated account, not SYSTEM. Configure credentials through the
Windows Services Log On UI, never chat, command arguments, or committed files.
The visible service is `GitHub Actions Runner (ksugahar-Radia.mdx1-ci)`
(or mdx2-ci). An alternate least-privilege account requires its own MATLAB
authentication and service-context acceptance; administrator privilege itself
is not the MATLAB requirement.

After runner reinstall/account changes, manually dispatch
`matlab-engine-diagnostic.yml`. It tests both actual runner accounts with a
dedicated Python Engine session, calculation and shutdown, with a 60-second
owned-process-tree timeout. Do not add this license-dependent diagnostic to
normal push/tag CI. SSH success or an open RDP MATLAB is not runner acceptance.
The initial SYSTEM failure and subsequent account-change success are recorded
in `validation_test/ci_runtime/matlab_engine_runner_20260908.json`.

mdx1 is scheduled for retirement in March 2027. Before retiring it, obtain the
replacement-target decision and update runner routing and release-quad target
contracts together; verify MATLAB authentication and the native build on the
replacement. Do not silently reduce the four-host gate or invent a replacement.

## Parallel WIP-Safe Editable Source

Do not stash, clean, reset, or rebase a shared LAB worktree just to release.
When `S:\Radia\01_GitHub` contains parallel work, create one clean release
worktree on the NAS and expose its two machine-local views before `all`,
`verify-editable`, and `done`:

```powershell
$env:RADIA_RELEASE_EDITABLE_REPO_LAB = "S:/Radia/release-quad/<release>"
$env:RADIA_RELEASE_EDITABLE_REPO_100 = "W:\00_CAE\Radia\release-quad\<release>"
python tools/release_quad.py all
```

The release worktree must contain the native build outputs needed by editable
installs. QUAD verifies its exact Git SHA and tracked-clean state on LAB and
100号機 before killing processes or installing anything. Repeat the same two
environment variables for `done`. The final `done` gate is non-mutating: it
requires the active LAB source to remain at the exact tracked-clean release SHA,
requires that SHA to equal the peeled `v<radia-version>` tag, verifies editable
metadata and actual import origins on LAB/100号機, and leaves those verified
pointers unchanged. After the canonical `01_GitHub` tree catches up with
published `main`, run `restore-editable` explicitly to return the development
tier to that tree.

## Rules

- Do not call a release done until `python tools/release_quad.py done`
  exits 0.
- Do not call a standalone radia-optuna release done until
  `python tools/release_quad.py optuna-done --wheel <path>` exits 0. Tag only
  afterward. Dispatch `release-radia-optuna.yml` with the exact candidate CI
  run ID and SHA256 emitted by QUAD; the workflow must publish that artifact to
  both PyPI and the GitHub Release without rebuilding it.
- Do not publish any GitHub Release containing the Radia Simulink library,
  MATLAB support files, or MEX assets until the complete four-machine gate
  passes for LAB, 100号機, mdx1, and mdx2.
- The Simulink gate hashes the exact ZIP and extracts it independently on each
  machine. A full library runs `verify_radia_simulink_release`; an IH preview
  runs `verify_radia_ih_release`. Rebuilding the ZIP invalidates the recorded
  gate state and requires all four checks again.
- Execute these checks through Python's MATLAB Engine, not a MATLAB `-batch`
  launcher. Each verifier owns a new Engine session and closes it on success
  or failure; do not attach to or close a user's shared session. The selected
  Python must have an Engine installation matching the requested MATLAB root.
  Keep the isolated worker timeout so a stalled startup is bounded without
  stopping unrelated MATLAB or solver processes.
- The `done` result is the authoritative publication gate; partial, failed,
  or manually waived machine checks do not authorize publication.
- Assemble and test the versioned Simulink package before publication. A full
  library package includes `radia_simulink_library.slx`; the standalone IH
  preview instead includes `radia_ih.slx` and only its native support files.
  Both forms include their applicable MEX assets, `manifest.json`, and
  `SHA256SUMS.txt`.
- This gate applies to every subsequent Simulink library revision as well as
  the initial release.
- Use `ssh mdx1` and `ssh mdx2` for the two compute targets. Follow
  the repository SSH policy: pipe a script into
  `ssh mdx1 'pwsh -ExecutionPolicy Bypass -Command -'`.
- Keep `packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/` and the
  checked `radia_mcp.meta` discovery catalog in sync when changing release or
  deploy knowledge. The live MCP registry is the authoritative tool inventory;
  do not restore a tracked generated `docs/TOOLS.md` gate.
