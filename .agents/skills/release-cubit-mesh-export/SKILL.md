---
name: release-cubit-mesh-export
description: Release cubit-mesh-export and deploy the same published artifacts to the editable development and wheel acceptance hosts, with native provenance, dual-host verification, handler acceptance and cleanup.
---

# release-cubit-mesh-export

cubit-mesh-export releases independently. Its deployment uses
`tools/release_cubit_dual.py`; solver releases must not install, repoint or
version-gate it. Cubit runs only on licensed Cubit hosts, not compute workers
that consume checked mesh artifacts.

## Host roles and private configuration

| Role | CLI target key | Deployment tier |
|---|---|---|
| Development/build/controller host | `100` | editable |
| Acceptance host | `lab` | wheel |

These are tool identifiers, not public hostnames or a physical host map.
Resolve machines, SSH aliases, users, local paths and durable evidence roots
from private operator configuration through management. Keep that mapping
outside tracked source, tests and release artifacts. Verify the current
hostname privately against the selected role before acting.

Run dual-host orchestration from the development host (`--run-from 100`).
Remote SSH tokens may not have interactive mapped drives or share credentials;
use target-local paths for remote work. Nested SSH with piped stdio from the
wheel host is not the supported controller route. Verify the configured alias
with `--ssh-alias-lab <wheel-alias>` / `--ssh-alias-100 <editable-alias>`;
resolution fails explicitly before remote action. Never silently pick another
machine. Host routing and privacy do not waive release gates.

Codex implements/validates and returns exact commits and evidence to Codex
management. Claude management owns review, push/PR/CI/merge/tag/publication and
deployment through the counterpart hub. Agree one executor per outward action.

## 1. Release candidate

1. Use a clean checkout from reviewed main on the build host. Update the
   package version in `packages/cubit-mesh-export/pyproject.toml`, the package
   `src/cubit_mesh_export/__init__.py`, and
   `src/cubit_plugin/cubit_mesh_export_pybind.cpp`, plus CHANGELOG. Commit source.
2. Use current `src/cubit_plugin/cubit_build.ps1 -Rebuild`. Commit the `.ccm`
   and `native_payloads.json`; its source commit identifies the native inputs.
   After rewritten commits, reconcile provenance and verify it again.
   After each native rebuild, inspect both the `.ccm` and curver `.pyd` with
   `dumpbin /dependents` and `dumpbin /imports`. Record the imported VC runtime
   DLLs/symbols and prove they resolve against the runtime shipped with Cubit
   2026.8 and Python 3.12 (the reviewed baseline is VC runtime 14.42 in both).
   Keep actual runtime versions, loaded paths and hashes in the evidence; the
   selected compiler/toolset alone does not prove runtime ABI compatibility.
   On both deployment host roles, require normal plugin loading, all six export
   formats, check-vol, curver import and actual curving, and `.vol` byte parity
   against the previous published release for the declared unchanged cases.
   Use isolated candidate runtimes and the shared handler-equivalence contract;
   do not replace the deployed release to probe compatibility. Any failure is a
   blocker: report it to management and never silently change the pinned
   toolset, substitute a backend or replace Cubit's bundled runtime DLLs.
3. Claude uploads the content-addressed pyd asset named by the manifest to
   the `binaries` release, downloads it back and verifies size/SHA-256. CI needs
   that exact asset. Use an operator-selected scratch path, not a fixed host path.
4. Validate `tests/test_cubit_mesh_quality.py`, `tests/test_cubit_installers.py`,
   `tests/test_release_cubit_dual.py`, `tests/test_native_build_provenance.py` and
   `packages/cubit-mesh-export/tests`. Missing required candidate binaries fail.
   For native changes run the owning geometric refit benchmark with an absolute
   `--plugin-dir` and compare outputs against the previous published release.
5. Hand the candidate to management for reciprocal review and CI/main merge.
   After exact-release main CI succeeds, Claude pushes the annotated tag
   `cubit-mesh-export-vN`; its workflow publishes to PyPI. Download the wheel
   from PyPI and match its SHA-256 with published metadata before deployment.

## 2. Release checkout and wheel

Use private configuration to set `$Repository`, `$ReleaseRoot`, `$Wheel`,
`$EvidenceWheel`, `$EvidenceEditable` and `$WheelAlias` in the controller's
PowerShell. All paths must be absolute and verified for the selected role.
Create the development release worktree at the exact release tag, copy the
hash-verified pyd into the owning package, and verify `_native_provenance.py`.
The wheel host needs no source checkout. Transported wheel bytes are rehashed,
installed without dependency changes and verified in a fresh interpreter.
The editable checkout's packaged files must equal those same published bytes.

## 3. Dual-host deploy

```powershell
python tools/release_cubit_dual.py --action deploy --run-from 100 `
  --wheel $Wheel --source-sha <tag-commit> --source-root-100 $ReleaseRoot `
  --ssh-alias-lab $WheelAlias `
  --evidence-lab $EvidenceWheel --evidence-100 $EvidenceEditable
```

Then run the same arguments with `--action done`. Both preflights finish
before either install. Target keys bind receipts to tiers: `lab` is wheel,
`100` is editable with `--all-users`. The controller imports installer helpers
from the selected verified payload; a prior install is unnecessary and an old
editable package must not shadow it. Receipts record exact identity, package
preservation, exports, MCP selftests and registration verification. `done`
rejects missing/stale/wrong-tier receipts. Archive operational evidence privately.

Normally running Cubit stops deployment. Only an explicit instruction for that
deployment permits `--force-close-cubit`; it terminates all detected
`coreform_cubit.exe`, `clarox.exe` and `cubitx.exe` processes on both targets and
records image/PID/user, with unsaved work loss. Never infer permission from a
release request or acceptance-policy change. Preserve other users' jobs.

## 4. Acceptance and post-release cleanup

For every release, Codex and Claude use the shared
[cubit-toolbar-acceptance](../cubit-toolbar-acceptance/SKILL.md) contract.
Local skills must reference this repository contract rather than carry
independent acceptance requirements.

- On the editable and wheel hosts, verify the release/import/native identity and
  installed startup/menu/toolbar registration, then all six production GUI
  handler/command-builder paths, licensed APREPRO/Python backend execution
  with the same options, and validated fresh output artifacts. Reuse exact
  existing release evidence when it satisfies the contract.
- Record `cubit-mesh-export.toolbar-acceptance.v2` and policy
  `aprepro-python-gui-handler-equivalence-2026-10-16`. Both hosts must pass.
  Unmeasured mouse dispatch, visible dialogs, Cancel, restart dispatch,
  normal-GUI automatic discovery and BBND GUI behavior stay `null`, with
  scope in separate fields. These are optional diagnostics, not additional
  release gates. Do not claim them from handler-equivalence evidence.
- Never manually import a missing toolbar to call normal registration PASS.
  Keep MCP Cubit execution batch/nographics. Do not create accounts, activate
  licenses or force-close user sessions merely to perform acceptance tests.
- Management assigns one cleanup executor within the human-authorized scope.
  Recover inputs, scripts, commands, outputs, successful/failed logs and
  runtime identity to durable private evidence; verify source/destination
  SHA-256 before deleting anything. Preserve the selected deployed checkout,
  wheel and evidence.
- Check exact absolute targets, active consumers, unique untracked/ignored
  work and escaping links/junctions. Close only owned sessions normally;
  delete approved scratch with literal paths and verify absence. Obsolete
  branch/worktree retirement additionally checks main inclusion and ownership.
  No blanket temp cleanup, unrelated deletion or implicit force-kill.
- Save an executor/timestamp/target/hash/precheck/outcome cleanup receipt.
  A rejected deletion is not retried through another tool or agent to bypass
  policy. Retained work needs an owner, reason and cleanup trigger.
- Codex hands verified results to its management hub. Claude management
  confirms acceptance and cleanup through the counterpart hub. Close the
  release only after actual confirmation, a final hash inventory check and
  cleanup completion (or explicitly report the separate retained cleanup item).
  A sent review packet or passing CI alone is not a closure receipt.
