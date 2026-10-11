---
name: release-quad
description: Two-machine CPU numerical Radia solver release gate. Use for release-quad, solver post-release deployment, Radia GitHub Release publication, or the numerical solver Definition of Done across LAB/102 and 100号機. radia-mcp and cubit-mesh-export use independent release-dual skills.
---

# release-quad

## Scope

`tools/release_quad.py` owns CPU Python/pybind11 `radia` solver QUAD acceptance.
Normal `radia` excludes mexcuda, MATLAB, LTspice and Simulink artifacts and
their CI gates; those belong to independent distribution releases. It does not publish, install, uninstall, repoint, reconnect,
or version-gate `radia-mcp` or `cubit-mesh-export`.

- `radia-mcp`: follow [package release checks](../../../packages/radia-mcp/CONTRIBUTING.md)
  and the [LAB/100 release-dual completion contract](../../../packages/radia-mcp/docs/operations/mcp-runtime-policy.md#release-completion).
- `cubit-mesh-export`: use `release-cubit-mesh-export` and LAB/100 release-dual.
- `radia-optuna`: use its exact-wheel candidate/done lane; do not install
  Radia, Cubit, or radia-mcp as a side effect.

Never restore an older checkout because its path was once canonical. Releases
and editable development move forward to an explicitly selected current source.

## Canonical solver commands

```powershell
python tools/release_quad.py preflight
python tools/release_quad.py phase8 --target lab,100
python tools/release_quad.py phase9
python tools/release_quad.py all
python tools/release_quad.py done
```

The retired `phase0` Cubit build and `restore-editable` routes are not part of
the solver workflow. Cubit native preparation belongs to the independent Cubit
release skill. Editable source advancement is explicit and forward-only.

### Retired Omega Override Gate

After the Kelvin source commits through `7ffda79d3` are integrated into main,
run `python tools/release_quad.py temp-shadows --apply`. This checks the exact
`C:\temp\radia-omega-test` tree on mdx1, mdx2, and hibino. It refuses removal
while Python/MATLAB processes are active, PYTHONPATH still names the tree, or
reparse points are present. It never terminates research processes or removes
other scratch directories. Unreachable hosts are unresolved, not clean. Unconfigured hibino and retired
mdx2 have not been checked; this standalone cleanup may fail on their aliases.
This standalone cleanup is not called by `all` or `done`.
Keep the JSON console report with the release evidence. hibino is checked for
this retired override only; it is not added to the four deployment targets.

## Solver machine policy

Decision 2026-10-01 supersedes the older two-editable-host recipe: LAB consumes
verified wheels; only 100 has a dedicated local editable development venv.
Keep its release runtime separate from in-progress development. `verify-editable`
checks LAB's wheel and 100's editable; `repoint` permits only 100.
100号機 owns CI/builds and isolated heavy numerical validation. 102号機 owns wheel acceptance.
The consolidated mdx host is reserved for idle timing evaluations; mdx1/mdx2/hibino are retired.


| Machine | Solver install tier | Solver release route |
|---|---|---|
| 102号機 | exact accepted Radia wheel | `phase8 --target lab` (wheel over SSH) |
| 100号機 | release runtime (machine-default Python): exact accepted Radia wheel; dedicated editable development venv `W:\00_CAE\Radia\environments\development` | `phase8 --target 100`: wheel into the release runtime, then editable into that venv only; `done`/`phase9` verify both separately |

100号機 owns isolated CI/builds and a separate editable development runtime.
102号機 tests fixed wheels with job-owned scratch. Formal timing belongs only
to otherwise idle mdx; compute timing is not a release deployment gate.


`phase9` compares only the solver version and the declared tracked solver-file
hashes. MCP and Cubit installations are independent observations, never solver
drift and never a QUAD blocker.

Deployments and final dependency checks address each target over SSH,
including LAB. Run the controller on a host with the configured `102` and `100` aliases; the controller's own Python is not a substitute for
a target runtime. A LAB-only controller restriction is no longer necessary.

## WIP-safe release sources

Never stash, reset, clean, or rebase a shared worktree for a release. Run the
controller from an independent clone (not a worktree of a shared checkout, whose
diverged local `main` would fail the main-sync gate). By default `done` verifies
that controller as the release source, so it must sit at the release tag. When
release tooling was repaired after the tag, keep the controller on current
`main` and pass a separate exact-tag checkout with `--release-source`: the
source must be tracked-clean at the tag commit, and the controller must be
tracked-clean, descend from that commit, and declare the same Radia version.
Neither side is retagged, and untracked tool copies are never substituted. The
LAB wheel needs no source tree. For 100号機's development venv, create one
tracked-clean release worktree containing the exact wheel's native payloads and
`src/radia/release_native_payloads.json`, and name its 100-local view:

```powershell
$env:RADIA_RELEASE_EDITABLE_REPO_100 = "W:\00_CAE\Radia\release-quad\<release>"
python tools/release_quad.py all
python tools/release_quad.py done --release-source <exact-tag checkout>
```

QUAD verifies exact SHA and tracked cleanliness before mutation. `done` is
non-mutating and leaves both runtimes unchanged. Later advance the development
venv's explicitly intended source to current `main` and verify it; do not use an
old-path restore operation.

## Optional mdx refresh

`phase8e` refreshes the installed Radia on SSH alias `mdx` only when explicitly
requested, using the controller's chosen released version. It is not called
by `all`, `phase9` or `done`, and is not a release acceptance target.

## Independent radia-optuna lane

```powershell
python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all
python tools/release_quad.py optuna-done --wheel <path>
```

The candidate is the exact successful-main wheel. Bind its SHA256, source SHA,
version, CI run, and all declared target results before tagging and publishing.
Do not run solver Phase 8 or install other distributions.

Each independent Optuna target runs MATLAB through its Engine worker. The runner's isolated venv installs `matlabengine` matching that host's
MATLAB version from PyPI and refuses an Engine bound to another MATLAB root.
Without a session it starts MATLAB only when the host has no MATLAB process
or shared Engine; otherwise it fails before launch. Name an existing shared
session per target with `--engine-session <target>=<name>` (repeatable). That
session is PID-checked, refused with loaded diagrams or Radia/Optuna MEX, keeps
its path and base workspace, and is not quit. The state records each target's
`engine_session`.

## Radia exact-artifact publication hold

Radia tag builds do not automatically publish to PyPI. Accept the exact tag-CI
wheel on the required solver validation hosts, commit the required machine-
readable evidence, and dispatch the release workflow with the exact CI run,
wheel SHA256, and acceptance commit. Promote the verified artifact without
rebuilding. Changed bytes require new acceptance.

CPU `ci-verify` requires the exact-SHA native release check and successful
shared/CPU checks. Independently owned GPU, MATLAB/LTspice, Optuna, MCP,
Cubit and Eqnedit64 workflows do not block CPU publication. Workflow ownership
is resolved from GitHub Actions run metadata; lookup failures fail closed.
GPU artifacts use `radia-gpu-release.yml` and their own acceptance evidence.

## Completion rules

- Do not call the numerical solver release complete until `done` exits 0.
- Do not call radia-optuna complete until `optuna-done` exits 0.
- Do not treat radia-mcp or Cubit client restart as a solver release gate.
- Do not infer live acceptance from editable metadata alone.
- Preserve other users' work and active jobs; never mass-kill MCP, Python,
  Cubit, MATLAB, Codex, or Claude processes.
- Keep this skill, `AGENTS.md`, `CLAUDE.md`, `tools/release_quad.py`, and
  `radia_mcp.radia_ngsolve.release_workflow` synchronized.
- Use live MCP catalog discovery; do not restore a tracked generated `docs/TOOLS.md` gate.
