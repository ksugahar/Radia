---
name: release-quad
description: Four-machine numerical Radia solver and Simulink release gate. Use for release-quad, solver post-release deployment, Radia GitHub Release publication, or the numerical solver Definition of Done across LAB, 100号機, mdx1, and mdx2. radia-mcp and cubit-mesh-export use independent release-dual skills.
---

# release-quad

## Scope

`tools/release_quad.py` owns the numerical `radia` solver and Radia Simulink
QUAD acceptance. It does not publish, install, uninstall, repoint, reconnect,
or version-gate `radia-mcp` or `cubit-mesh-export`.

- `radia-mcp`: use `release-radia-mcp` and LAB/100 release-dual.
- `cubit-mesh-export`: use `release-cubit-mesh-export` and LAB/100 release-dual.
- `radia-optuna`: use its exact-wheel candidate/done lane; do not install
  Radia, Cubit, or radia-mcp as a side effect.

Never restore an older checkout because its path was once canonical. Releases
and editable development move forward to an explicitly selected current source.

## Canonical solver commands

```powershell
python tools/release_quad.py preflight
python tools/release_quad.py phase8 --target lab,100
python tools/release_quad.py phase8e
python tools/release_quad.py phase9
python tools/release_quad.py simulink-candidate --package <zip> --target all
python tools/release_quad.py all
python tools/release_quad.py done --simulink-package <zip>
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
other scratch directories. Unreachable hosts are unresolved, not clean.
`all` runs this cleanup first; `done` repeats the read-only absence check.
Keep the JSON console report with the release evidence. hibino is checked for
this retired override only; it is not added to the four deployment targets.

## Solver machine policy

| Machine | Solver install tier | Solver release route |
|---|---|---|
| LAB | verified current editable | `phase8 --target lab` |
| 100号機 | verified current editable over SSH | `phase8 --target 100` |
| mdx1 | exact accepted Radia wheel | `phase8e` |
| mdx2 | exact accepted Radia wheel | `phase8e` |

LAB is the development host. 100号機 is the student-facing release/usage host;
keep acceptance there to installation, import, and necessary application smoke.
mdx1/mdx2 are the self-hosted CI and preflight pool. hibino is a computation
host, not a QUAD acceptance target.

`phase9` compares only the solver version and the declared tracked solver-file
hashes. MCP and Cubit installations are independent observations, never solver
drift and never a QUAD blocker.

## WIP-safe editable sources

Never stash, reset, clean, or rebase a shared worktree for a release. If the
ordinary LAB tree contains parallel work, create one tracked-clean current
release worktree on durable storage and expose its two host-local views:

```powershell
$env:RADIA_RELEASE_EDITABLE_REPO_LAB = "S:/Radia/release-quad/<release>"
$env:RADIA_RELEASE_EDITABLE_REPO_100 = "W:\00_CAE\Radia\release-quad\<release>"
python tools/release_quad.py all
python tools/release_quad.py done --simulink-package <zip>
```

The worktree contains the solver-native artifacts required by the editable
install. QUAD verifies exact SHA and tracked cleanliness before mutation. `done`
is non-mutating and leaves verified pointers unchanged. Later advance the
explicitly intended editable source to current `main` and verify it; do not use
an old-path restore operation.

## Independent radia-optuna lane

```powershell
python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all
python tools/release_quad.py optuna-done --wheel <path>
```

The candidate is the exact successful-main wheel. Bind its SHA256, source SHA,
version, CI run, and all declared target results before tagging and publishing.
Do not run solver Phase 8 or install other distributions.

## Radia exact-artifact publication hold

Radia tag builds do not automatically publish to PyPI. Accept the exact tag-CI
wheel on the required solver validation hosts, commit the required machine-
readable evidence, and dispatch the release workflow with the exact CI run,
wheel SHA256, and acceptance commit. Promote the verified artifact without
rebuilding. Changed bytes require new acceptance.

This does not replace QUAD `done` or the Simulink candidate gate. A full
Simulink package includes `radia_simulink_library.slx`, support files, standalone
MEX handles, runtime DLLs, `manifest.json`, and `SHA256SUMS.txt`. Verify the exact
ZIP independently on LAB, 100, mdx1, and mdx2 through verified MATLAB Engine
sessions. Rebuilding the ZIP invalidates all prior candidate evidence.

Check MATLAB processes, shared Engine names, and the official MCP connection
separately before acceptance. Reuse an appropriate existing session explicitly:
`simulink-candidate --package <zip> --target 100 --engine-session 100=<name>`.
The verifier checks its PID, refuses loaded Simulink diagrams or existing
Radia/Optuna MEX handles, restores the borrowed path and environment, and does
not quit it. Without an explicit session it starts MATLAB only when no MATLAB
process or shared Engine exists; inaccessible existing sessions are not a
reason to launch a substitute. Close only sessions owned by this operation.

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
