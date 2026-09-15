# Radia MCP maintenance contract

MCP is experimental development tooling. Developers may edit source and repoint
editable installations without a separate deployment approval or mandatory
snapshot. Coordinate overlapping work and preserve active CAD/MATLAB jobs.
Follow the [shared runtime policy](operations/mcp-runtime-policy.md); installation
and fresh-process checks do not establish what existing clients have loaded.

## Standard client configuration

The maintained baseline is `radia-meta`, `radia-build123d`,
`radia-gmsh`, `radia-analysis`, `radia-motion`, and `radia-publication`.
Additional servers are retained. Gmsh is for post-processing, not solver meshing.
Cubit MCP is owned by `cubit-mesh-export`, not this baseline. Its owner installs
and configures `python -s -m cubit_mesh_export.mcp.server`. Existing external
entries and their permissions are preserved. Retired `radia_mcp.cubit.server`
or `maintenance serve cubit` launchers cause a read-only conflict report with
the replacement module; migrate them through the owning package before retrying.
Do not delete a client's Cubit entry or silently re-enable it during migration.
The baseline includes tools that can execute code; existing client approval
policies still apply. Installation does not grant permission for arbitrary runs.

Install the small optional TOML editor with `radia-mcp[maintenance]`.
Run with the intended absolute Python executable and `-s` to isolate user-site
packages. Do not rely on whichever `python` an unrelated shell resolves.

```powershell
python -s -m radia_mcp.maintenance config C:\Users\NAME\.codex\config.toml
python -s -m radia_mcp.maintenance config C:\Users\NAME\.claude.json
# Review the JSON plan, stop concurrent configuration editors, then add --apply.
```

The command plans by default. `--apply` backs up existing bytes alongside the
file, preserves its existing ACL by writing in place, checks intervening edits,
and verifies the written bytes. A crash during writing is recoverable from the
backup; this is not an atomic transaction against other client writers. Protect
the parent directory and backup with the same access restrictions as the client
configuration. No secrets or environment values are included in plan reports.
TOML comments, nested tool policies, disabled entries, environment and working
directory settings survive migration. JSON formatting may change once.
Malformed input, custom wrappers, remote transports, and non-default launch
arguments on baseline names block the entire apply. Review them explicitly;
never silently expand a restricted profile or remove its approval policy.
Different aliases and project-scoped overrides need operator review; the command
does not delete or merge them. The same command is a byte-preserving no-op on
its second successful invocation.

Standard launchers call `python -s -m radia_mcp.maintenance serve <catalog-key>`.
The catalog restricts module selection. Arguments following the key go to the
server, e.g. `serve radia-motion --profile ih`. The launcher calls `main()`
directly, including servers without a module execution guard. It does not print
diagnostic JSON on the MCP stdout stream. SDK/server lifecycle owns shutdown;
maintenance never kills unrelated MATLAB, Python, or client processes.

## Repeatable editable update

radia-mcp uses **release-dual: LAB and 100 only**, independently of Radia's
solver release. Both keep editable installs. Do not deploy to hibino/mdx1/mdx2
or invoke the full QUAD installer for this package-only update. Compute-runner
CI remains isolated; existing excluded-host installations are left untouched.

The default is **update source -> reconnect -> check one affected live tool**.
For release-dual, LAB supplies the live-client acceptance. Verify the editable
installation and fresh import on both hosts; existing 100 clients can update
at their next normal restart without blocking release completion. Do not force
all-user restarts merely to close the checklist. See the
[release completion contract](operations/mcp-runtime-policy.md#release-completion).

1. Update the usual editable development checkout with reviewed changes. Check
   its actual path and coordinate overlapping edits; preserve others' WIP.
   Do not blindly pull/reset a dirty tree or switch to every temporary worktree.
   Pure Python edits at the same root normally need no reinstall.
2. Reconnect the affected client using supported controls and the mcp-reconnect
   skill, at a boundary without active work. Known-safe compatible hot reload
   is an optional shortcut. If automatic reconnect is unavailable, request one
   targeted manual Restart; do not kill Python, CAD, MATLAB or whole clients.
3. Check the live source and one harmless affected tool through that original
   client; check discovery too when tool names/schemas changed. Report the
   observed source and result briefly. Unknown evidence stays unverified;
   a Git merge or pip success alone is not live update completion.

No new daemon, watcher, generation database or mandatory receipt is needed.
The existing tools below are exceptions and diagnostics, not compulsory steps
for every edit. Default to the affected client/user, not every host/account.

### When more is needed

- **Source relocation or dependencies/package metadata changed:** use the
  intended Python to run
  `python -s -m pip install -e '<checkout>/packages/radia-mcp[maintenance]'`.
  Do not uninstall first. Check the exit code and actual import path, then
  reconnect. Locked entry points mean a pending update, not success.
- **Installation/launch changes or uncertain source:** run
  `python -s -m radia_mcp.maintenance doctor --expected-root
  <checkout>/packages/radia-mcp/src/radia_mcp --expected-version
  <selected-version> --expected-commit <actual-full-SHA>`.
  Record relevant uncommitted changes too. For launch/transport checks run
  `python -s tools/smoke_mcp_stdio.py --server <catalog-key>` from the package
  directory. Neither new process verifies an existing client; follow with the
  live check above. Investigate failures rather than bypassing them.
- **Client configuration changed:** plan/apply only the affected JSON/TOML
  settings, preserving disabled servers and access policies. Configuration
  migration is not required for ordinary source edits.
- **Immediate multi-user reconnection explicitly requested:** check each named user's launch and
  live client. Administrator success is not evidence for another user. LAB and
  100 have separate path namespaces; on 100 use its local `W:` paths, not LAB's
  `S:` paths or UNC. Do not copy credentials or alter unrelated site packages.
- **Busy CAD/MATLAB or native changes:** defer affected busy sessions. Native
  binaries require rebuilding, their numerical checks and a fresh loading
  process; editable installation does not rebuild or replace loaded binaries.

Fix failed MCP updates forward, never reinstall an old version as recovery.
Preserve others' work and report pending/manual/failed verification plainly.
Before removing an obsolete source, check consumers and unique work. Release
wheel tests remain isolated and retain independent acceptance gates.

## Behavioral acceptance, not just connectivity

Run existing tests, not duplicated maintenance versions. From the package root:

```powershell
python -s -m pytest tests/test_maintenance.py tests/test_capability_packs.py tests/test_paper_writing_review_regressions.py --junitxml=C:/temp/radia-mcp-maintenance.xml
```

Install maintenance/document extras as appropriate. A missing-dependency skip
is **not** evidence for that feature. Check JUnit counts. Cubit MCP tests and
licensed meshing acceptance belong to cubit-mesh-export, not this lane. For changed
Gmsh execution, run `tests/test_gmsh_post_guards.py` with the Gmsh extra. Real
solver/CAD/license acceptance stays in the corresponding validation lane.
Normal PR CI remains impact-scoped; this operator lane is not an all-server
test requirement for every documentation change.

## MathWorks implementation comparison

Inspected the official `matlab/matlab-mcp-server`, formerly MCP Core Server,
at commit `2e44b0ac789f43083cc29f5dce9079244b081a4c` (2026-09-07).

| Official source | Adopted Radia principle |
| --- | --- |
| [main.go](https://github.com/matlab/matlab-mcp-server/blob/2e44b0ac789f43083cc29f5dce9079244b081a4c/cmd/matlab-mcp-server/main.go) | Thin entry point; domain server owns execution. |
| [modeselector.go](https://github.com/matlab/matlab-mcp-server/blob/2e44b0ac789f43083cc29f5dce9079244b081a4c/internal/adaptors/application/modeselector/modeselector.go) | Separate serve, version and operator modes; validate before domain imports. |
| [config.go](https://github.com/matlab/matlab-mcp-server/blob/2e44b0ac789f43083cc29f5dce9079244b081a4c/internal/adaptors/application/config/config.go) | Parse then validate; redact sensitive configuration from diagnostics. |
| [lifecyclesignaler.go](https://github.com/matlab/matlab-mcp-server/blob/2e44b0ac789f43083cc29f5dce9079244b081a4c/internal/adaptors/application/lifecyclesignaler/lifecyclesignaler.go) | Explicit lifecycle ownership; no broad process killing. |

We do not copy the Go service graph, MATLAB watchdog, telemetry, session manager,
or fixed MATLAB shutdown timeout. Generic MATLAB operations remain with the
official server; Radia supplies only its domain layer and maintenance boundary.
Codex setting names and reconnect behavior follow the
[official MCP documentation](https://learn.chatgpt.com/docs/extend/mcp?surface=cli).

See [known issues](maintenance-known-issues.md) for remaining acceptance work.
