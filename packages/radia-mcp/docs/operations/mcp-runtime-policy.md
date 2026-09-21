# Shared MCP runtime policy

MCP is experimental development tooling. This policy supersedes the former
snapshot-freeze and per-edit deployment-approval requirements for MCP packages.
It does not relax numerical solver, native ABI or release acceptance.

## Development

### Lightweight editable updates (2026-09-15)

The routine is **update source -> reconnect -> check one affected live tool**.
Editable installation selects a directory; it does not advance Git, rebuild
native binaries, or refresh objects already loaded by a server.

- Keep the usual maintained development checkout as the editable target. Advance
  that source with reviewed changes while preserving others' WIP; do not repoint
  installations to every temporary review/build worktree or blindly pull/reset
  a dirty checkout. Intentional experiments may still select another source.
- Pure Python edits at the same root normally need no pip invocation. Run
  `pip install -e` with the intended interpreter only when relocating the source
  or changing dependencies/package metadata (including entry points). Do not
  uninstall first as a routine step.
- Reconnect only affected clients through supported controls. A known-safe
  compatible hot reload is an optional shortcut, not a prerequisite. If no
  supported automatic control exists, ask for one targeted manual Restart.
- Confirm the live source and one harmless affected tool through the original
  client. Include tool discovery when names or schemas changed. A short result
  in the task is enough; no new daemon, watcher, generation database, mandatory
  receipt file, or per-call Git/pip check is required.

Use existing doctor/stdio checks for installation or launch changes and for
unknown or contradictory evidence, not a full all-user audit on every edit.
Default to the affected client/user; expand only to explicitly requested targets.
Busy CAD/MATLAB work is deferred, never interrupted for a routine update. Native
changes still need rebuilding, numerical acceptance and a fresh process that
loads the binary. These exceptions do not make ordinary MCP edits a release gate.

### Forward-only updates (2026-09-15)

Always advance the maintained MCP source. Never revert to, reinstall or
redistribute an older MCP version as a recovery or deployment strategy.
Correct failed updates forward and repeat the relevant tests and live-client
checks. Until then, report the update as incomplete; do not present an old
runtime as the current deployment or bypass failed verification.

Old versions are not rollback reserves. Keep superseded files only when they
serve a specific debugging investigation, with source commit/version, purpose
and removal condition recorded. They must not become editable, fallback or
distribution targets. Prefer minimal logs, hashes and reproductions over whole
old installations. Routine updates do not require recoverable backup copies.

Remove obsolete copies after checking active consumers and reviewing unique
work. This does not authorize deleting others' changes or interrupting active
jobs. Git history remains historical evidence, not permission to deploy an old
version. Numerical solver/native release acceptance remains a separate contract.

Developers, including Codex and Claude Code, may directly edit live MCP source
and change its editable installation source with `pip install -e`. Routine
experiments require no dedicated branch, immutable snapshot or separate
deployment approval. Improve workflows through use and focused tests.

Coordinate overlapping source edits and changes to the same interpreter rather
than racing another developer. Preserve uncommitted work and active CAD/MATLAB
jobs. An MCP change does not require reinstalling unrelated solver packages.
Use an isolated environment for release wheel installation tests so those tests
do not accidentally replace the selected development installation.

## Verify What Changed

After repointing, check the selected interpreter and actual module import path.
Mapped-drive and UNC paths may identify the same source; compare resolved roots.
An intentional editable-source change is not automatically drift to repair.

Keep three observations separate:

- Installed registration: the source selected for an interpreter.
- Fresh-process import: what a newly started process actually resolves.
- Live client: the code and tool schema already held by that client's server.

Package-version equality and on-disk hashes do not establish loaded-code
identity. Use registration-time provenance when exposed and a harmless changed
tool call. Report contradictory evidence as unverified, with a mixed-generation
reason; one verified client does not establish all clients are current.

The intended source per host, interpreter and package is a recorded fact, not
an inference. `python tools/release_quad.py repoint` (and Phase 8 of a release)
writes it to `%ProgramData%\Radia\editable-intent.json`
(`RADIA_EDITABLE_INTENT_FILE`) together with the previous pointer, commit,
actor, time and reason, and appends every change to a log beside it.
`verify-editable` and `tools/verify_lab_editable.py` compare installations with
that record. A package without a record is UNVERIFIED, not drift, and no tool
proposes a repair target for it: record the current pointer if it is intended
(`repoint --record-current --reason ...`) or move it explicitly. `repoint` does
not uninstall first and does not stop processes. A pushed ref is required only
for formal handoff or completion evidence (`repoint --require-pushed`), not for
routine MCP development.

## Reload And Reconnect

Editing source or running pip does not refresh existing Python objects or tool
schemas. Reload compatible code or reconnect affected clients through supported
controls as needed; use the mcp-reconnect skill. Preserve active jobs, and defer
disruptive operations on busy clients. Never mass-kill Python, MATLAB, Cubit or
client processes merely to refresh MCP.

The `_reload_code` MCP tool is registered at server startup only when installed
distribution metadata confirms an editable installation. Set
`RADIA_MCP_HOT_RELOAD=0` in the server environment to disable its registration;
setting it to `1` does not enable it for a wheel or unknown installation.
After changing this setting, reconnect the affected client through supported
controls. This is a tool-exposure policy, not a restriction on source edits or
direct Python development helpers. Registration metadata alone still does not
prove the live imported source: keep the three provenance observations separate.

Reload refreshes or removes already registered tools only. A matching name prefix
does not authorize a new callable as a public tool. Reconnect to run normal
server registration and policy checks when adding tools.

The wheel verification lane runs every catalog server's stdio status contract
and `--selftest` from isolated Python outside the checkout. It requires actual
module paths under the installation root, noneditable registration, and no
exposed reload tool; a failure or timeout is not a successful skip.

Before deleting a superseded source tree, check that no active consumer still
needs it. Do not automatically restore an older tree because its path was once
canonical. Report source integration, installation and live-client verification
separately, with host/client and observation time. When not checked, say so.

## Release completion

For radia-mcp release-dual, completion requires passing the package release
checks, verified publication, verified editable updates/fresh imports on LAB
and 100, and LAB client reconnection with live source and harmless-tool checks.
The LAB check must cover affected servers/contracts, not an unrelated tool.

Existing clients on 100 may load the update at their next normal restart.
Report these as `next-launch-pending`; they do not block publication or release
completion and must not be described as live-verified. Failed installation or
fresh import on 100 still blocks deployment completion. Do not wait for every
student's live tool call or force a restart solely to close a release checklist.
Immediate all-user reconnection is a separate, explicitly requested operation.
This distinction does not weaken wheel, dependency, or numerical acceptance.

## Scope

radia-mcp is independently versioned and distributed. Its release-dual targets
are LAB and 100 only, using editable installations like cubit-mesh-export.
Do not deploy radia-mcp to hibino, mdx1 or mdx2, or run the full Radia QUAD
installer for an MCP-only update. Isolated CI/wheel tests on compute runners
are tests, not host deployment. Existing installations on excluded hosts are
not automatically removed. MCP updates do not wait for a Radia solver release.
Release-dual names this two-host scope, not a new release_dual.py command.

`cae-mcp-core` is retired and must not return as a dependency, import namespace,
wheel payload, or hidden shared release gate. radia-mcp owns its internal shared
runtime under `radia_mcp.common` and `radia_mcp._shared`; cubit-mesh-export owns
its independent equivalents and Cubit MCP. Cross-package handoffs use explicit
artifacts and public contracts rather than importing either distribution as a
private foundation. A future shared distribution requires an explicit ownership
and release decision, not restoration of `cae-mcp-core`.

These permissions apply to MCP development, including editable repointing.
Numerical solvers and native binaries still require their independent tests,
provenance and release checks. Do not turn MCP experimentation into an excuse
to waive those checks or overwrite another task's changes.
