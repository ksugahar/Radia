# Shared MCP runtime policy

MCP is experimental development tooling. This policy supersedes the former
snapshot-freeze and per-edit deployment-approval requirements for MCP packages.
It does not relax numerical solver, native ABI or release acceptance.

## Development

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

## Reload And Reconnect

Editing source or running pip does not refresh existing Python objects or tool
schemas. Reload compatible code or reconnect affected clients through supported
controls as needed; use the mcp-reconnect skill. Preserve active jobs, and defer
disruptive operations on busy clients. Never mass-kill Python, MATLAB, Cubit or
client processes merely to refresh MCP.

Before deleting a superseded source tree, check that no active consumer still
needs it. Do not automatically restore an older tree because its path was once
canonical. Report source integration, installation and live-client verification
separately, with host/client and observation time. When not checked, say so.

## Scope

These permissions apply to MCP development, including editable repointing.
Numerical solvers and native binaries still require their independent tests,
provenance and release checks. Do not turn MCP experimentation into an excuse
to waive those checks or overwrite another task's changes.
