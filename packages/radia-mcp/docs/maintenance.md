# Radia MCP maintenance contract

MCP is experimental development tooling. Developers may edit source and repoint
editable installations without a separate deployment approval or mandatory
snapshot. Coordinate overlapping work and preserve active CAD/MATLAB jobs.
Follow the [shared runtime policy](operations/mcp-runtime-policy.md); installation
and fresh-process checks do not establish what existing clients have loaded.

## Standard client configuration

The maintained baseline is `radia-meta`, `radia-build123d`, `radia-cubit`,
`radia-gmsh`, `radia-analysis`, `radia-motion`, and `radia-publication`.
Additional servers are retained. Gmsh is for post-processing, not solver meshing.
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

1. Select the development source and intended interpreter. Routine experiments
   require no dedicated branch, immutable snapshot or separate approval. Record
   the source path and actual commit plus uncommitted changes when reporting
   results; a commit alone does not identify modified source.
2. Inventory all explicitly named human users, both clients and project-scoped
   settings. LAB and 100 have separate executable/path namespaces. On 100 use
   its local `W:` path, not LAB's mapped `S:` or a UNC path.
3. Coordinate changes to the same source or interpreter with other developers.
   Preserve their uncommitted changes and active jobs. Do not automatically
   restore a source merely because it used to be canonical.
4. Using the intended Python, run
   `python -s -m pip install -e '<checkout>/packages/radia-mcp[maintenance]'`.
   Check the exit code. Locked Windows entry points are a blocked update, not
   success: arrange client shutdown/retry; do not kill all Python or manually
   fabricate package metadata. Reload compatible code or reconnect the affected
   client when needed; source edits do not refresh loaded objects automatically.
5. Run `doctor --expected-root <checkout>/packages/radia-mcp/src/radia_mcp
   --expected-version <selected-version> --expected-commit <actual-full-SHA>` through
   `python -s -m radia_mcp.maintenance`. A nonzero result blocks acceptance.
   This checks the new process, not already-running sessions. Version alone
   cannot distinguish editable changes; record the commit and source path too.
6. Plan/apply each user's JSON and TOML configurations, record conflicts, and
   verify access under the intended user. An administrator's import test is not
   proof that another user's launch works. Do not change that user's unrelated
   Python site packages or copy credentials to enable impersonation.
7. Run the existing fresh-process real-transport probe
   `python -s tools/smoke_mcp_stdio.py --server <catalog-key>` from the package
   directory for each distinct launch configuration. It verifies initialize,
   tools/list, status, schema annotations and loaded-source provenance using
   the standard launcher. Separately reload/reconnect affected clients as needed
   using supported controls and the mcp-reconnect skill. Verify live discovery,
   loaded-source evidence and a harmless changed call through each original
   client. Neither doctor nor this new stdio process verifies an existing client.

Recovery: select the intended known-working source in coordination with affected
developers; do not reset another worktree or automatically restore an old path.
Restore only affected settings from backups when appropriate, preserving ACLs.
Repeat fresh-process and live-client checks. Before removing an old source,
verify that no active consumer needs it. Release wheel tests use isolated
environments and retain their independent publication/acceptance gates.

## Behavioral acceptance, not just connectivity

Run existing tests, not duplicated maintenance versions. From the package root:

```powershell
python -s -m pytest tests/test_maintenance.py tests/test_capability_packs.py tests/test_paper_writing_review_regressions.py tests/test_cubit_webcut_conformal_hex_gate.py --junitxml=C:/temp/radia-mcp-maintenance.xml
```

Install maintenance/document extras as appropriate. A missing-dependency skip
is **not** evidence for that feature. Check JUnit counts. The Cubit suite tests
the conformal gate contract; it is not a live licensed meshing run. For changed
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
