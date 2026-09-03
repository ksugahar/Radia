# MCP Server Design Pattern

This document adapts maintainability patterns from the public MathWorks MATLAB
MCP Server source to `radia-mcp`.  It is a design transfer, not source-code
vendoring.

The reference points are the official server's
[annotation presets](https://github.com/matlab/matlab-mcp-server/blob/v0.12.0/pkg/tools/annotations.go),
[structured-output wrapper](https://github.com/matlab/matlab-mcp-server/blob/v0.12.0/internal/adaptors/mcp/tools/basetool/withstructuredcontent.go),
[logger factory](https://github.com/matlab/matlab-mcp-server/blob/v0.12.0/internal/adaptors/logger/factory.go),
[session manager](https://github.com/matlab/matlab-mcp-server/blob/v0.12.0/internal/adaptors/globalmatlab/sessionmanager/sessionmanager.go),
and [real workflow tests](https://github.com/matlab/matlab-mcp-server/blob/v0.12.0/tests/system/workflow_suite_test.go).

`radia-mcp` is a suite of many small MCP servers.  That makes the public tool
surface powerful, but it also means every new tool needs a predictable shape.

## Core Pattern

- Keep MCP handlers thin: validate inputs, call an importable helper, return a
  compact structured result.
- Keep solver, CAD, meshing, notebook, and file-system boundaries behind small
  adapters.
- Put reusable domain logic in package modules, not directly in tool functions.
- Prefer manifest or file-path boundaries for large inputs and outputs.
- Keep fast tests small and deterministic; put heavy validation in the
  validation lane.
- Do not require optional solver backends merely to import an MCP server.
- Add public-boundary lint before claiming a tool or knowledge card is
  release-ready.

## Fleet Runtime Contract

Every server calls `register_status_tool` after registering its domain tools.
That one composition point now applies these rules to the complete fleet:

- every tool has a human-readable title and all four MCP annotation hints;
- explicit annotations are preserved, recognizable operations receive one of
  four shared presets, and ambiguous operations default to destructive;
- every tool records its annotation source and structured/unstructured output
  mode in metadata, while `<server>_status` reports a live aggregate audit;
- status is a structured object with the loaded module path, SHA-256 at
  registration, current on-disk SHA-256 and a changed-since-registration flag,
  plus Python executable, MCP SDK version, distribution version/location, and
  editable-install source URL when present;
- all tool calls pass one idempotent JSONL logging choke point. Logs retain
  tool name, timing, outcome, and argument type/length metadata, never argument
  values. `RADIA_MCP_CALL_LOG=0` disables it fleet-wide; a server-specific
  `RADIA_MCP_<SERVER>_CALL_LOG` setting takes precedence;
- `tools/smoke_mcp_stdio.py` exercises initialize, tools/list, and the status
  tools/call over a real stdio connection. Change-scoped CI probes affected
  servers, common-runtime changes probe all servers, and a built wheel must
  prove that its loaded module comes from the clean install location.

This common layer is control-plane infrastructure only. Solver state remains
inside domain adapters (for example the Gmsh session boundary), and the
existing core/full tool profiles continue to prevent a single giant public
tool surface. A grouped operation catalog is read-only; its dynamic runner is
classified as a non-destructive write because the selected operation may run a
solver or create an artifact.

## Python Mapping

| Public MCP pattern | `radia-mcp` equivalent |
|---|---|
| One boring entry point | `[project.scripts] mcp-server-*` entries in `pyproject.toml` |
| Public tool definition layer | `<domain>/server.py` or explicit tool registration module |
| Internal adapters | `<domain>` helpers for CAD, mesh, solver, file, or notebook boundaries |
| Domain services | Importable helper modules and curated knowledge cards |
| Guides | `docs/`, `skills/`, and notebook-oriented usage material |
| Fake/system tests | Fast pytest fixtures plus separate validation artifacts |

## New Tool Checklist

Before adding a new MCP tool:

1. Define the public name, title, and one-sentence use case.
2. Decide whether the tool is read-only, idempotent, or potentially expensive.
   Add an explicit annotation when name/description inference would be
   ambiguous; conservative defaulting is a safety net, not the target state.
3. Put calculations or parsing in a helper that can be unit-tested without MCP.
4. Make the tool return a compact dictionary or text block with stable keys.
   Prefer a typed dictionary/Pydantic result so FastMCP advertises an output
   schema. Text remains appropriate for long human-readable knowledge cards.
5. Add at least one fast unit test.
6. Add or update public docs only with scrubbed, reproducible information.
7. Run policy lint and the narrowest relevant pytest.
8. Run `python tools/smoke_mcp_stdio.py --server <catalog-name>` for a changed
   server boundary.

## Adapter Checklist

When a tool touches an external backend:

- keep import-time behavior backend-free;
- provide a preflight/status function before execution;
- avoid starting long-running GUI or solver sessions from ordinary knowledge
  queries;
- record version, run date, elapsed time, and dominant timing stages for
  validation artifacts;
- keep public examples driven by analytic/open references or stored scrubbed
  fixtures.

## Why This Matters

The public MCP API becomes easier to review when each tool has the same shape:
thin handler, tested helper, explicit adapter boundary, and a small stable
result contract.  That is the part of the official MCP server construction
worth copying across `radia-mcp`.
