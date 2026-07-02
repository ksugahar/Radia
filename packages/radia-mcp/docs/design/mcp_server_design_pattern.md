# MCP Server Design Pattern

This document adapts maintainability patterns from the public MathWorks MATLAB
MCP Server source to `radia-mcp`.  It is a design transfer, not source-code
vendoring.

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
3. Put calculations or parsing in a helper that can be unit-tested without MCP.
4. Make the tool return a compact dictionary or text block with stable keys.
5. Add at least one fast unit test.
6. Add or update public docs only with scrubbed, reproducible information.
7. Run policy lint and the narrowest relevant pytest.

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
