# Radia MCP Architecture & Naming Policy

## 1. Two-Track Management Philosophy
- **Deterministic Solvers (`src/radia`, C++/CUDA, AMS, PEEC-fin)**: Strict Shogi-style engineering — clean worktrees, CI, exact-SHA validation, reciprocal review.
- **MCP Tooling (`packages/radia-mcp`)**: Agile Go-style engineering — edit live editable source, run scoped verification, keep the tool surface thin and easy to rewrite. Do not add worktree/ceremony overhead for routine MCP maintenance.

## 2. Three-Layer Naming Convention
1. **Standalone External Software MCPs**: `{software}-mcp`
   - Examples: `build123d-mcp`, `onshape-mcp`, `gmsh-mcp`, `blender-mcp`, `cubit-mcp`, `optuna-mcp`, `matlab-mcp`, `jmag-mcp`, `femtet-mcp`, `comsol-mcp`, `cst-mcp`, `opera-mcp`, `femm-mcp`, `ltspice-mcp`, `emdlab-mcp`, `elf-magic-mcp`, `solidworks-mcp`.
2. **Software Knowledge Inside `radia-mcp`**: `radia-mcp.{software}`
   - Examples: `radia-mcp.mathematica`, `radia-mcp.ngsolve`, `radia-mcp.matlab`.
   - Avoids collision with official/external `{software}-mcp` servers while clarifying Radia-side integration knowledge.
3. **Radia Owned Technologies & Domains**: `radia-{technology}`
   - Examples: `radia-ih`, `radia-mor`, `radia-kelvin-transformation`, `radia-motor`, `radia-maglev`, `radia-electromagnet`, `radia-pcb`, `radia-matrix-solvers`.

## 3. Core Server Design (`mcp-core-server` Pattern)
Follow MathWorks' `matlab-mcp-core-server` design principles for every MCP server:
- **Minimal, orthogonal core tools** (target ~5–10 core `@mcp.tool`s per server: execute, inspect, validate, status).
- **Lazy initialization**: fast startup; load heavy engines/solvers only on tool invocation.
- **Separation of execution and knowledge/gates**: keep `server.py` lean; move validation-catalog wrappers (`_validation_tools.py`) and domain knowledge out of `server.py`.

## 4. On-Demand Commercial MCP Loading
- **Global always-on**: Core Radia & open CAE servers only (fast startup, low context overhead).
- **On-demand**: Commercial solver and document MCPs (`commercial-mcp-servers.toml` / plugins / per-project `.mcp.json`) loaded only when needed.
