# radia-mcp

MCP (Model Context Protocol) servers for the Radia CAE ecosystem.

## Install

```bash
pip install radia-mcp
```

## Servers

| Server | Command | Purpose |
|--------|---------|---------|
| Radia | `mcp-server-radia` | Radia MMM/MSC/PEEC usage, linting |
| NGSolve | `mcp-server-ngsolve` | NGSolve FEM/BEM usage, linting |
| Cubit | `mcp-server-cubit` | Coreform Cubit meshing, linting |
| GMSH | `mcp-server-gmsh` | GMSH post-processing, linting |

## Claude Code Configuration

Add to your Claude Code MCP settings:

```json
{
  "mcpServers": {
    "radia": {"command": "mcp-server-radia"},
    "ngsolve": {"command": "mcp-server-ngsolve"},
    "cubit": {"command": "mcp-server-cubit"},
    "gmsh": {"command": "mcp-server-gmsh"}
  }
}
```

For local development (no pip install needed):

```json
{
  "mcpServers": {
    "radia": {"command": "python", "args": ["src/radia/mcp_server/radia/server.py"]}
  }
}
```
