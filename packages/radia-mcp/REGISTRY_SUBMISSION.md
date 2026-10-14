# MCP registry submission metadata for `radia-mcp`

This file describes the current `radia-mcp` distribution. Registry listings
must be refreshed from this file and the live `mcp-server-radia-meta` catalog;
do not copy historical fixed server counts or retired entry points.

## One-line pitch

> MCP guidance and executable workflows for electromagnetic CAE, NGSolve,
> build123d, Gmsh, optimization, scientific writing, and the Radia ecosystem.

## Name

`radia-mcp`

## Description

MCP servers for electromagnetic CAE and research workflows: Radia/NGSolve
analysis, mixed-omega formulations, open-boundary methods, motion coupling,
optimization, CAD authoring with build123d, Gmsh inspection, validation,
publication figures, bibliography, and paper/grant preparation. Cubit MCP is
distributed separately by `cubit-mesh-export` and can be installed as an
optional integration.

## Keywords / tags

```
electromagnetics, cae, fem, bem, ngsolve, radia, optimization,
scientific-computing, model-context-protocol, build123d, gmsh,
motion-coupling, mixed-omega, technical-writing
```

## Categories

- Scientific computing / CAE / FEM / BEM
- Electromagnetic engineering
- Optimization and autonomous workflows
- CAD and mesh workflow integration
- Research communication

## URLs

| Field | Value |
|---|---|
| PyPI | https://pypi.org/project/radia-mcp/ |
| Source | https://github.com/ksugahar/Radia/tree/main/packages/radia-mcp |
| Install | `pip install radia-mcp` |
| License | BSD-3-Clause |

## Discovery and representative entry points

The wheel contains a catalog of focused MCP servers. The catalog is discovered
at runtime rather than documented as a fixed count.

```text
mcp-server-radia-meta       # authoritative catalog and cross-server routing
mcp-server-radia-ngsolve    # Radia/NGSolve methods and validated workflows
mcp-server-radia-design     # design and optimization capability pack
mcp-server-radia-motion     # motion-coupled electromagnetic workflows
mcp-server-build123d        # Python-native parametric CAD
mcp-server-gmsh             # Gmsh inspection and post-processing
mcp-server-paper-writing    # papers, grants, slides, and figures
```

Start with `mcp-server-radia-meta` when selecting a domain server. The complete
entry-point list is the `[project.scripts]` table in `pyproject.toml`; the live
catalog is returned by `radia_mcp_overview()`.

## MCP client configuration

```json
{
  "mcpServers": {
    "radia-meta": {
      "command": "mcp-server-radia-meta"
    },
    "radia-ngsolve": {
      "command": "mcp-server-radia-ngsolve"
    },
    "radia-design": {
      "command": "mcp-server-radia-design"
    }
  }
}
```

## Optional Cubit integration

Cubit ownership is deliberately outside this distribution:

```bash
pip install cubit-mesh-export
mcp-server-cubit
```

`cubit-mesh-export` owns `cubit_mesh_export.mcp`, the `mcp-server-cubit`
entry point, plugin binaries, Cubit-specific tests, and Cubit release CI.
`radia-mcp` may call that independently installed package for optional
STEP-to-Cubit workflows, but its own wheel must not contain the retired
`radia_mcp.cubit` namespace or a Cubit entry point.

## Maintainer contact

Kengo Sugahara (ksugahar@ele.kindai.ac.jp)
Kindai University, Sugahara Laboratory
