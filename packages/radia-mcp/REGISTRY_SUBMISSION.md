# MCP registry submission metadata for `radia-mcp`

This file describes the current `radia-mcp` distribution. Registry listings
must be refreshed from this file and the live `radia-mcp-meta` catalog;
do not copy historical fixed server counts or retired entry points.

## One-line pitch

> MCP guidance and executable workflows for electromagnetic CAE, NGSolve,
> build123d, Gmsh, optimization, scientific writing, and the Radia ecosystem.

## Name

`radia-mcp`

## Description

MCP servers for electromagnetic CAE and research workflows: Radia/NGSolve
analysis, mixed-omega formulations, open-boundary methods, application packs,
optimization, CAD authoring with build123d, Gmsh inspection and validation.
Cubit MCP is
distributed separately by `cubit-mesh-export` and can be installed as an
optional integration.

## Keywords / tags

```
electromagnetics, cae, fem, bem, ngsolve, radia, optimization,
scientific-computing, model-context-protocol, build123d, gmsh,
motion-coupling, mixed-omega
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
| License | Apache-2.0 |

## Discovery and representative entry points

The wheel contains a catalog of focused MCP servers. The catalog is discovered
at runtime rather than documented as a fixed count.

```text
radia-mcp-meta           # authoritative catalog and cross-server routing
radia-mcp-analysis       # field formulations, numerical methods, verification
radia-mcp-optimize       # optimization, topology and stream-function coil design
radia-mcp-learning       # surrogate models, PINN and data assimilation
radia-mcp-ih             # induction heating
radia-mcp-maglev         # magnetic levitation and motion
radia-mcp-electromagnet  # electromagnets and accelerator magnets
radia-mcp-motor          # electric machines
radia-mcp-pcb            # PCB coils, wireless power, PEEC
radia-mcp-acoustic       # acoustic FEM/BEM routes
radia-mcp-rna            # reluctance networks and hysteresis
build123d-mcp            # Python-native parametric CAD
gmsh-mcp                 # Gmsh inspection and post-processing
```

Every server runs as `radia-mcp serve <name>`. Start with `radia-mcp-meta` when
selecting a server; the live catalog is returned by `radia_mcp_overview()`.

## MCP client configuration

```json
{
  "mcpServers": {
    "radia-mcp-meta": {"command": "radia-mcp", "args": ["serve", "radia-mcp-meta"]},
    "radia-mcp-analysis": {"command": "radia-mcp", "args": ["serve", "radia-mcp-analysis"]},
    "radia-mcp-optimize": {"command": "radia-mcp", "args": ["serve", "radia-mcp-optimize"]}
  }
}
```

## Optional Cubit integration

Cubit ownership is deliberately outside this distribution:

```bash
pip install cubit-mesh-export
```

`cubit-mesh-export` owns `cubit_mesh_export.mcp`, its Cubit MCP entry point, plugin binaries, Cubit-specific tests, and Cubit release CI.
`radia-mcp` may call that independently installed package for optional
STEP-to-Cubit workflows, but its own wheel must not contain the retired
`radia_mcp.cubit` namespace or a Cubit entry point.

## Maintainer contact

Kengo Sugahara (ksugahar@ele.kindai.ac.jp)
Kindai University, Sugahara Laboratory


Owner-authored code is licensed under Apache-2.0 from the next releases.
See LICENSE, NOTICE and THIRD_PARTY_NOTICES.md. Third-party code and assets
retain their original licenses; previously published artifacts are unchanged.
