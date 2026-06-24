# MCP registry submission metadata for `radia-mcp`

Copy-paste these fields into the submission forms:
- [glama.ai](https://glama.ai/mcp/servers) — auto-indexes from GitHub; submit via "Add a server"
- [pulsemcp.com/submit](https://www.pulsemcp.com/submit)
- [mcp.so/submit](https://mcp.so/submit)
- [modelcontextprotocol/servers](https://github.com/modelcontextprotocol/servers) — community listing, PR against README

---

## One-line pitch

> First-and-only public MCP server suite for Coreform Cubit, Gmsh,
> build123d, and the Radia CAE ecosystem. Ships `cubit_mesh_auto`
> (batch-validated scheme ladder + live GUI replay),
> `cubit_exec_safely` (auto-checkpoint dry-run), universal
> STEP-to-hex backend (`any_step_to_cubit_hex`) for any upstream CAD
> MCP, plus auto-scraped community knowledge (Coreform forum +
> gmsh/gmsh GitLab issues + StackOverflow).

## Name

`radia-mcp`

## Description (long, ≤ 512 chars)

MCP servers for the Radia CAE ecosystem: standalone Cubit + build123d
+ gmsh-post via Plan A (persistent Cubit GUI + PySide6 QTimer + file-
drop IPC). First-and-only public MCP for Gmsh and Cubit. Includes
cubit_mesh_auto (scheme-ladder batch validation → live GUI replay),
cubit_exec_safely (checkpoint + dry-run + apply), gmsh_post suite
(MSH v4.1 inspect/validate/convert + quality + post-view writer),
build123d with auto-generated API reference + 13 templates + lint,
CadQuery interop, and universal STEP→hex backend for any CAD MCP.

## Keywords / tags

```
mesh, meshing, cae, fem, simulation, cubit, coreform, gmsh,
build123d, cadquery, opencascade, occt, hex-mesh, tet-mesh,
post-processing, msh, step, netgen, ngsolve, radia, magnet,
electromagnetism, scientific-computing
```

## Categories

- CAD / 3D modeling
- Scientific computing / CAE / FEM
- Mesh generation
- File conversion

## URLs

| Field | Value |
|---|---|
| PyPI | https://pypi.org/project/radia-mcp/ |
| Source | https://github.com/ksugahar/Radia |
| Install | `pip install radia-mcp` (core) / `pip install radia-mcp[full]` (everything) |
| License | BSD-3-Clause |

## Entry points (3 MCP servers shipped in one wheel)

```
mcp-server-cubit         # Coreform Cubit (mesh generation, hex/tet, plugins)
mcp-server-build123d     # Python-native parametric CAD (OCCT)
mcp-server-gmsh-post     # Gmsh post-processing (v4.1 inspect/validate/convert,
                         # quality, boundary, $NodeData/$ElementData writers)
```

## Installation snippet (for "MCP client config" section)

```json
{
  "mcpServers": {
    "cubit": {
      "command": "mcp-server-cubit"
    },
    "build123d": {
      "command": "mcp-server-build123d"
    },
    "gmsh-post": {
      "command": "mcp-server-gmsh-post"
    }
  }
}
```

## Highlights (for a feature bullet list)

- **First-and-only** public MCP for Gmsh, Cubit, and mesh generators
  more broadly (no NETGEN/TetGen/Triangle/MeshLab/Pointwise MCPs
  exist anywhere else, as of 2026-04-20).
- **Batch-validated then committed to live GUI** — `cubit_mesh_auto`
  tries a scheme ladder headlessly, replays only the winning recipe
  in the live Cubit GUI (user sees success path only).
- **Safety layer** — `cubit_exec_safely` auto-checkpoints live state
  to `.cub5`, dry-runs in batch, applies only on success; silent-
  error detection via `cubit.get_error_count()`.
- **MSH v4.1 standardization** — `gmsh_post_convert` lifts any older
  .msh to v4.1 + spec-compliance validator + `$NodeData` /
  `$ElementData` / view-from-CSV writers.
- **Deep knowledge** — auto-generated API references (build123d:
  142 classes / 65 functions / 1673 lines; gmsh: 651 function
  entries / 2008 lines) + curated cookbooks (Plane/Axis/Location,
  Builder vs Algebra rosetta, $NodeData usage, physical groups) +
  live-scraped examples (Coreform forum, build123d GitHub issues +
  Discussions, gmsh/gmsh GitLab issues, SciComp.SE / SO).
- **Universal CAD-MCP mesh backend** — `any_step_to_cubit_hex` lets
  any upstream CAD MCP (FreeCAD, OpenSCAD, Blender, Onshape,
  AutoCAD, KiCad, CadQuery, …) hand off a STEP and get a hex mesh
  back via the same ladder.

## Screenshot / demo suggestion (optional)

A 30-second demo of `build123d_to_cubit_hex(template="helix_coil")`:
    script → STEP → batch ladder win on `scheme auto` → live Cubit
GUI shows 1668 hex / 0 tet / 3780 nodes in ~15s. Great for the
"hero visual" on registry pages.

## Maintainer contact

Kengo Sugahara (ksugahar@ele.kindai.ac.jp)
Kindai University, CAE research group
