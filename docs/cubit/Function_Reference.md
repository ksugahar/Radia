# Function Reference

Reference documentation for Cubit Mesh Export APREPRO commands and Python API.

## APREPRO Commands (Recommended)

Native Cubit commands registered by the Cubit Mesh Export plugin. No Python import needed.
All commands are available in journal files (.jou) and the Cubit command line.

### Mesh Export Commands

| Command | Format | Orders | HO Method | Doc |
|---------|--------|--------|-----------|-----|
| `export netgen "f.vol" order N` | Netgen .vol (+ .vol.json) | 1-5 | NetgenCurver + ACIS | [export_NetgenMesh](export_NetgenMesh.md) |
| `export gmsh "f.msh" order N` | Gmsh v4.2 | 1-3 | NetgenCurver + ACIS | [export_Gmsh](export_Gmsh.md) |
| `export nastran_bdf "f.bdf" order N` | Nastran BDF | 1-2 | NetgenCurver + ACIS | [export_Nastran](export_Nastran.md) |
| `export vtk "f.vtk" order N` | VTK Legacy | 1-2 | NetgenCurver + ACIS | [export_vtk](export_vtk.md) |
| `export meg "f.meg"` | ELF/MAGIC MEG | 1 | — | [export_meg](export_meg.md) |
| `export femeem "dir"` | FEMEEM (Gifu Univ.) | 1 (tet only) | — | [export_femeem](export_femeem.md) |

> **IMPORTANT**: The plugin's mesh exporters are `export netgen / gmsh /
> vtk / femeem / meg / nastran_bdf`. Cubit
> has a built-in `export nastran` (different format, no high-order support),
> so the plugin uses the distinct `nastran_bdf` keyword. No alternate spelling
> is registered.

### Build & Installation

The plugin `radia_cubit_export.dll` is built from
`packages/radia-mcp/native/cubit_vol_export` against the fork Netgen:

```bash
python packages/radia-mcp/native/cubit_vol_export/build.py \
  --cubit-sdk "<Cubit>/cmake" --netgen-dir "<venv>/Lib/site-packages/netgen/cmake" --out <plugin_dir>
# The packaged Export Netgen menu (radia_mcp.cubit.gui) is the sole GUI entry point.
```

Installation: `python -m radia_mcp.cubit.install --plugin-dir <plugin_dir>`

---

## Command Details

### export netgen

```
export netgen "filename.vol" [order <1-5>] [overwrite]
```

Exports mesh with high-order curving (CallbackGeometry + ACIS projection).
Produces companion JSON (.vol.json) with CAD reference values for consistency checks.

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 2 | Curve order (1=linear, 2-5=high-order via NetgenCurver) |
| overwrite | off | Overwrite existing file |

### export gmsh

```
export gmsh "filename.msh" [order <1-3>] [dimension <2|3>] [overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | Element order (1-3; wedge limited to order 2) |
| dimension | 3 | 2D or 3D mode |
| overwrite | off | Overwrite existing file |

Order 4-5 not supported (use `export netgen`).

### export nastran_bdf

```
export nastran_bdf "filename.bdf" [order <1|2>] [dimension <2|3>] [nopyramid] [overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | Element order (1=CTETRA/CHEXA, 2=CTETRA10/CHEXA20) |
| dimension | 3 | 2D (CTRIA3/CQUAD4) or 3D |
| nopyramid | off | Convert pyramids to degenerate hex (for importers without CPYRAM support) |

Blocks become PSOLID/PSHELL properties, sidesets become collision-free PSHELL
properties, and nodesets become SET1 cards. MAT cards are deliberately omitted;
assign real physical materials in the receiving application.

### export vtk

```
export vtk "filename.vtk" [order <1|2>] [dimension <2|3>] [overwrite]
```

VTK Legacy ASCII format. Cell types: TET(10/24), HEX(12/25), WEDGE(13/26), PYRAMID(14/27), TRI(5/22), QUAD(9/23).

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | Element order (1 or 2) |
| dimension | 3 | 2D or 3D mode |

### export meg

```
export meg "filename.meg" [threed|twod|axisymmetric] [labels "1:MMB,2:MWL,..."] [overwrite]
```

ELF/MAGIC MEG format. Block names define ELF element type prefixes (MMB, MWL, MCO, etc.).
Pyramids exported as degenerate 8-node hex. Nodesets/sidesets named `SPACE` become MGR2 spatial nodes.

| Parameter | Default | Description |
|-----------|---------|-------------|
| threed | yes | 3D analysis (DIM=T) |
| twod | | 2D planar (DIM=K, z=0) |
| axisymmetric | | Axisymmetric (DIM=R, y=0) |
| labels | | Per-block prefix override (`blockID:PREFIX,...`) |

### export femeem

```
export femeem "dirname" [scale <value>] [overwrite]
```

FEMEEM format (Gifu Univ. 3D FEM). Tet-only, 1st order.
Creates directory with `in.dat`, `sin.dat.B`, `sina.dat`, and `d3`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| scale | 1.0 | Coordinate scale factor |

## Python API

### In-Memory Curving

Curved meshes come from the Cubit command:

```python
cubit.cmd('export netgen "model.vol" order 3 overwrite')
```

The plugin curves the mesh in-process after Cubit has supplied linear
mesh arrays and geometry callbacks.

### check_consistency (radia.vol_check)

```bash
check-vol model.vol                         # Full .vol gate; sidecar optional
check-vol model.vol --strict-labels         # Enforce canonical label names
check-vol model.vol --contract labels.json  # Application/mode label contract
check-vol model.vol --report-json run/vol_check.json
```

```python
from radia.vol_check import check_consistency, check_label_contract
```

The checker runs after export and before solver/Simulink initialization. It
always validates mesh loading, label relations, and the curved NGSolve mapping.
When `model.vol.json` exists it is auto-discovered for CAD volume/area/length
and mesh metadata comparison; `--json` makes an explicit sidecar mandatory.
Material constants remain part of the application DesignSpec/configuration and
are not inferred from `.vol` label strings.

---

## GUI Menu Structure

```
Menu bar: ... Export Netgen  Help

Export Netgen (PySide6):
  Export Netgen (.vol)...
```

- **Export Netgen**: PySide6 dialog with settings persistence (`%APPDATA%/radia-mcp/export_settings.json`)
- Solver orchestration is intentionally outside the Cubit plugin. Downstream
  applications consume the exported mesh through their own process boundary.
- **Mesh p-convergence demo**: a saved notebook under `docs/`; it is not an engineering design panel or Cubit menu action.

---

## Export Format Comparison

| Feature | netgen | gmsh | nastran | vtk | meg | femeem |
|---------|--------|------|---------|-----|-----|--------|
| Max order | 5 | 3 | 2 | 2 | 1 | 1 |
| HO method | NetgenCurver | NetgenCurver | NetgenCurver | NetgenCurver | — | — |
| Tet | yes | yes | yes | yes | yes | yes |
| Hex | yes | yes | yes | yes | yes | no |
| Wedge | yes | yes (o2) | yes | yes | yes | no |
| Pyramid | yes | yes | yes | yes | degen hex | no |
| Labels | block+sideset | block+sideset | block | block | block+sideset | block |
| Companion | .vol.json | .geo | — | — | — | d3 |

---

## Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| "Interrupt Detected" on AddPoint | ABI mismatch (full Netgen DLL) | Rebuild with compact_netgen |
| HEX20 has 8 nodes in .msh | edge_ho_nodes_ bug (fixed 2026-04-05) | Update ccm |
| `export nastran` wrong format | Using Cubit built-in | Use `export nastran_bdf` |
| cp932 UnicodeDecodeError | Non-ASCII in .py | Use ASCII only + encoding='utf-8' |
| ccm size < 400 KB | Old full-Netgen build | Rebuild with compact_netgen (~600 KB) |
| GMSH order 4-5 fails | Not supported in GMSH export | Use `export netgen` |
| Wedge order 3 in GMSH | GMSH limitation | Falls back to linear; use netgen |
