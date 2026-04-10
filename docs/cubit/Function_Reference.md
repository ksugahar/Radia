# Function Reference

Reference documentation for the Radia Cubit plugin APREPRO commands and Python API.

## APREPRO Commands (Recommended)

Native Cubit commands registered by the Radia plugin. No Python import needed.
All commands are available in journal files (.jou) and the Cubit command line.

### Mesh Export Commands

| Command | Format | Orders | Block Required |
|---------|--------|--------|---------------|
| `radia_export netgen "f.vol" order N` | Netgen .vol (+ .vol.json) | 1-5 | No |
| `radia_export gmsh "f.msh" order N version 2` | Gmsh v2.2 | 1-4 | No |
| `radia_export gmsh "f.msh" order N version 4` | Gmsh v4.1 | 1-4 | No |
| `radia_export nastran "f.bdf" order N` | Nastran BDF | 1-2 | No |
| `radia_export vtk "f.vtk" order N` | VTK Legacy | 1-2 | No |

### Coil Generation Command

| Command | Description |
|---------|-------------|
| `coil "script.py"` | Generate coil STEP from CoilBuilder script + import |
| `coil "script.py" output "path.step"` | Custom output path |
| `coil "script.py" noimport` | Generate STEP without importing |

> **IMPORTANT**: Use `radia_export`, NOT `export`.
> Cubit has built-in `export nastran` and `export abaqus` commands with different
> formats and no high-order support. `export gmsh` does NOT exist in Cubit — only
> `radia_export gmsh` is available. Using `export gmsh` will fail with
> "is not a valid type of file to be exported".

### Build & Installation

The plugin is built with **compact_netgen** (static link, no nglib.dll dependency):

```bash
# Recommended: compact_netgen (no ABI mismatch risk)
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCubit_DIR="C:/Program Files/Coreform Cubit 2025.3/cmake" \
  -DNETGEN_SRC_DIR="C:/netgen_build/netgen_fork" \
  src/cubit_plugin

cmake --build . --target radia_cubit_ccm   # APREPRO commands (plugins/)
cmake --build . --target radia_cubit_ccl   # Qt5 GUI menu (bin/)
```

Installation: `pip install cubit-mesh-export && cubit-plugin-install`

---

## Command Details

### radia_export netgen

```
radia_export netgen "filename.vol" [order <1-5>] [overwrite]
```

Exports mesh with high-order curving (CallbackGeometry + ACIS projection).
Produces companion JSON (.vol.json) with CAD reference values for consistency checks.

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | Curve order (1=linear, 2-5=high-order via NetgenCurver) |
| overwrite | off | Overwrite existing file |

### radia_export gmsh

```
radia_export gmsh "filename.msh" [order <1-4>] [version <2|4>] [dimension <2|3>] [overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | Element order (1-4, order 3+ requires NetgenCurver) |
| version | 2 | GMSH format (2=v2.2, 4=v4.1) |
| dimension | 3 | 2D or 3D mode |

### radia_export nastran

```
radia_export nastran "filename.bdf" [order <1|2>] [dimension <2|3>] [nopyramid] [overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | Element order (1=CTETRA/CHEXA, 2=CTETRA10/CHEXA20) |
| dimension | 3 | 2D (CTRIA3/CQUAD4) or 3D |
| nopyramid | off | Convert pyramids to degenerate hex (JMAG compatible) |

### radia_export vtk

```
radia_export vtk "filename.vtk" [order <1|2>] [dimension <2|3>] [overwrite]
```

VTK Legacy format. Cell types: TET(10), HEX(12), WEDGE(13), PYRAMID(14), TRI(5), QUAD(9).

### coil

```
coil "script.py" [output "path.step"] [noimport]
```

Generates coil STEP via external Python 3.12 subprocess (CoilBuilder).
The script must define `build_coil()` returning a `CoilBuilder` instance.

Requires: Python 3.12 with NGSolve/OCC. Set `RADIA_PYTHON` env var to override.

---

## Python API

### extract_curved_mesh (cubit-mesh-export package)

```python
from cubit_mesh_export import extract_curved_mesh
ng_mesh = extract_curved_mesh(cubit, order=3)
```

Returns `netgen.meshing.Mesh` with high-order curving. Requires Cubit running.

### check_consistency (cubit-mesh-export package)

```bash
check-vol model.vol                         # CLI
check-vol model.vol --json model.vol.json   # With companion JSON
```

```python
from cubit_mesh_export.check import check_consistency  # API
```

---

## GUI Menu Structure

```
Menu bar: ... Export Mesh  Help  Solve

Export Mesh (C++ .ccl):        Solve (Python):
  Netgen Vol (.vol)...           Radia-NGSolve...
  GMSH...                        Generate Coil...
  Nastran BDF...                 --------
  VTK...                         Reload Panels
  --------
  Mesh Evaluation...
```

- **Export Mesh**: Qt5 dialogs with settings persistence (`AppData/Roaming/Radia/export_settings.json`)
- **Solve**: Python subprocess to external Python 3.12 (Cubit embeds Python 3.10)
- **Generate Coil**: Calls `coil` APREPRO command via file dialog

---

## Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| "Interrupt Detected" on AddPoint | ABI mismatch (full Netgen DLL) | Rebuild with compact_netgen |
| HEX20 has 8 nodes in .msh | edge_ho_nodes_ bug (fixed 2026-04-05) | Update ccm |
| `export nastran` wrong format | Using Cubit built-in | Use `radia_export nastran` |
| cp932 UnicodeDecodeError | Non-ASCII in .py | Use ASCII only + encoding='utf-8' |
| ccm size < 400 KB | Old full-Netgen build | Rebuild with compact_netgen (~600 KB) |
