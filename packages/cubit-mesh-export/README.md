# cubit-mesh-export

High-order curved mesh export from [Coreform Cubit](https://coreform.com/products/coreform-cubit/) to [NGSolve](https://ngsolve.org)/[Netgen](https://github.com/NGSolve/netgen).

## Features

- **Cubit plugin**: `export netgen/gmsh/nastran/vtk` APREPRO commands + Export Mesh GUI
- **Arbitrary-order curving** (order 1-5) via ACIS geometry projection
- **Label preservation**: material (block), boundary (sideset), edge (BBND)
- **Companion JSON**: CAD reference values for Volume/Area/Length consistency checking
- **Standalone checker**: verify mesh quality without Cubit (`check-vol` CLI)

## Install

```bash
pip install cubit-mesh-export
cubit-plugin-install
```

The second command deploys the plugin binaries (.ccm, .ccl, .pyd) to your Coreform Cubit installation.

### Upgrade

```bash
pip install --upgrade cubit-mesh-export
cubit-plugin-install
```

Always run `cubit-plugin-install` after upgrading to deploy the new binaries.

## Usage

### Cubit commands (recommended)

After installation, these commands are available in Cubit:

```
export netgen "model.vol" order 3 overwrite     # NGSolve FEM (.vol)
export gmsh "model.msh" order 2 version 2 overwrite  # GMSH v2.2
export gmsh "model.msh" order 2 version 4 overwrite  # GMSH v4.1
export radia_nastran "model.bdf" order 2 overwrite    # Nastran BDF
export vtk "model.vtk" order 2 overwrite              # VTK Legacy
```

Or use the **Export Mesh** menu in the Cubit GUI.

### Python API (for scripting / verification)

```python
import netgen       # must import before cubit
import ngsolve
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd('open "model.cub5"')

from cubit_mesh_export import extract_curved_mesh
ng_mesh = extract_curved_mesh(cubit, order=3)
ng_mesh.Save("model.vol")
```

### Check (does NOT require Cubit)

```bash
check-vol model.vol                         # basic check
check-vol model.vol --json model.vol.json   # compare vs CAD values
```

```python
from cubit_mesh_export.check import check_consistency
results = check_consistency("model.vol")
```

## Part of the Radia project

Source: [github.com/ksugahar/Radia](https://github.com/ksugahar/Radia)
