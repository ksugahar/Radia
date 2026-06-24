# Cubit Mesh Export Examples

Example scripts for the radia Cubit C++ plugin (`cubit_mesh_export.ccm`), organized by export format.

## Folder Structure

| Folder | Format | Description |
|--------|--------|-------------|
| [gmsh/](gmsh/) | Gmsh (.msh) | v4.1, order 1-3 |
| [nastran/](nastran/) | Nastran BDF (.bdf) | CTETRA/CHEXA/CTRIA, order 1-2 |
| [other_formats/](other_formats/) | Custom formats | FreeFEM, ANSYS CDB, Lukas 2D |

## Validation examples

| Example | Shows |
|---|---|
| [`validation_vol_surface_closure.py`](validation_vol_surface_closure.py) | Netgen `.vol` boundary orientation, vector-area closure, and boundary/tet volume agreement |
| [`validation_vol_fem_bem_topology.py`](validation_vol_fem_bem_topology.py) | First-order FEM/BEM topology view: closed surface edge manifold, Euler characteristic, compact scalar-BEM nodes, RWG-to-HCurl edge trace |
| [`validation_vol_tet_quality.py`](validation_vol_tet_quality.py) | Netgen `.vol` tetrahedron quality: edge ratio, inradius, circumradius, radius-ratio quality, and optional real Cubit export evaluation |

```powershell
python validation_vol_surface_closure.py
python validation_vol_fem_bem_topology.py
python validation_vol_tet_quality.py
```

Each folder contains:
- Example Python scripts (using `cubit.cmd('export ...')`)
- Pre-generated sample output files
- README.md with format-specific documentation

## Running Examples

```powershell
python <folder>/<script>.py
```

Requires Cubit installed with the radia plugin (`radia-setup`).

## Export Commands (C++ .ccm plugin)

```
export gmsh "file.msh" [order {1|2}] [dimension {2|3}] [overwrite]
export jmag_nastran "file.bdf" [order {1|2}] [dimension {2|3}] [nopyramid] [overwrite]
export vtk "file.vtk" [order {1|2}] [dimension {2|3}] [overwrite]
```

## See Also

- [src/cubit_plugin/](../../src/cubit_plugin/) - C++ plugin source
- [tests/cubit/](../../tests/cubit/) - Automated tests (42/42 PASS)
