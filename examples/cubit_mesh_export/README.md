# Cubit Mesh Export Examples

Example scripts for the radia Cubit C++ plugin (`cubit_mesh_export.ccm`), organized by export format.

## Folder Structure

| Folder | Format | Description |
|--------|--------|-------------|
| [gmsh/](gmsh/) | Gmsh (.msh) | v4.1, order 1-3 |
| [nastran/](nastran/) | Nastran BDF (.bdf) | CTETRA/CHEXA/CTRIA, order 1-2 |
| [other_formats/](other_formats/) | Custom formats | FreeFEM, ANSYS CDB, Lukas 2D |

Each folder contains:
- Example Python scripts (using `cubit.cmd('export ...')`)
- Pre-generated sample output files
- README.md with format-specific documentation

## Running Examples

```bash
python <folder>/<script>.py
```

Requires Cubit installed with the radia plugin (`radia-setup`).

## Export Commands (C++ .ccm plugin)

```
radia_export gmsh "file.msh" [order {1|2}] [dimension {2|3}] [overwrite]
radia_export nastran "file.bdf" [order {1|2}] [dimension {2|3}] [nopyramid] [overwrite]
radia_export vtk "file.vtk" [order {1|2}] [dimension {2|3}] [overwrite]
```

## See Also

- [src/cubit_plugin/](../../src/cubit_plugin/) - C++ plugin source
- [tests/cubit/](../../tests/cubit/) - Automated tests (42/42 PASS)
