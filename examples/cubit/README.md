# Examples

Example scripts for `cubit_mesh_export` module, organized by export format.

## Folder Structure

| Folder | Format | Description |
|--------|--------|-------------|
| [exodus/](exodus/) | Exodus II (.exo) | Cubit's native format, full feature support |
| [gmsh/](gmsh/) | Gmsh (.msh) | v2.2 and v4.1 formats |
| [meg/](meg/) | MEG (.meg) | ELF/MAGIC electromagnetic solver |
| [nastran/](nastran/) | Nastran BDF (.bdf) | Structural analysis |
| [netgen/](netgen/) | Netgen (.vol) | NGSolve/Netgen integration |
| [vtk/](vtk/) | VTK (.vtk/.vtu) | ParaView visualization |

Each folder contains:
- Example Python scripts
- Pre-generated sample output files
- README.md with format-specific documentation

## Running Examples

```bash
# Using Cubit's Python
"${CUBIT_PATH:-C:/Program Files/Coreform Cubit 2025.3/bin}/python3/python.exe" <folder>/<script>.py

# Using system Python (netgen examples require ngsolve)
python netgen/netgen_sphere_example.py
```

## See Also

- [docs/cubit/](../../docs/cubit/) - Full API documentation
- [cubit_mesh_export.py](../../src/radia/cubit_mesh_export.py) - Main module
